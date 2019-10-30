#include "io.h"
#include "molecule.h"
#include "geometry.h"
#include "rotate.h"
#include "bonds.h"
#include <fstream>
#include <iostream>
#include "spdlog/spdlog.h"
#include "aux/eigensupport.h"

namespace Faunus {

// global instances available throughout Faunus.
std::vector<MoleculeData> molecules;
std::vector<ReactionData> reactions;

MoleculeData::MoleculeData(const std::string &name, ParticleVector particles,
                           const BasePointerVector<Potential::BondData> &bonds)
    : name(name), bonds(bonds) {
    for (auto particle : particles) {
        atoms.push_back(particle.id);
    }
    if (!particles.empty()) {
        conformations.push_back(particles);
    }
}

int &MoleculeData::id() { return _id; }

const int &MoleculeData::id() const { return _id; }

ParticleVector MoleculeData::getRandomConformation(Geometry::GeometryBase &geo, ParticleVector otherparticles) {
    assert(inserter != nullptr);
    return (*inserter)(geo, otherparticles, *this);
}

void MoleculeData::loadConformation(const std::string &file, bool keep_positions, bool keep_charges) {
    ParticleVector v;
    if (loadStructure(file, v, false, keep_charges)) {
        if (keep_positions == false)
            Geometry::cm2origo(v.begin(), v.end()); // move to origo
        conformations.push_back(v);
        for (auto &p : v) // add atoms to atomlist
            atoms.push_back(p.id);
    }
    if (v.empty())
        throw std::runtime_error("Structure " + file + " not loaded. Filetype must be .aam/.pqr/.xyz");
}

void MoleculeData::setInserter(std::shared_ptr<MoleculeInserter> ins) { inserter = ins; }

MoleculeData::MoleculeData() { setInserter(std::make_shared<RandomInserter>()); }

void to_json(json &j, const MoleculeData &a) {
    j[a.name] = {
        {"id", a.id()},
        {"atomic", a.atomic},
        {"rigid", a.rigid},
        {"compressible", a.compressible},
        {"activity", a.activity / 1.0_molar},
    };
    j[a.name].update(a.json_cfg);

    if(a.inserter != nullptr) {
        a.inserter->to_json(j[a.name]);
    }
    j[a.name]["atoms"] = json::array();
    for (auto id : a.atoms)
        j[a.name]["atoms"].push_back(atoms.at(id).name);
    if (!a.bonds.empty()) {
        j[a.name]["bondlist"] = a.bonds;
    }
    if (!a.exclusions.empty()) {
        j[a.name]["exclusionlist"] = a.exclusions;
    }
}

void MoleculeData::createMolecularConformations(const json &val) {
    assert(val.is_object());

    std::string traj = val.value("traj", std::string());
    if (traj.empty())
        return;

    conformations.clear(); // remove all previous conformations
    // read tracjectory w. conformations from disk
    FormatPQR::load(traj, conformations.vec);
    if (not conformations.empty()) {
        // create atom list
        atoms.clear();
        atoms.reserve(conformations.vec.front().size());
        for (auto &p : conformations.vec.front()) // add atoms to atomlist
            atoms.push_back(p.id);

        // center mass center for each frame to origo assuming whole molecules
        if (val.value("trajcenter", false)) {
            faunus_logger->debug("Centering conformations in trajectory file {}", traj);
            for (auto &p : conformations.vec) // loop over conformations
                Geometry::cm2origo(p.begin(), p.end());
        }

        // set default uniform weight
        std::vector<float> w(conformations.size(), 1);
        conformations.setWeight(w);

        // look for weight file
        std::string weightfile = val.value("trajweight", std::string());
        if (not weightfile.empty()) {
            std::ifstream f(weightfile.c_str());
            if (f) {
                w.clear();
                w.reserve(conformations.size());
                double _val;
                while (f >> _val)
                    w.push_back(_val);
                if (w.size() == conformations.size())
                    conformations.setWeight(w);
                else
                    throw std::runtime_error("Number of weights does not match conformations.");
            } else
                throw std::runtime_error("Weight file " + weightfile + " not found.");
        }
    } else
        throw std::runtime_error("Trajectory " + traj + " not loaded or empty.");
} // done handling conformations

// ============ NeighboursGenerator ============

NeighboursGenerator::NeighboursGenerator(const BondVector &bonds) {
    createBondMap(bonds);
}

void NeighboursGenerator::createBondMap(const BondVector &bonds) {
    for (auto &bond : bonds.find<Potential::StretchData>()) {
        auto add_bond = [this, &bond](int i, int j) -> void {
            auto atom_emplaced = bond_map.emplace(bond->index[i], AtomList{}); // find or create empty vector
            atom_emplaced.first->second.push_back(bond->index[j]); // atom_emplaced.<iterator>-><vector>.push_back()
        };
        const int head_ndx = 0, tail_ndx = 1;
        add_bond(head_ndx, tail_ndx); // add the bond to the map in both head-tail
        add_bond(tail_ndx, head_ndx); // and tail-head directions
    }
}

void NeighboursGenerator::generatePaths(int bonded_distance) {
    if(paths.empty()) {
        // starting with a path of length 0 bonds which are just separated atoms
        paths.emplace_back();
        for (auto &bonded_atoms : bond_map) {
          paths.back().emplace(AtomList{bonded_atoms.first});
        }
    }
    // continue wherever the last path generation ended
    for(int distance = paths.size(); distance <= bonded_distance; ++distance) {
        paths.emplace_back();
        const auto &base_paths = paths.rbegin()[1]; // last but one: a set of paths of length (distance - 1)
        for (auto &base_path : base_paths) {
            // try to extend the path with every atom bonded to the tail
            for(auto appending : bond_map.at(base_path.back())) {
                // avoid loops: the appending atom shall not be anywhere in the path yet
                if(find(base_path.begin(), base_path.end(), appending) == base_path.end()) {
                    auto path(base_path);
                    path.push_back(appending);
                    paths.back().emplace(path);
                }
            }
        }
    }
}

void NeighboursGenerator::generatePairs(AtomPairList &pairs, int bond_distance) {
    generatePaths(bond_distance);
    std::set<std::pair<int, int>> pairs_set; // use set to prevent duplicities
    for (auto &paths_set : paths) {
        for (auto &path : paths_set) {
            auto excluded_pair = std::make_pair(path.front(), path.back());
            if (excluded_pair.first != excluded_pair.second) {
                if (excluded_pair.first > excluded_pair.second) {
                    std::swap(excluded_pair.first, excluded_pair.second);
                }
                pairs_set.emplace(excluded_pair);
            }
        }
    }
    // copy the set to the output vector
    std::copy(pairs_set.begin(), pairs_set.end(), std::back_inserter(pairs));
}

// ============ MoleculeBuilder ============

void MoleculeBuilder::from_json(const json &j, MoleculeData &molecule) {
    if (is_used) {
        throw std::logic_error("MoleculeBuilder cannot be used twice");
    }
    is_used = true;
    try {
        if (j.is_object() == false || j.size() != 1) {
            throw ConfigurationError("invalid json");
        }
        auto j_molecule_it = j.items().begin(); // a persistent copy of iterator needed in clang
        molecule.name = molecule_name = (*j_molecule_it).key();
        auto &j_properties = (*j_molecule_it).value();
        molecule.id() = j_properties.value("id", molecule.id());
        molecule.atomic = j_properties.value("atomic", molecule.atomic);
        molecule.rigid = j_properties.value("rigid", molecule.rigid);
        molecule.compressible = j_properties.value("compressible", molecule.compressible);
        molecule.activity = j_properties.value("activity", molecule.activity / 1.0_molar) * 1.0_molar;

        readCompoundValues(j_properties);
        for (auto particle : particles) {
            molecule.atoms.push_back(particle.id);
        }
        if (!particles.empty()) {
            molecule.conformations.push_back(particles);
        }
        molecule.createMolecularConformations(j_properties); // fixme do better
        molecule.setInserter(createInserter(j_properties));
        molecule.bonds = bonds;
        molecule.exclusions = ExclusionsVicinity::create(particles.size(), exclusion_pairs);

        // todo better if these values have to be stored at all
        try {
            auto structure = j_properties.at("structure");
            if (structure.is_string()) {
                molecule.json_cfg["structure"] = structure;
            }
        } catch (json::out_of_range &) {}
        molecule.json_cfg["keepcharges"] = j_properties.value("keepcharges", true);

        // at this stage all given keys should have been accessed or "spend". If any are
        // left, an exception will be thrown.
        //        if (! val.empty()) {
        //            throw std::runtime_error("unused key(s):\n"s + val.dump() + usageTip["moleculelist"]);
        //        }
    } catch (std::exception &e) {
        throw std::runtime_error("JSON->molecule " + molecule_name + ": " + e.what());
    }
}

std::shared_ptr<MoleculeInserter> MoleculeBuilder::createInserter(const json &j_properties) {
    auto inserter_ptr = std::make_shared<RandomInserter>();
    inserter_ptr->from_json(j_properties);
    return inserter_ptr;
}

void MoleculeBuilder::readCompoundValues(const json &j) {
    auto is_atomic = j.value("atomic", false);
    if(is_atomic) {
        readAtomic(j);
    } else {
        readParticles(j);
        readBonds(j);
        if(isFasta(j)) {
            readFastaBonds(j);
        }
        readExclusions(j);
    }
}

void MoleculeBuilder::readAtomic(const json &j_properties) {
    auto j_atoms = j_properties.value("atoms", json::array());
    particles.reserve(j_atoms.size());
    for (auto atom_id : j_atoms) {
        auto atom_it = findName(atoms, atom_id.get<std::string>());
        if(atom_it == atoms.end()) {
            faunus_logger->error("Unknown atom '{}' in molecule '{}'", atom_id, molecule_name);
            throw ConfigurationError("unknown atom in atomic molecule");
        }
        particles.emplace_back(*atom_it);
    }
}

void MoleculeBuilder::readParticles(const json &j_properties) {
    auto j_structure_it = j_properties.find("structure");
    if (j_structure_it != j_properties.end()) {
        bool read_charges = j_properties.value("keepcharges", true);
        MoleculeStructureReader structure_reader(read_charges);
        structure_reader.readJson(particles, *j_structure_it);
    } else {
        // allow virtual molecules :-/
        // shall we rather try to fallback on readAtomic()?
        // throw ConfigurationError("structure of the molecule not given");
    }
}

void MoleculeBuilder::readBonds(const json &j_properties) {
    bonds = j_properties.value("bondlist", bonds);

    // assert that all bonds are *internal*
    for (auto &bond : bonds) {
        for (int i : bond->index) {
            if (i >= particles.size() || i < 0) {
                throw ConfigurationError("bonded atom index " + std::to_string(i) + " out of range");
            }
        }
    }
}

void MoleculeBuilder::readFastaBonds(const json &j_properties) {
    auto &j_structure = j_properties.at("structure");
    Potential::HarmonicBond bond; // harmonic bond
    bond.from_json(j_structure);  // read 'k' and 'req' from json
    for (int i = 1; i < particles.size(); ++i) {
        bond.index = {i - 1, i};
        bonds.push_back<Potential::HarmonicBond>(bond.clone());
    }
}

void MoleculeBuilder::readExclusions(const json &j_properties) {
    auto excluded_neighbours_distance = j_properties.value("excluded_neighbours", 0);
    if(excluded_neighbours_distance > 0) {
        NeighboursGenerator generator(bonds);
        generator.generatePairs(exclusion_pairs, excluded_neighbours_distance);
    }
    for (auto j_exclusion_pair : j_properties.value("exclusionlist", json::array())) {
        if (!j_exclusion_pair.is_array() || j_exclusion_pair.size() != 2) {
            throw ConfigurationError("unrecognized molecule's exclusion format");
        }
        exclusion_pairs.emplace_back(j_exclusion_pair[0].get<int>(), j_exclusion_pair[1].get<int>());
    }
}

bool MoleculeBuilder::isFasta(const json &j_properties) {
    auto j_structure_it = j_properties.find("structure");
    bool is_fasta = (j_structure_it != j_properties.end() && j_structure_it->find("fasta") != j_structure_it->end());
    return is_fasta;
}

// ============ MoleculeStructureReader ============

void MoleculeStructureReader::readJson(ParticleVector &particles, const json &j) {
    if (j.is_string()) {
        auto filename = j.get<std::string>();
        readFile(particles, filename);
    } else if (j.is_array()) {
        readArray(particles, j);
    } else if (j.is_object() && j.find("fasta") != j.end()) {
        readFasta(particles, j);
    } else {
        throw ConfigurationError("unrecognized structure format");
    }
}

void MoleculeStructureReader::readFile(ParticleVector &particles, const std::string &filename) {
    faunus_logger->info("Reading molecule configuration from file: {}", filename);
    auto success = Faunus::loadStructure(filename, particles, false, read_charges);
    if (!success) {
        throw ConfigurationError("unable to open structure file");
    }
}

void MoleculeStructureReader::readArray(ParticleVector &particles, const json &j_particles) {
    particles.reserve(j_particles.size());
    for (auto &j_particle_wrap : j_particles) {
        if (!j_particle_wrap.is_object() || j_particle_wrap.size() != 1) {
            throw ConfigurationError("unrecognized molecule's atom format");
        }
        auto j_particle_it = j_particle_wrap.items().begin(); // a persistent copy of iterator needed in clang
        auto atom_it = findName(atoms, (*j_particle_it).key());
        if (atom_it == atoms.end()) {
            faunus_logger->error("An unknown atom '{}' in the molecule.", (*j_particle_it).key());
            throw ConfigurationError("unknown atom in molecule");
        }
        Point pos = (*j_particle_it).value();
        particles.emplace_back(*atom_it, pos);
    }
}

void MoleculeStructureReader::readFasta(ParticleVector &particles, const json &j_fasta) {
    if (j_fasta.find("fasta") == j_fasta.end()) {
        throw ConfigurationError("invalid FASTA format");
    }
    std::string fasta = j_fasta.at("fasta").get<std::string>();
    Potential::HarmonicBond bond; // harmonic bond
    bond.from_json(j_fasta);      // read 'k' and 'req' from json
    particles = Faunus::fastaToParticles(fasta, bond.req);
}

void from_json(const json &j, MoleculeData &a) {
    MoleculeBuilder builder;
    builder.from_json(j, a);
}

void from_json(const json &j, std::vector<MoleculeData> &v) {
    v.reserve(v.size() + j.size());
    for (auto &i : j) {
        v.push_back(i);
        v.back().id() = v.size() - 1; // id always match vector index
    }
}

// ============ ExclusionsSimple ============

ExclusionsSimple ExclusionsSimple::create(int atoms_cnt, const std::vector<std::pair<int, int>> &pairs) {
    ExclusionsSimple exclusions(atoms_cnt);
    exclusions.add(pairs);
    return exclusions;
}

ExclusionsSimple::ExclusionsSimple(int size)
    : size(size), excluded_pairs(std::make_shared<std::vector<unsigned char>>()) {
    faunus_logger->log(size < 1000 ? spdlog::level::trace : spdlog::level::warn,
                       "creating exclusion matrix {}×{} for {} atoms", size, size, size);
    excluded_pairs->resize(size * size, false);
}

void ExclusionsSimple::add(const std::vector<std::pair<int, int>> &exclusions) {
    for(auto pair : exclusions) {
        add(pair.first, pair.second);
    }
}

inline void ExclusionsSimple::add(int i, int j) {
    if(i >= size || j >= size) {
        throw std::range_error("exclusion index out of range");
    }
//    if (i > j) {
//        std::swap(i, j);
//    }
    any_exclusions = true;
    excluded_pairs->at(i * size + j) = excluded_pairs->at(j * size + i) = true;
}

void from_json(const json &j, ExclusionsSimple &exclusions) {
    auto &exclusion_list = j;
    if (exclusion_list.is_array()) {
        for (auto exclusion_it = exclusion_list.begin(); exclusion_it != exclusion_list.end(); ++exclusion_it) {
            if (exclusion_it->is_array() && exclusion_it->size() == 2) {
                exclusions.add((*exclusion_it)[0], (*exclusion_it)[1]);
            } else {
                mcloop_logger->warn("Ignoring unknown exclusion: {}", exclusion_it->dump());
            }
        }
    } else {
        faunus_logger->warn("Unknown exclusions format: array expected.");
    }
}

void to_json(json &j, const ExclusionsSimple &exclusions) {
    j = json::array();
    for (auto m = 0; m < exclusions.size; ++m) {
        for (auto n = m + 1; n < exclusions.size; ++n) {
            if (exclusions.isExcluded(m, n)) {
                j.push_back({m, n});
            }
        }
    }
}

// ============ ExclusionsVicinity ============

ExclusionsVicinity ExclusionsVicinity::create(int atoms_cnt, const std::vector<std::pair<int, int>> &pairs) {
    auto distance = [](int i, int j) -> int { return j > i ? j - i : i - j; };
    int max_neighbours_distance = 0;
    for (auto pair : pairs) {
        auto neighbours_distance = distance(pair.first, pair.second);
        if (neighbours_distance > max_neighbours_distance) {
            max_neighbours_distance = neighbours_distance;
        }
    }
    ExclusionsVicinity exclusions(atoms_cnt, max_neighbours_distance);
    exclusions.add(pairs);
    return exclusions;
}

ExclusionsVicinity::ExclusionsVicinity(int atoms_cnt, int max_difference)
    : atoms_cnt(atoms_cnt), max_bond_distance(max_difference),
      excluded_pairs(std::make_shared<std::vector<unsigned char>>()) {
    faunus_logger->log(atoms_cnt * max_difference < 1'000'000 ? spdlog::level::trace : spdlog::level::warn,
                       "creating exclusion matrix {}×{} for {} atoms within distance {}", atoms_cnt, max_difference,
                       atoms_cnt, max_difference);

    excluded_pairs->resize(atoms_cnt * max_difference, false);
}

void ExclusionsVicinity::add(int i, int j) {
    if (i > j) {
        std::swap(i, j);
    }
    if (j >= atoms_cnt || j - i > max_bond_distance) {
        throw std::range_error("exclusion indexes out of range");
    }
    excluded_pairs->at(toIndex(i, j)) = true;
}

void ExclusionsVicinity::add(const std::vector<std::pair<int, int>> &pairs) {
    for (auto pair : pairs) {
        add(pair.first, pair.second);
    }
}

void to_json(json &j, const ExclusionsVicinity &exclusions) {
    j = json::array();
    for (auto it = exclusions.excluded_pairs->begin(); it != exclusions.excluded_pairs->end(); ++it) {
        if (*it) {
            int n = std::distance(exclusions.excluded_pairs->begin(), it);
            j.push_back(exclusions.fromIndex(n));
        }
    }
}

void from_json(const json &j, MoleculeInserter &inserter) { inserter.from_json(j); }
void to_json(json &j, const MoleculeInserter &inserter) { inserter.to_json(j); }

ParticleVector RandomInserter::operator()(Geometry::GeometryBase &geo, const ParticleVector &, MoleculeData &mol) {
    int cnt = 0;
    QuaternionRotate rot;
    bool containerOverlap; // true if container overlap detected

    if (std::fabs(geo.getVolume()) < 1e-20)
        throw std::runtime_error("geometry has zero volume");

    ParticleVector v = mol.conformations.get(); // get random, weighted conformation
    conformation_ndx = mol.conformations.index; // lastest index

    do {
        if (cnt++ > max_trials)
            throw std::runtime_error("Max. # of overlap checks reached upon insertion.");

        if (mol.atomic) {       // insert atomic species
            for (auto &i : v) { // for each atom type id
                if (rotate) {
                    rot.set(2 * pc::pi * random(), ranunit(random));
                    i.rotate(rot.first, rot.second);
                }
                geo.randompos(i.pos, random);
                i.pos = i.pos.cwiseProduct(dir) + offset;
                geo.boundary(i.pos);
            }
        } else {                  // insert molecule
            if (keep_positions) {        // keep original positions (no rotation/trans)
                for (auto &i : v) // ...but let's make sure it fits
                    if (geo.collision(i.pos))
                        throw std::runtime_error("Error: Inserted molecule does not fit in container");
            } else {
                Point cm;                                        // new mass center position
                geo.randompos(cm, random);                       // random point in container
                cm = cm.cwiseProduct(dir);                       // apply user defined directions (default: 1,1,1)
                Geometry::cm2origo(v.begin(), v.end());          // translate to origin
                rot.set(random() * 2 * pc::pi, ranunit(random)); // random rot around random vector
                if (rotate) {
                    Geometry::rotate(v.begin(), v.end(), rot.first);
                    assert(Geometry::massCenter(v.begin(), v.end()).norm() < 1e-6); // cm shouldn't move
                }
                for (auto &i : v) {
                    i.pos += cm + offset;
                    geo.boundary(i.pos);
                }
            }
        }

        if (v.empty())
            throw std::runtime_error("Nothing to load/insert for molecule '"s + mol.name + "'");

        // check if molecules / atoms fit inside simulation container
        containerOverlap = false;
        if (allow_overlap == false) {
            for (auto &i : v) {
                if (geo.collision(i.pos)) {
                    containerOverlap = true;
                    break;
                }
            }
        }
    } while (containerOverlap);
    return v;
}

void RandomInserter::from_json(const json &j) {
    dir = j.value("insdir", dir);
    offset = j.value("insoffset", offset);
    rotate = j.value("rotate", rotate);
    keep_positions = j.value("keeppos", keep_positions);
}

void RandomInserter::to_json(json &j) const {
    j["insdir"] = dir;
    j["insoffset"] = offset;
    j["rotate"] = rotate;
    j["keeppos"] = keep_positions;
}

bool Conformation::empty() const {
    return positions.empty() && charges.empty();
}

ParticleVector &Conformation::toParticleVector(ParticleVector &p) const {
    assert(not p.empty() and not empty());
    // copy positions
    if (positions.size() == p.size()) {
        for (size_t i = 0; i < p.size(); ++i) {
            p[i].pos = positions[i];
        }
    }
    // copy charges
    if (charges.size() == p.size()) {
        for (size_t i = 0; i < p.size(); ++i) {
            p[i].charge = charges[i];
        }
    }
    return p;
}

bool ReactionData::empty(bool forward) const {
    if (forward)
        if (canonic)
            if (N_reservoir <= 0)
                return true;
    return false;
}
std::vector<int> ReactionData::participatingMolecules() const {
    std::vector<int> v;
    v.reserve(_reagid_m.size() + _prodid_m.size());
    for (auto i : _reagid_m)
        v.push_back(i.first);
    for (auto i : _prodid_m)
        v.push_back(i.first);
    return v;
}

bool ReactionData::containsMolecule(int molid) const {
    return _reagid_m.count(molid) > 0 || _prodid_m.count(molid) > 0;
}

const ReactionData::Tmap &ReactionData::Molecules2Add(bool forward) const { return (forward) ? _prodid_m : _reagid_m; }
const ReactionData::Tmap &ReactionData::Atoms2Add(bool forward) const { return (forward) ? _prodid_a : _reagid_a; }

void from_json(const json &j, ReactionData &a) {

    if (j.is_object() == false || j.size() != 1)
        throw std::runtime_error("Invalid JSON data for ReactionData");

    for (auto &m : Faunus::molecules)
        for (auto &a : Faunus::atoms)
            if (m.name == a.name)
                throw std::runtime_error("Molecules and atoms must have different names");

    for (auto it = j.begin(); it != j.end(); ++it) {
        a.name = it.key();
        auto &val = it.value();
        a.canonic = val.value("canonic", false);
        a.neutral = val.value("neutral", false);
        if (val.count("lnK") == 1)
            a.lnK = val.at("lnK").get<double>();
        else if (val.count("pK") == 1)
            a.lnK = -std::log(10) * val.at("pK").get<double>();
        a.pK = -a.lnK / std::log(10);
        a.N_reservoir = val.value("N_reservoir", a.N_reservoir);

        // get pair of vectors containing reactant and product species
        auto process = parseProcess(a.name);
        a._reag = process.first;
        a._prod = process.second;

        for (auto &name : a._reag) {                // loop over reactants
            auto pair = a.findAtomOrMolecule(name); // {iterator to atom, iterator to mol.}
            if (pair.first != atoms.end()) {
                a.swap = true; // if the reaction involves atoms, identify it as swap move
            }
        }

        for (auto &name : a._reag) {                // loop over reactants
            auto pair = a.findAtomOrMolecule(name); // {iterator to atom, iterator to mol.}
            if (pair.first != atoms.end()) {
                if (pair.first->implicit) {
                    // implicit reagent? K is multiplied by its activity
                    a.lnK += std::log(pair.first->activity / 1.0_molar);
                } else {
                    a._reagid_a[pair.first->id()]++;
                }
            }
            if (pair.second != molecules.end()) {
                a._reagid_m[pair.second->id()]++;
                if (pair.second->activity > 0) {
                    // explicit reagent?
                    // its activity is not part of K?
                    // K is divided by its activity
                    a.lnK -= std::log(pair.second->activity / 1.0_molar);
                }
            }
        }

        for (auto &name : a._prod) { // loop over products
            auto pair = a.findAtomOrMolecule(name);
            if (pair.first != atoms.end()) {
                if (pair.first->implicit) {
                    // implicit product? K is divided by its activity
                    a.lnK -= std::log(pair.first->activity / 1.0_molar);
                } else {
                    a._prodid_a[pair.first->id()]++;
                }
            }
            if (pair.second != molecules.end()) {
                a._prodid_m[pair.second->id()]++;
                if (pair.second->activity > 0) {
                    // explicit product?
                    // its activity is not part of K?
                    // K is multiplied by its activity
                    a.lnK += std::log(pair.second->activity / 1.0_molar);
                }
            }
        }
    }
}

void to_json(json &j, const ReactionData &a) {
    j[a.name] = {{"pK", a.pK},
                 {"pK'", -a.lnK / std::log(10)},
                 //{"canonic", a.canonic }, {"N_reservoir", a.N_reservoir },
                 {"neutral", a.neutral},
                 {"products", a._prod},
                 {"reactants", a._reag}};
} //!< Serialize to JSON object

} // namespace Faunus
