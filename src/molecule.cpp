#include "io.h"
#include "molecule.h"
#include "geometry.h"
#include "rotate.h"
#include "bonds.h"
#include <fstream>
#include "spdlog/spdlog.h"
#include "aux/eigensupport.h"

namespace Faunus {

// global instances available throughout Faunus.
std::vector<MoleculeData> molecules;
std::vector<ReactionData> reactions;

int &MoleculeData::id() { return _id; }

const int &MoleculeData::id() const { return _id; }

ParticleVector MoleculeData::getRandomConformation(Geometry::GeometryBase &geo, ParticleVector otherparticles) {
    assert(inserterFunctor != nullptr);
    return inserterFunctor(geo, otherparticles, *this);
}

void MoleculeData::loadConformation(const std::string &file, bool keepcharges) {
    ParticleVector v;
    if (loadStructure(file, v, false, keepcharges)) {
        if (keeppos == false)
            Geometry::cm2origo(v.begin(), v.end()); // move to origo
        conformations.push_back(v);
        for (auto &p : v) // add atoms to atomlist
            atoms.push_back(p.id);
    }
    if (v.empty())
        throw std::runtime_error("Structure " + structure + " not loaded. Filetype must be .aam/.pqr/.xyz");
}

void MoleculeData::setInserter(const MoleculeData::TinserterFunc &ifunc) { inserterFunctor = ifunc; }

MoleculeData::MoleculeData() { setInserter(RandomInserter()); }

void to_json(json &j, const MoleculeData &a) {
    j[a.name] = {{"activity", a.activity / 1.0_molar},
                 {"atomic", a.atomic},
                 {"id", a.id()},
                 {"insdir", a.insdir},
                 {"insoffset", a.insoffset},
                 {"keeppos", a.keeppos},
                 {"keepcharges", a.keepcharges},
                 {"bondlist", a.bonds},
                 {"rigid", a.rigid},
                 {"rotate", a.rotate}};
    if (not a.structure.empty())
        j[a.name]["structure"] = a.structure;

    j[a.name]["atoms"] = json::array();
    for (auto id : a.atoms)
        j[a.name]["atoms"].push_back(atoms.at(id).name);
}

void MoleculeData::createMolecularConformations(SingleUseJSON &val) {
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

/**
 * @todo make more readable be splitting into lambdas (since c++ function in function is impossible)
 */
void from_json(const json &j, MoleculeData &a) {
    try {
        if (j.is_object() == false || j.size() != 1)
            throw std::runtime_error("invalid json");
        for (auto it : j.items()) {
            a.name = it.key();
            SingleUseJSON val = it.value(); // keys are deleted after access
            a.insoffset = val.value("insoffset", a.insoffset);
            a.activity = val.value("activity", a.activity) * 1.0_molar;
            a.keeppos = val.value("keeppos", a.keeppos);
            a.keepcharges = val.value("keepcharges", a.keepcharges);
            a.atomic = val.value("atomic", a.atomic);
            a.insdir = val.value("insdir", a.insdir);
            a.bonds = val.value("bondlist", a.bonds);
            a.rigid = val.value("rigid", a.rigid);
            a.rotate = val.value("rotate", true);
            a.id() = val.value("id", a.id());

            if (a.atomic) {
                // read `atoms` list of atom names and convert to atom id's
                for (auto &i : val.at("atoms").get<std::vector<std::string>>()) {
                    auto it = findName(atoms, i);
                    if (it == atoms.end())
                        throw std::runtime_error("unknown atoms in 'atoms'\n");
                    a.atoms.push_back(it->id());
                }
                assert(!a.atoms.empty());
                assert(a.bonds.empty() && "bonds undefined for atomic groups");

                // generate config
                ParticleVector v;
                v.reserve(a.atoms.size());
                for (auto id : a.atoms) {
                    Particle _p;
                    _p = atoms.at(id);
                    v.push_back(_p);
                }
                if (!v.empty())
                    a.conformations.push_back(v);
            }      // done handling atomic groups
            else { // molecular groups
                if (val.count("structure") > 0) {
                    json _struct = val["structure"s];

                    // `structure` is a file name
                    if (_struct.is_string()) // structure from file
                        a.loadConformation(_struct.get<std::string>(), a.keepcharges);

                    else if (_struct.is_object()) {
                        // `structure` is a fasta sequence
                        if (_struct.count("fasta")) {
                            Potential::HarmonicBond bond; // harmonic bond
                            bond.from_json(_struct);      // read 'k' and 'req' from json
                            std::string fasta = _struct.at("fasta").get<std::string>();
                            auto v = Faunus::fastaToParticles(fasta, bond.req);
                            if (not v.empty()) {
                                a.conformations.push_back(v);
                                for (auto &p : v)
                                    a.atoms.push_back(p.id);
                                // connect all atoms with harmonic bonds
                                for (int i = 0; i < (int)v.size() - 1; i++) {
                                    bond.index = {i, i + 1};
                                    a.bonds.push_back(bond.clone());
                                }
                            }
                        } // end of fasta handling
                    }

                    // `structure` is a list of atom positions
                    else if (_struct.is_array()) { // structure is defined inside json
                        std::vector<Particle> v;
                        a.atoms.clear();
                        v.reserve(_struct.size());
                        for (auto &m : _struct) {
                            if (m.is_object()) {
                                if (m.size() == 1) {
                                    for (auto &i : m.items()) {
                                        auto it = findName(atoms, i.key());
                                        if (it == atoms.end())
                                            throw std::runtime_error("unknown atoms in 'structure'");
                                        v.push_back(*it);         // set properties from atomlist
                                        v.back().pos = i.value(); // set position
                                        a.atoms.push_back(it->id());
                                    }
                                }
                            }
                        }
                        if (v.empty())
                            throw std::runtime_error("invalid 'structure' format");
                        a.conformations.push_back(v);
                    } // end of position parser
                }     // end of `structure`

            } // done handling molecular groups

            a.createMolecularConformations(val);

            // pass information to inserter
            auto ins = RandomInserter();
            ins.dir = a.insdir;
            ins.rotate = a.rotate;
            ins.offset = a.insoffset;
            ins.keeppos = a.keeppos;
            a.setInserter(ins);

            // assert that all bonds are *internal*
            for (auto &bond : a.bonds) {
                for (int i : bond->index) {
                    if (i >= a.atoms.size() || i < 0)
                        throw std::runtime_error("bonded atom index " + std::to_string(i) + " out of range");
                }
            }
            // at this stage all given keys should have been accessed or "spend". If any are
            // left, an exception will be thrown.
            if (not val.empty())
                throw std::runtime_error("unused key(s):\n"s + val.dump() + usageTip["moleculelist"]);
        }
    } catch (std::exception &e) {
        throw std::runtime_error("JSON->molecule: " + a.name + ": " + e.what());
    }
}

void from_json(const json &j, std::vector<MoleculeData> &v) {
    v.reserve(v.size() + j.size());
    for (auto &i : j) {
        v.push_back(i);
        v.back().id() = v.size() - 1; // id always match vector index
    }
}

ParticleVector RandomInserter::operator()(Geometry::GeometryBase &geo, const ParticleVector &, MoleculeData &mol) {
    int cnt = 0;
    QuaternionRotate rot;
    bool containerOverlap; // true if container overlap detected

    if (std::fabs(geo.getVolume()) < 1e-20)
        throw std::runtime_error("geometry has zero volume");

    ParticleVector v = mol.conformations.get(); // get random, weighted conformation
    confindex = mol.conformations.index; // lastest index

    do {
        if (cnt++ > maxtrials)
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
            if (keeppos) {        // keep original positions (no rotation/trans)
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
        if (allowoverlap == false) {
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
bool Conformation::empty() const {
    if (positions.empty())
        if (charges.empty())
            return true;
    return false;
}

ParticleVector &Conformation::toParticleVector(ParticleVector &p) const {
    assert(not p.empty() and not empty());
    // copy positions
    if (positions.size() == p.size())
        for (size_t i = 0; i < p.size(); i++)
            p[i].pos = positions[i];
    // copy charges
    if (charges.size() == p.size())
        for (size_t i = 0; i < p.size(); i++)
            p[i].charge = charges[i];
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
    if (_reagid_m.count(molid) == 0)
        if (_prodid_m.count(molid) == 0)
            return false;
    return true;
}
const ReactionData::Tmap &ReactionData::Molecules2Add(bool forward) const { return (forward) ? _prodid_m : _reagid_m; }
const ReactionData::Tmap &ReactionData::Atoms2Add(bool forward) const { return (forward) ? _prodid_a : _reagid_a; }

void from_json(const json &j, ReactionData &a) {

    if (j.is_object() == false || j.size() != 1)
        throw std::runtime_error("Invalid JSON data for ReactionData");

    for (auto &m : Faunus::molecules)
        for (auto &a : Faunus::atoms)
            if (m.name == a.name)
                throw std::runtime_error("Molecules and atoms nust have different names");

    for (auto it = j.begin(); it != j.end(); ++it) {
        a.name = it.key();
        auto &val = it.value();
        a.canonic = val.value("canonic", false);
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
                 {"products", a._prod},
                 {"reactants", a._reag}};
} //!< Serialize to JSON object

} // namespace Faunus
