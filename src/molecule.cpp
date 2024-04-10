#include "io.h"
#include "molecule.h"
#include "geometry.h"
#include "rotate.h"
#include "bonds.h"
#include "spdlog/spdlog.h"
#include "aux/eigensupport.h"
#include <functional>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <utility>
#include <doctest/doctest.h>

namespace Faunus {

TEST_SUITE_BEGIN("Molecule");

// global instances available throughout Faunus.
std::vector<MoleculeData> molecules; //!< List of molecule types
std::vector<ReactionData> reactions; //!< List of reactions

MoleculeData::MoleculeData(const std::string& name, const ParticleVector& particles,
                           const BasePointerVector<pairpotential::BondData>& bonds)
    : name(name)
    , bonds(bonds) {
    for (const auto& particle : particles) {
        atoms.push_back(particle.id);
    }
    if (!particles.empty()) {
        conformations.push_back(particles);
    }
}

MoleculeData::index_type& MoleculeData::id() { return _id; }

const MoleculeData::index_type& MoleculeData::id() const { return _id; }

bool MoleculeData::isImplicit() const { return implicit; }

bool MoleculeData::isMolecular() const { return !atomic; }

bool MoleculeData::isAtomic() const { return atomic; }

ParticleVector MoleculeData::getRandomConformation(Geometry::GeometryBase& geo, const ParticleVector& otherparticles) {
    assert(inserter != nullptr);
    return (*inserter)(geo, *this, otherparticles);
}

void MoleculeData::setInserter(std::shared_ptr<MoleculeInserter> ins) { inserter = std::move(ins); }

MoleculeData::MoleculeData() { setInserter(std::make_shared<RandomInserter>()); }

void to_json(json &j, const MoleculeData &a) {
    j[a.name] = {
        {"id", a.id()},
        {"atomic", a.atomic},
        {"rigid", a.rigid},
        {"compressible", a.compressible},
        {"implicit", a.isImplicit()},
        {"activity", a.activity / 1.0_molar},
    };
    j[a.name].update(a.json_cfg);

    if(a.inserter != nullptr) {
        a.inserter->to_json(j[a.name]);
    }
    j[a.name]["atoms"] = json::array();
    for (const auto id : a.atoms) {
        j[a.name]["atoms"].push_back(atoms.at(id).name);
    }
    if (!a.bonds.empty()) {
        j[a.name]["bondlist"] = a.bonds;
    }
    if (!a.exclusions.empty()) {
        j[a.name]["exclusionlist"] = a.exclusions;
    }
}

void MoleculeData::createMolecularConformations(const json &j) {
    assert(j.is_object());
    if (auto trajfile = j.value("traj", ""s); not trajfile.empty()) {
        conformations.clear();                                  // remove all previous conformations
        if (j.contains("structure")) {
            throw ConfigurationError("`structure` and `traj` are mutually exclusive");
        }
        try {
            PQRTrajectoryReader::loadTrajectory(trajfile, conformations.data); // read traj. from disk
            if (bool read_charges = j.value("keepcharges", true); read_charges) {
                faunus_logger->debug("charges in {} preferred over `atomlist` values", trajfile);
            } else {
                faunus_logger->debug("replacing all charges from {} with atomlist values", trajfile);
                auto all_particles = conformations.data | ranges::cpp20::views::join;
                Faunus::applyAtomDataCharges(all_particles.begin(), all_particles.end());
            }
            // make sure all conformations are initially placed in the center of the box which
            // in order to reduce PBC problems when inserting later on
            for (ParticleVector& conformation : conformations.data) {
                Geometry::translateToOrigin(conformation.begin(), conformation.end());
            }
        } catch (std::exception& e) {
            throw ConfigurationError("error loading {}: {}", trajfile, e.what());
        }
        if (conformations.empty()) {
            throw ConfigurationError("{} not loaded or empty", trajfile);
        }
        faunus_logger->debug("{} conformations loaded from {}", conformations.size(), trajfile);

        // create atom list
        atoms.clear();
        atoms.reserve(conformations.data.front().size());
        for (const Particle& particle : conformations.data.front()) { // add atoms to atomlist
            atoms.push_back(particle.id);
        }

        // center mass center for each frame to origo assuming whole molecules
        if (j.value("trajcenter", false)) {
            faunus_logger->info("centering conformations in {}", trajfile);
            for (auto& particles : conformations.data) { // loop over conformations
                Geometry::translateToOrigin(particles.begin(), particles.end());
            }
        }

        setConformationWeights(j);
    }
}
void MoleculeData::setConformationWeights(const json& j) {
    std::vector<float> weights(conformations.size(), 1.0); // default uniform weight

    if (auto filename = j.value("trajweight", ""s); !filename.empty()) {
        std::ifstream stream(filename);
        if (!stream) {
            throw ConfigurationError("{} not found", filename);
        }
        weights.clear();
        weights.reserve(conformations.size());
        float weight = 1.0;
        while (stream >> weight) {
            weights.push_back(weight);
        }
        stream.close();
        if (weights.size() != conformations.size()) {
            throw ConfigurationError("{} conformation weights found while expecting {}", weights.size(),
                                     conformations.size());
        }
        faunus_logger->info("{} weights loaded from {}", weights.size(), filename);
    }
    conformations.setWeight(weights.begin(), weights.end());
}

TEST_CASE("[Faunus] MoleculeData") {
    //    json j = R"(
    //            { "moleculelist": [
    //                { "B": {"activity":0.2, "atomic":true, "insdir": [0.5,0,0], "insoffset": [-1.1, 0.5, 10],
    //                "atoms":["a"] } }, { "A": { "atomic":false, "structure": [{"a": [0.1, 0.1, 0.1]}] } }
    //            ]})"_json;
    json j = R"([
                { "B": {"activity": 0.2, "atomic": true, "atoms": ["a"] } },
                { "A": { "atomic": false, "structure": [{"a": [0.1, 0.1, 0.1]}] } }
            ])"_json;

    SUBCASE("Unknown atom") {
        atoms.clear();
        CHECK_THROWS_AS_MESSAGE(j.get<decltype(molecules)>(), std::runtime_error,
                                "JSON->molecule B: unknown atom in atomic molecule");
    }
    SUBCASE("Construction") {
        using doctest::Approx;
        json ja = R"([{"a": {}}])"_json;
        atoms = ja.get<decltype(atoms)>();
        molecules = j.get<decltype(molecules)>(); // fill global instance

        CHECK_EQ(molecules.size(), 2);
        CHECK_EQ(molecules.front().id(), 0);
        CHECK_EQ(molecules.front().name, "B"); // alphabetic order in std::map
        CHECK_EQ(molecules.front().atomic, true);
        CHECK_EQ(molecules.back().id(), 1);
        CHECK_EQ(molecules.back().name, "A"); // alphabetic order in std::map
        CHECK_EQ(molecules.back().atomic, false);

        MoleculeData m = json(molecules.front()); // moldata --> json --> moldata

        CHECK_EQ(m.name, "B");
        CHECK_EQ(m.id(), 0);
        CHECK_EQ(m.activity, Approx(0.2_molar));
        CHECK_EQ(m.atomic, true);
        // CHECK(m.insdir == Point(0.5, 0, 0));
        // CHECK(m.insoffset == Point(-1.1, 0.5, 10));
    }
}

// ============ NeighboursGenerator ============

NeighboursGenerator::NeighboursGenerator(const BondVector &bonds) {
    createBondMap(bonds);
}

void NeighboursGenerator::createBondMap(const BondVector &bonds) {
    for (auto& bond : bonds.find<pairpotential::StretchData>()) {
        auto add_bond = [this, &bond](auto i, auto j) {
            auto atom_emplaced = bond_map.emplace(bond->indices.at(i), AtomList{}); // find or create empty vector
            atom_emplaced.first->second.push_back(
                bond->indices.at(j)); // atom_emplaced.<iterator>-><vector>.push_back()
        };
        const auto head_ndx = 0;
        const auto tail_ndx = 1;
        add_bond(head_ndx, tail_ndx); // add the bond to the map in both head-tail
        add_bond(tail_ndx, head_ndx); // and tail-head directions
    }
}

void NeighboursGenerator::generatePaths(int bonded_distance) {
    if(paths.empty()) {
        // starting with a path of length 0 bonds which are just separated atoms
        paths.emplace_back();
        for (const auto& bonded_atoms : bond_map) {
            paths.back().emplace(AtomList{bonded_atoms.first});
        }
    }
    // continue wherever the last path generation ended
    for (int distance = (int)paths.size(); distance <= bonded_distance; ++distance) {
        paths.emplace_back();
        const auto base_paths = paths.rbegin()[1]; // last but one: a set of paths of length (distance - 1)
        for (const auto& base_path : base_paths) {
            // try to extend the path with every atom bonded to the tail
            for (const auto& appending : bond_map.at(base_path.back())) {
                // avoid loops: the appending atom shall not be anywhere in the path yet
                if (std::find(base_path.begin(), base_path.end(), appending) == base_path.end()) {
                    auto path(base_path);
                    path.push_back(appending);
                    paths.back().emplace(path);
                }
            }
        }
    }
}

void NeighboursGenerator::generatePairs(AtomPairList& pairs, int bond_distance) {
    using namespace ranges::cpp20::views;
    generatePaths(bond_distance);
    auto is_valid_pair = [](const auto& pair) { return pair.first != pair.second; };
    auto make_ordered_pair = [](const auto& path) {
        return path.front() < path.back() ? AtomPair(path.front(), path.back()) : AtomPair(path.back(), path.front());
    };
    auto excluded_pairs = // temp. conversion to std::set is avoid duplicates
        paths | join | transform(make_ordered_pair) | filter(is_valid_pair) | ranges::to<std::set<AtomPair>>;
    pairs.reserve(pairs.size() + excluded_pairs.size());
    std::copy(excluded_pairs.begin(), excluded_pairs.end(), std::back_inserter(pairs));
}

TEST_CASE("NeighboursGenerator") {
    BasePointerVector<pairpotential::BondData> bonds;
    std::vector<std::pair<int, int>> pairs;

    SUBCASE("Linear Chain") {
        // decamer connected with harmonic bonds
        const auto mer = 10;
        for (auto i = 0; i < mer - 1; ++i) {
            const auto j = i + 1;
            bonds.emplace_back<pairpotential::HarmonicBond>(0.0, 0.0, std::vector<int>{i, j});
        }

        const int distance = 3;
        NeighboursGenerator generator(bonds);
        generator.generatePairs(pairs, distance);
        CHECK_EQ(pairs.size(), (mer - 1) + (mer - 2) + (mer - 3));
        auto pairs_matched = [pairs]() {
            int match_cnt = 0;
            for (auto dist = 1; dist <= distance; ++dist) {
                for (auto n = 0; n < mer - dist; ++n) {
                    if (std::find(pairs.begin(), pairs.end(), std::make_pair(n, n + dist)) != pairs.end()) {
                        ++match_cnt;
                    }
                }
            }
            return match_cnt;
        };
        CHECK_EQ(pairs_matched(), pairs.size());
    }

    SUBCASE("Cycle") {
        // cyclic hexamer connected with harmonic bonds
        const auto mer = 6;
        for (auto i = 0; i < mer; ++i) {
            auto const j = (i + 1) % mer;
            bonds.emplace_back<pairpotential::HarmonicBond>(0., 0., std::vector<int>{i, j});
        }

        const auto distance = 3; // up to dihedrals
        NeighboursGenerator generator(bonds);
        generator.generatePairs(pairs, distance);
        CHECK_EQ(pairs.size(), mer + mer + mer / 2); // 1-4 pairs in the cyclic hexamer are double degenerated
        auto pairs_matched = [pairs]() {
            int match_cnt = 0;
            for (auto dist = 1; dist <= distance; ++dist) {
                for (auto n = 0; n < mer; ++n) {
                    auto i = n;
                    auto j = (n + dist) % mer;
                    if (i > j) {
                        std::swap(i, j);
                        if (j - i == 3) {
                            continue; // skip the pair doubles in the cyclic hexamer, e.g., 1-4 and 4-1
                        }
                    }
                    if (std::find(pairs.begin(), pairs.end(), std::make_pair(i, j)) != pairs.end()) {
                        ++match_cnt;
                    }
                }
            }
            return match_cnt;
        };
        CHECK_EQ(pairs_matched(), pairs.size());
    }

    SUBCASE("Branched") {
        // isopentane like structure
        std::vector<std::vector<int>> bonds_ndxs = {{0, 1}, {1, 2}, {1, 3}, {3, 4}};
        for (const auto& bond_ndxs : bonds_ndxs) {
            bonds.emplace_back<pairpotential::HarmonicBond>(0., 0., bond_ndxs);
        }
        NeighboursGenerator generator(bonds);
        generator.generatePairs(pairs, 2);
        CHECK_EQ(pairs.size(), 4 + 4);
    }

    SUBCASE("Harmonic and FENE") {
        bonds.emplace_back<pairpotential::HarmonicBond>(0., 0., std::vector<int>{1, 2});
        bonds.emplace_back<pairpotential::FENEBond>(0., 0., std::vector<int>{2, 3});
        bonds.emplace_back<pairpotential::FENEBond>(0., 0., std::vector<int>{4, 5});
        bonds.emplace_back<pairpotential::HarmonicBond>(0., 0., std::vector<int>{4, 5});
        NeighboursGenerator generator(bonds);
        generator.generatePairs(pairs, 7);
    }
}

// ============ MoleculeBuilder ============

void MoleculeBuilder::from_json(const json &j, MoleculeData &molecule) {
    if (is_used) {
        throw std::logic_error("MoleculeBuilder cannot be used twice");
    }
    is_used = true;
    try {
        const auto& [key, j_properties] = jsonSingleItem(j);
        molecule.name = molecule_name = key;
        molecule.id() = j_properties.value("id", molecule.id());
        molecule.atomic = j_properties.value("atomic", molecule.atomic);
        molecule.rigid = j_properties.value("rigid", molecule.rigid);
        molecule.compressible = j_properties.value("compressible", molecule.compressible);
        molecule.activity = j_properties.value("activity", molecule.activity / 1.0_molar) * 1.0_molar;
        molecule.implicit = j_properties.value("implicit", false);
        if (molecule.implicit and molecule.atomic) {
            throw std::runtime_error("atomic molecules cannot be implicit");
        }

        readCompoundValues(j_properties);
        for (const auto& particle : particles) {
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
        if (auto it = j_properties.find("structure"); it != j_properties.end()) {
            if (it->is_string()) {
                molecule.json_cfg["structure"] = it->get<std::string>();
            }
        }
        molecule.json_cfg["keepcharges"] = j_properties.value("keepcharges", true);

        // at this stage all given keys should have been accessed or "spend". If any are
        // left, an exception will be thrown.
        //        if (! val.empty()) {
        //            throw std::runtime_error("unused key(s):\n"s + val.dump() + usageTip["moleculelist"]);
        //        }
    } catch (std::exception &e) {
        usageTip.pick("moleculelist");
        throw ConfigurationError("molecule '{}': {}", molecule_name, e.what());
    }
}

std::shared_ptr<MoleculeInserter> MoleculeBuilder::createInserter(const json& j) {
    return std::make_shared<RandomInserter>(j);
}

void MoleculeBuilder::readCompoundValues(const json &j) {
    auto is_atomic = j.value("atomic", false);
    if (is_atomic) {
        readAtomic(j);
    } else {
        readParticles(j);
        readBonds(j);
        if (isFasta(j)) {
            readFastaBonds(j);
        }
        readExclusions(j);
    }
}

void MoleculeBuilder::readAtomic(const json& j) {
    if (const auto atomlist = j.find("atoms"); atomlist != j.end()) {
        if (!atomlist->is_array()) {
            throw ConfigurationError("`atoms` must be an array");
        }
        try {
            particles.reserve(atomlist->size());
            for (const auto& atomname : *atomlist) {
                const auto& atomdata = findAtomByName(atomname.get<std::string>());
                particles.emplace_back(atomdata);
            }
        } catch (std::exception& e) {
            throw ConfigurationError("molecule '{}': {}", molecule_name, e.what());
        }
    }
}

void MoleculeBuilder::readParticles(const json& j) {
    auto structure = j.find("structure");
    if (structure != j.end()) {
        bool read_charges = j.value("keepcharges", true);
        MoleculeStructureReader structure_reader(read_charges);
        structure_reader.readJson(particles, *structure);
        if (j.value("ensphere", false)) {
            particles = Geometry::mapParticlesOnSphere(particles);
        }
        if (j.value("to_disk", false)) {
            PQRWriter().save(molecule_name + "-initial.pqr", particles.begin(), particles.end(), Point(0, 0, 0));
        }
    } else {
        // allow virtual molecules :-/
        // shall we rather try to fallback on readAtomic()?
        // throw ConfigurationError("structure of the molecule not given");
    }
}

void MoleculeBuilder::readBonds(const json& j) {
    using namespace ranges::cpp20;
    bonds = j.value("bondlist", bonds);
    auto is_invalid_index = [size = particles.size()](auto& i) { return i >= size || i < 0; };
    auto indices = bonds | views::transform(&pairpotential::BondData::indices) | views::join;
    if (any_of(indices, is_invalid_index)) {
        throw ConfigurationError("bonded index out of range");
    }
}

void MoleculeBuilder::readFastaBonds(const json& j) {
    if (particles.size() < 2) {
        return;
    }
    bonds.vec.reserve(bonds.size() + particles.size() - 1);
    pairpotential::HarmonicBond bond;
    bond.from_json(j.at("structure")); // read 'k' and 'req' from json
    for (int i = 1; i < (int)particles.size(); ++i) {
        bond.indices = {i - 1, i};
        bonds.push_back(bond.clone());
    }
}

void MoleculeBuilder::readExclusions(const json& j) {
    auto excluded_neighbours_distance = j.value("excluded_neighbours", 0);
    if (excluded_neighbours_distance > 0) {
        NeighboursGenerator generator(bonds);
        generator.generatePairs(exclusion_pairs, excluded_neighbours_distance);
    }
    for (auto pair : j.value("exclusionlist", json::array())) {
        if (!pair.is_array() || pair.size() != 2) {
            throw ConfigurationError("unrecognized molecule's exclusion format");
        }
        exclusion_pairs.emplace_back(pair[0].get<int>(), pair[1].get<int>());
    }
}

bool MoleculeBuilder::isFasta(const json& j) {
    auto structure = j.find("structure");
    return structure != j.end() && structure->find("fasta") != structure->end();
}

TEST_CASE("[Faunus] MoleculeBuilder") {
    atoms = R"([{"ALA": {}}, {"GLY": {}}, {"NTR": {}}, {"CTR": {}}])"_json.get<decltype(atoms)>();
    auto j = R"(
        {"peptide": { "excluded_neighbours": 2,
        "structure": {"fasta": "nAAAAGGc", "k": 3, "req": 7},
        "exclusionlist": [[3,6]],
        "bondlist": [{"harmonic": {"index": [0, 7], "k": 3, "req": 7}}] }}
    )"_json;

    auto molecule = j.get<MoleculeData>();

    REQUIRE_EQ(molecule.atoms.size(), 8);
    REQUIRE_EQ(molecule.bonds.size(), 7 + 1);
    CHECK(molecule.isPairExcluded(0, 1));
    CHECK(molecule.isPairExcluded(0, 2));
    CHECK_FALSE(molecule.isPairExcluded(0, 3));
    CHECK(molecule.isPairExcluded(7, 5));
    CHECK_FALSE(molecule.isPairExcluded(7, 4));
    CHECK(molecule.isPairExcluded(3, 6));
    CHECK(molecule.isPairExcluded(6, 3));
    CHECK_FALSE(molecule.isPairExcluded(3, 7));
    CHECK(molecule.isPairExcluded(7, 1));
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

void MoleculeStructureReader::readFile(ParticleVector& particles, const std::string& filename) const {
    faunus_logger->info("Reading molecule configuration from file: {}", filename);
    particles = Faunus::loadStructure(filename, read_charges); // throws if nothing is loaded!
}

void MoleculeStructureReader::readArray(ParticleVector& particles, const json& j) {
    particles.reserve(j.size());
    for (const auto& particle : j) {
        const auto& [atomname, position] = jsonSingleItem(particle);
        particles.emplace_back(findAtomByName(atomname), position.get<Point>());
    }
}

void MoleculeStructureReader::readFasta(ParticleVector& particles, const json& j) {
    if (auto it = j.find("fasta"); it != j.end()) {
        pairpotential::HarmonicBond bond; // harmonic bond
        bond.from_json(j);            // read 'k' and 'req' from json

        auto fasta = it->get<std::string>(); // fasta sequence or filename

        // is `fasta` a valid filename?
        if ("fasta" == fasta.substr(fasta.find_last_of('.') + 1)) {
            if (auto stream = std::ifstream(fasta)) {
                fasta = CoarseGrainedFastaFileReader(0.0).loadSequence(stream);
                stream.close();
            } else {
                throw ConfigurationError("could not open fasta file: {}", fasta);
            }
        }
        particles = Faunus::fastaToParticles(fasta, bond.equilibrium_distance);
        faunus_logger->debug("fasta sequence parsed with {} letters", particles.size());
    } else {
        throw ConfigurationError("invalid FASTA format");
    }
}
MoleculeStructureReader::MoleculeStructureReader(bool read_charges) : read_charges(read_charges) {}

void from_json(const json &j, MoleculeData &a) {
    MoleculeBuilder builder;
    builder.from_json(j, a);
}

size_t MoleculeData::numConformations() const { return conformations.data.size(); }

void from_json(const json &j, std::vector<MoleculeData> &v) {
    v.reserve(v.size() + j.size());
    for (const auto& i : j) {
        v.push_back(i);
        v.back().id() = static_cast<MoleculeData::index_type>(v.size() - 1); // id always match vector index
    }
}

TEST_CASE("[Faunus] MoleculeStructureReader") {
    using doctest::Approx;
    MoleculeStructureReader sf;
    ParticleVector particles;

    SUBCASE("[Faunus] Structure JSON") {
        atoms = R"([{"Na": {}}, {"Cl": {}}, {"M": {}}])"_json.get<decltype(atoms)>();
        auto j = R"([ {"Na": [0,0,0]}, {"Cl": [1,0,0]}, {"M": [0,4.2,0]} ])"_json;
        particles.clear();
        sf.readJson(particles, j);
        CHECK_EQ(particles.size(), 3);
        CHECK_EQ(particles.front().id, 0);
        CHECK_EQ(particles.back().pos, Point{0, 4.2, 0});
    }

    SUBCASE("[Faunus] Structure FASTA") {
        atoms = R"([{"ALA": {}}, {"GLY": {}}, {"NTR": {}}, {"CTR": {}}])"_json.get<decltype(atoms)>();
        auto j = R"({"fasta": "nAAAAGGc", "k": 3, "req": 7})"_json;
        particles.clear();
        sf.readJson(particles, j);
        CHECK_EQ(particles.size(), 8);
        CHECK_EQ(particles[4].id, 0);
        CHECK_EQ(particles[5].id, 1);
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

void ExclusionsSimple::add(const std::vector<AtomPair>& exclusions) {
    for (const auto &pair : exclusions) {
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
    excluded_pairs->at(i * size + j) = excluded_pairs->at(j * size + i) = static_cast<char_bool>(true);
}

void from_json(const json &j, ExclusionsSimple &exclusions) {
    auto is_invalid_pair = [](auto& pair) { return !pair.is_array() || pair.size() != 2; };
    if (!j.is_array() || std::any_of(j.begin(), j.end(), is_invalid_pair)) {
        throw ConfigurationError("invalid exclusion: {}", j.dump());
    }
    std::for_each(j.begin(), j.end(), [&](auto& pair) { exclusions.add(pair[0], pair[1]); });
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

ExclusionsVicinity ExclusionsVicinity::create(int atoms_cnt, const std::vector<AtomPair>& pairs) {
    auto max_neighbours_distance = 0;
    auto compare = [](auto &a, auto &b) { return std::abs(a.first - a.second) < std::abs(b.first - b.second); };
    if (auto it = std::max_element(pairs.begin(), pairs.end(), compare); it != pairs.end()) {
        max_neighbours_distance = std::abs(it->first - it->second);
    }
    ExclusionsVicinity exclusions(atoms_cnt, max_neighbours_distance);
    exclusions.add(pairs);
    return exclusions;
}

ExclusionsVicinity::ExclusionsVicinity(int atoms_cnt, int max_difference)
    : atoms_cnt(atoms_cnt), max_bond_distance(max_difference),
      excluded_pairs(std::make_shared<std::vector<char_bool>>()) {
    faunus_logger->log(atoms_cnt * max_difference < 1'000'000 ? spdlog::level::trace : spdlog::level::warn,
                       "creating exclusion matrix {}×{} for {} atoms within distance {}", atoms_cnt, max_difference,
                       atoms_cnt, max_difference);

    excluded_pairs->resize(atoms_cnt * max_difference, static_cast<char_bool>(false));
}

void ExclusionsVicinity::add(int i, int j) {
    if (i > j) {
        std::swap(i, j);
    }
    if (j >= atoms_cnt || j - i > max_bond_distance) {
        throw std::range_error("exclusion indexes out of range");
    }
    excluded_pairs->at(toIndex(i, j)) = static_cast<unsigned char>(true);
}

void ExclusionsVicinity::add(const std::vector<AtomPair>& pairs) {
    for (const auto& pair : pairs) {
        add(pair.first, pair.second);
    }
}

void to_json(json &j, const ExclusionsVicinity &exclusions) {
    j = json::array();
    for (auto it = exclusions.excluded_pairs->begin(); it != exclusions.excluded_pairs->end(); ++it) {
        if (static_cast<bool>(*it)) {
            const auto n = std::distance(exclusions.excluded_pairs->begin(), it);
            j.push_back(exclusions.fromIndex(n));
        }
    }
}

TEST_CASE("[Faunus] ExclusionsVicinity") {
    std::vector<ExclusionsVicinity::AtomPair> pairs{{0, 1}, {1, 2}, {1, 3}, {6, 7}};
    auto exclusions = ExclusionsVicinity::create(10, pairs);

    CHECK(exclusions.isExcluded(1, 3));
    CHECK(exclusions.isExcluded(6, 7));
    CHECK_FALSE(exclusions.isExcluded(2, 3));
    CHECK_FALSE(exclusions.isExcluded(4, 5));
    CHECK_FALSE(exclusions.isExcluded(8, 9));
}

void from_json(const json &j, MoleculeInserter &inserter) { inserter.from_json(j); }
void to_json(json &j, const MoleculeInserter &inserter) { inserter.to_json(j); }

/**
 * @param geo Geometry to use for PBC and container overlap check
 * @param molecule Molecular type to insert
 * @param ignored_other_particles Other particles in the system (ignored for this inserter!)
 * @return Inserted particle vector
 */
ParticleVector RandomInserter::operator()(const Geometry::GeometryBase &geo, MoleculeData &molecule,
                                          [[maybe_unused]] const ParticleVector &ignored_other_particles) {
    auto particles = molecule.conformations.sample(random.engine); // random, weighted conformation
    if (particles.empty()) {
        throw std::runtime_error("nothing to insert for molecule '"s + molecule.name + "'");
    }

    QuaternionRotate rotator;
    auto container_overlap = [&geo](const auto& particle) { return geo.collision(particle.pos); };

    if (keep_positions) { // keep given positions
        if (ranges::cpp20::none_of(particles, container_overlap)) {
            return particles;
        }
        throw std::runtime_error("inserted molecule does not fit in container");
    }

    for (auto attempts = 0; attempts < max_trials; attempts++) {
        if (molecule.atomic) {
            translateRotateAtomicGroup(geo, rotator, particles);
        } else {
            translateRotateMolecularGroup(geo, rotator, particles);
        }
        if (allow_overlap || ranges::cpp20::none_of(particles, container_overlap)) {
            return particles;
        }
    };
    throw std::runtime_error("Max. # of overlap checks reached upon insertion.");
}

void RandomInserter::translateRotateMolecularGroup(const Geometry::GeometryBase& geo, QuaternionRotate& rotator,
                                                   ParticleVector& particles) const {
    Geometry::translateToOrigin(particles.begin(), particles.end()); // translate to origin
    if (rotate) {
        rotator.set(2.0 * pc::pi * random(), randomUnitVector(random)); // random rot around random vector
        Geometry::rotate(particles.begin(), particles.end(), rotator.getQuaternion());
        assert(Geometry::massCenter(particles.begin(), particles.end()).norm() < 1e-6); // cm shouldn't move
    }
    Point new_mass_center;
    geo.randompos(new_mass_center, random);                       // random point in container
    new_mass_center = new_mass_center.cwiseProduct(dir) + offset; // add defined dirs (default: 1,1,1)
    Geometry::translate(particles.begin(), particles.end(), new_mass_center, geo.getBoundaryFunc());
}

void RandomInserter::translateRotateAtomicGroup(const Geometry::GeometryBase& geo, QuaternionRotate& rotator,
                                                ParticleVector& particles) const {
    for (auto& particle : particles) { // for each atom type id
        if (rotate) {                  // internal rotation of atomic particles
            rotator.set(2.0 * pc::pi * random(), randomUnitVector(random));
            particle.rotate(rotator.getQuaternion(), rotator.getRotationMatrix());
        }
        geo.randompos(particle.pos, random);
        particle.pos = particle.pos.cwiseProduct(dir) + offset;
        geo.boundary(particle.pos);
    }
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
    j["allow overlap"] = allow_overlap;
}

bool Conformation::empty() const {
    return positions.empty() && charges.empty();
}

/**
 * @param particles Destination particle vector
 * @throw If there's a size mismatch with destination
 *
 * Overwrites positions and charges; remaining properties are untouched
 */
void Conformation::copyTo(ParticleVector &particles) const {
    if (positions.size() != particles.size() or charges.size() != particles.size()) {
        throw std::runtime_error("conformation size mismatch for positions and/or charges");
    }
    /*
     * Note how `std::ranges::transform` accepts data member pointers like `&Particle::pos`
     * which are internally accessed using `std::invoke`.
     */
    auto particle_positions = particles | ranges::cpp20::views::transform(&Particle::pos);
    std::copy(positions.begin(), positions.end(), particle_positions.begin());

    auto particle_charges = particles | ranges::cpp20::views::transform(&Particle::charge);
    std::copy(charges.begin(), charges.end(), particle_charges.begin());
}

TEST_CASE("[Faunus] Conformation") {
    Conformation conformation;
    CHECK(conformation.empty());

    conformation.positions.push_back({1, 2, 3});
    conformation.charges.push_back(0.5);
    CHECK(not conformation.empty());

    ParticleVector particles;
    CHECK_THROWS(conformation.copyTo(particles));
    particles.resize(1);
    conformation.copyTo(particles);
    CHECK(particles[0].pos == Point(1, 2, 3));
    CHECK(particles[0].charge == 0.5);
}

ReactionData::Direction ReactionData::getDirection() const { return direction; }

void ReactionData::setDirection(ReactionData::Direction new_direction) {
    if (new_direction != direction) {
        direction = new_direction;
    }
}

void ReactionData::setRandomDirection(Random& random) {
    auto direction = static_cast<ReactionData::Direction>((char)random.range(0, 1)); // random direction
    setDirection(direction);
}

ReactionData::AtomicAndMolecularPair ReactionData::getProducts() const {
    if (direction == Direction::RIGHT) {
        return {right_atoms, right_molecules};
    }
    return {left_atoms, left_molecules};
}

ReactionData::AtomicAndMolecularPair ReactionData::getReactants() const {
    if (direction == Direction::RIGHT) {
        return {left_atoms, left_molecules};
    }
    return {right_atoms, right_molecules};
}

/**
 * @return Pair of sets with id's for atoms and molecules participating in the reaction (i.e. reactants & products)
 */
std::pair<std::set<int>, std::set<int>> ReactionData::participatingAtomsAndMolecules() const {
    auto extract_keys = [](const auto &map, std::set<int> &keys) { // map keys --> set
        std::transform(map.begin(), map.end(), std::inserter(keys, keys.end()), [](auto &pair) { return pair.first; });
    };
    std::set<int> atomic;
    std::set<int> molecular;
    auto reactants = getReactants();
    auto products = getProducts();
    extract_keys(reactants.first, atomic);     // copy atomic reactants to `atomic`.
    extract_keys(reactants.second, molecular); // copy molecular reactants to `molecular`
    extract_keys(products.first, atomic);      // same for products, merge into same set
    extract_keys(products.second, molecular);  // same for products, merge into same set
    return {atomic, molecular};
}

void ReactionData::reverseDirection() {
    if (direction == Direction::RIGHT) {
        setDirection(Direction::LEFT);
    } else {
        setDirection(Direction::RIGHT);
    }
}

void from_json(const json &j, ReactionData &a) {
    if (!j.is_object() || j.size() != 1) {
        throw ConfigurationError("Invalid JSON data for ReactionData");
    }

    a.direction = ReactionData::Direction::RIGHT;

    for (const auto& molecule : Faunus::molecules) {
        for (const auto& atom : Faunus::atoms) {
            if (molecule.name == atom.name) {
                throw ConfigurationError("Molecules and atoms must have different names");
            }
        }
    }

    for (const auto& [key, val] : j.items()) {
        a.reaction_str = key; // reaction string, e.g. "A + B = C"
        a.only_neutral_molecules = val.value("neutral", false);
        if (val.contains("lnK")) {
            a.lnK_unmodified = val.at("lnK").get<double>();
        } else if (val.contains("pK")) {
            a.lnK_unmodified = -std::numbers::ln10 * val.at("pK").get<double>();
        } else {
            a.lnK_unmodified = 0.0;
        }
        a.lnK = a.lnK_unmodified;

        // helper function used to parse and register atom and molecule names; updates lnK
        auto registerNames = [&](auto& names, auto&& atom_map, auto& mol_map, double sign) {
            for (const auto& atom_or_molecule_name : names) { // loop over species on reactant side (left)
                auto [atom, molecule] = a.findAtomOrMolecule(atom_or_molecule_name);
                if (atom != Faunus::atoms.end()) { // atomic reactants
                    if (atom->implicit) {          // if atom is implicit, multiply K by its activity
                        if (atom->activity > 0.0) {
                            a.lnK += sign * std::log(atom->activity / 1.0_molar);
                        }
                    } else {
                        assert(std::fabs(atom->activity) <= pc::epsilon_dbl);
                        atom_map[atom->id()]++; // increment stoichiometric number
                    }
                } else if (molecule != Faunus::molecules.end()) { // molecular reactants (incl. "atomic" groups)
                    mol_map[molecule->id()]++;                    // increment stoichiometric number
                    if (molecule->activity > 0) { // assume activity not part of K -> divide by activity
                        a.lnK -= sign * std::log(molecule->activity / 1.0_molar);
                    }
                } else {
                    assert(false); // we should never reach here
                }
            }
        }; // end of lambda function

        std::tie(a.left_names, a.right_names) = parseReactionString(a.reaction_str); // lists of species
        registerNames(a.left_names, a.left_atoms, a.left_molecules, 1.0);            // reactants
        registerNames(a.right_names, a.right_atoms, a.right_molecules, -1.0);        // products
    }
}

void to_json(json &j, const ReactionData &reaction) {
    ReactionData a = reaction;
    // we want lnK to show for LEFT-->RIGHT direction
    a.setDirection(ReactionData::Direction::RIGHT);
    j[a.getReactionString()] = {{"lnK", a.lnK_unmodified},
                                {"pK", -a.lnK_unmodified / std::numbers::ln10},
                                {"swap_move", a.containsAtomicSwap()},
                                {"neutral", a.only_neutral_molecules},
                                {"pK'", a.freeEnergy() / std::numbers::ln10}};
}

std::pair<decltype(Faunus::atoms)::const_iterator, decltype(Faunus::molecules)::const_iterator>
ReactionData::findAtomOrMolecule(const std::string& atom_or_molecule_name) {
    auto atom_iter = findName(Faunus::atoms, atom_or_molecule_name);
    auto molecule_iter = findName(Faunus::molecules, atom_or_molecule_name);
    if (molecule_iter == Faunus::molecules.end() and atom_iter == Faunus::atoms.end()) {
        throw std::runtime_error("unknown species '" + atom_or_molecule_name + "'");
    }
    return {atom_iter, molecule_iter};
}

/**
 * @return Reaction free energy in the current direction, i.e. products minus reactants (kT)
 */
double ReactionData::freeEnergy() const { return (direction == Direction::RIGHT) ? -lnK : lnK; }

const std::string& ReactionData::getReactionString() const { return reaction_str; }

bool ReactionData::containsAtomicSwap() const {
    // If exactly one atomic reactant and one atomic products, it's a swap move!
    return (left_atoms.size() == 1 && right_atoms.size() == 1);
}

const ReactionData::MapFilter ReactionData::is_implicit_group = [](const auto& pair) {
    return Faunus::molecules.at(pair.first).isImplicit();
};

const ReactionData::MapFilter ReactionData::not_implicit_group = [](const auto& pair) {
    return !Faunus::molecules.at(pair.first).isImplicit();
};

const ReactionData::MapFilter ReactionData::not_implicit_atom = [](const auto& pair) {
    return !Faunus::atoms.at(pair.first).implicit;
};

const ReactionData::MapFilter ReactionData::is_atomic_group = [](const auto& pair) {
    return Faunus::molecules.at(pair.first).isAtomic();
};

const ReactionData::MapFilter ReactionData::is_molecular_group = [](const auto& pair) {
    return Faunus::molecules.at(pair.first).isMolecular();
};

TEST_CASE("[Faunus] ReactionData") {
    using doctest::Approx;
    json j = R"(
            {
                "atomlist" :
                    [ {"a": { "r":1.1 } } ],
                "moleculelist": [
                    { "A": { "atomic":false, "activity":0.2 } },
                    { "B": { "atomic":true, "atoms":["a"] } }
                ],
                "reactionlist": [
                    {"A = B": {"lnK":-10.051, "canonic":true, "N":100 } }
                ]
            } )"_json;

    Faunus::atoms = j["atomlist"].get<decltype(atoms)>();
    molecules = j["moleculelist"].get<decltype(molecules)>(); // fill global instance

    auto &r = reactions; // reference to global reaction list
    r = j["reactionlist"].get<decltype(reactions)>();

    CHECK(r.size() == 1);
    CHECK(r.front().getReactionString() == "A = B");
    CHECK(r.front().freeEnergy() == Approx(10.051 + std::log(0.2)));
}

void MoleculeInserter::from_json(const json &) {}
void MoleculeInserter::to_json(json &) const {}

TEST_SUITE_END();

UnknownMoleculeError::UnknownMoleculeError(std::string_view molecule_name)
    : GenericError("unknown molecule: '{}'", molecule_name) {}

/**
 * @throw if molecule not found
 */
MoleculeData& findMoleculeByName(std::string_view name) {
    const auto result = findName(Faunus::molecules, name);
    if (result == Faunus::molecules.end()) {
        throw UnknownMoleculeError(name);
    }
    return *result;
}

/**
 * @param process_string Reaction such as "A + B = C"
 * @return Pair with reactants (first) and products (second) as atomic or molecule names
 * @throws std::runtime_error if unknown atom or molecule, or syntax error
 */
std::pair<std::vector<std::string>, std::vector<std::string>> parseReactionString(const std::string& process_string) {
    using Tvec = std::vector<std::string>;
    Tvec names; // vector of atom/molecule names
    std::string atom_or_molecule_name;
    std::istringstream iss(process_string);
    while (iss >> atom_or_molecule_name) { // stream all words into vector
        names.push_back(atom_or_molecule_name);
    }

    names.erase(std::remove(names.begin(), names.end(), "+"), names.end());

    auto it = std::find(names.begin(), names.end(), "=");
    if (it == names.end()) {
        throw std::runtime_error("products and reactants must be separated by ' = '");
    }
    return {Tvec(names.begin(), it), Tvec(it + 1, names.end())};
}

/**
 * @param j Json object with a `molecules` list of strings
 * @return vector of corresponding molecule ids (sorted)
 * @throw if unknown molecule name or not a list of strings
 */
std::vector<MoleculeData::index_type> parseMolecules(const json& j) {
    if (!j.is_array()) {
        throw ConfigurationError("array of molecules expected");
    }
    auto molids = Faunus::names2ids(Faunus::molecules, j.get<std::vector<std::string>>());
    std::sort(molids.begin(), molids.end());
    return molids;
}

} // namespace Faunus
