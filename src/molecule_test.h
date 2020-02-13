#include "molecule.h"

namespace Faunus {

using doctest::Approx;

TEST_SUITE_BEGIN("Molecule");

TEST_CASE("[Faunus] Structure") {
    MoleculeStructureReader sf;
    ParticleVector particles;

    SUBCASE("[Faunus] Structure JSON") {
        atoms = R"([{"Na": {}}, {"Cl": {}}, {"M": {}}])"_json.get<decltype(atoms)>();
        auto j = R"([ {"Na": [0,0,0]}, {"Cl": [1,0,0]}, {"M": [0,4.2,0]} ])"_json;
        particles.clear();
        sf.readJson(particles, j);
        CHECK_EQ(particles.size(), 3);
        CHECK_EQ(particles.front().id, 0);
        CHECK_EQ(particles.back().pos, Point {0, 4.2, 0});
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

TEST_CASE("NeighboursGenerator") {
    BasePointerVector<Potential::BondData> bonds;
    std::vector<std::pair<int, int>> pairs;

    SUBCASE("Linear Chain") {
        // decamer connected with harmonic bonds
        const auto mer = 10;
        for (auto i = 0; i < mer - 1; ++i) {
            const auto j = i + 1;
            bonds.emplace_back<Potential::HarmonicBond>(0., 0., std::vector<int>{i, j});
        }

        const int distance = 3;
        NeighboursGenerator generator(bonds);
        generator.generatePairs(pairs, distance);
        CHECK_EQ(pairs.size(), (mer - 1) + (mer - 2) + (mer - 3));
        auto pairs_matched = [pairs]() -> int {
            int match_cnt = 0;
            for (auto dist = 1; dist <= distance; ++dist) {
                for (auto n = 0; n < mer - dist; ++n) {
                    if( std::find(pairs.begin(), pairs.end(), std::make_pair(n, n + dist)) !=
                    pairs.end()) {
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
            bonds.emplace_back<Potential::HarmonicBond>(0., 0., std::vector<int>{i, j});
        }

        const auto distance = 3; // up to dihedrals
        NeighboursGenerator generator(bonds);
        generator.generatePairs(pairs, distance);
        CHECK_EQ(pairs.size(), mer + mer + mer / 2); // 1-4 pairs in the cyclic hexamer are double degenerated
        auto pairs_matched = [pairs]() -> int {
          int match_cnt = 0;
            for (auto dist = 1; dist <= distance; ++dist) {
                for (auto n = 0; n < mer; ++n) {
                    auto i = n;
                    auto j = (n + dist) % mer;
                    if (i > j) {
                        std::swap(i, j);
                        if(j - i == 3) {
                            continue; // skip the pair doubles in the cyclic hexamer, e.g., 1-4 and 4-1
                        }
                    }
                    if(std::find(pairs.begin(), pairs.end(), std::make_pair(i, j)) != pairs.end()) {
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
        for (auto bond_ndxs : bonds_ndxs) {
            bonds.emplace_back<Potential::HarmonicBond>(0., 0., bond_ndxs);
        }
        NeighboursGenerator generator(bonds);
        generator.generatePairs(pairs, 2);
        CHECK_EQ(pairs.size(), 4 + 4);
    }

    SUBCASE("Harmonic and FENE") {
        bonds.emplace_back<Potential::HarmonicBond>(0., 0., std::vector<int> {1,2});
        bonds.emplace_back<Potential::FENEBond>(0., 0., std::vector<int> {2,3});
        bonds.emplace_back<Potential::FENEBond>(0., 0., std::vector<int> {4,5});
        bonds.emplace_back<Potential::HarmonicBond>(0., 0., std::vector<int> {4,5});
        NeighboursGenerator generator(bonds);
        generator.generatePairs(pairs, 7);
    }
}

TEST_CASE("[Faunus] ExclusionsVicinity") {
    std::vector<std::pair<int, int>> pairs{{0,1}, {1,2}, {1,3}, {6,7}};
    auto exclusions = ExclusionsVicinity::create(10, pairs);

    CHECK(exclusions.isExcluded(1,3));
    CHECK(exclusions.isExcluded(6,7));
    CHECK_FALSE(exclusions.isExcluded(2,3));
    CHECK_FALSE(exclusions.isExcluded(4,5));
    CHECK_FALSE(exclusions.isExcluded(8,9));
}

TEST_CASE("[Faunus] MoleculeData") {
//    json j = R"(
//            { "moleculelist": [
//                { "B": {"activity":0.2, "atomic":true, "insdir": [0.5,0,0], "insoffset": [-1.1, 0.5, 10], "atoms":["a"] } },
//                { "A": { "atomic":false, "structure": [{"a": [0.1, 0.1, 0.1]}] } }
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
        //CHECK(m.insdir == Point(0.5, 0, 0));
        //CHECK(m.insoffset == Point(-1.1, 0.5, 10));
    }
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

TEST_CASE("[Faunus] Conformation") {
    ParticleVector p(1);
    Conformation c;
    CHECK(c.empty());

    c.positions.push_back({1, 2, 3});
    c.charges.push_back(0.5);
    CHECK(not c.empty());

    c.toParticleVector(p);

    CHECK(p[0].pos == Point(1, 2, 3));
    CHECK(p[0].charge == 0.5);
}

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
    CHECK(r.front().reaction == "A = B");
    CHECK(r.front().lnK == Approx(-10.051 - std::log(0.2)));
}

TEST_SUITE_END();
} // namespace Faunus
