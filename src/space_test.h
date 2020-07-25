#include "space.h"

namespace Faunus {

using doctest::Approx;
TEST_SUITE_BEGIN("Space");

TEST_CASE("[Faunus] Change") {
    Change change;
    CHECK(not change);
    change.dV = true;
    CHECK(not change.empty());
    CHECK(change);
    change.clear();
    CHECK(change.empty());
}

TEST_CASE("[Faunus] Space") {
    Space spc1;
    spc1.geo = R"( {"type": "sphere", "radius": 1e9} )"_json;

    // check molecule insertion
    Faunus::molecules.at(0).atomic = false;
    REQUIRE_EQ(Faunus::molecules.at(0).atomic, false);
    Faunus::molecules.at(0).atoms.resize(2);
    CHECK(Faunus::molecules.at(0).atoms.size() == 2);

    Faunus::atoms.resize(2);
    CHECK(Faunus::atoms.at(0).mw == 1);

    Particle a;
    a.id = 0;
    a.pos.setZero();
    ParticleVector p(2, a);
    CHECK(p[0].traits().mw == 1);
    p[0].pos.x() = 2;
    p[1].pos.x() = 3;
    spc1.push_back(0, p); // insert molecular group
    CHECK(spc1.p.size() == 2);
    CHECK(spc1.groups.size() == 1);
    CHECK(spc1.groups.back().isMolecular());
    CHECK(spc1.groups.front().id == 0);
    CHECK(spc1.groups.front().cm.x() == doctest::Approx(2.5));

    // check `positions()`
    CHECK(&spc1.positions()[0] == &spc1.p[0].pos);

    SUBCASE("getActiveParticles") {
        // add three groups to space
        Space spc;
        spc.geo = R"( {"type": "sphere", "radius": 1e9} )"_json;
        Particle a;
        a.pos.setZero();
        a.id = 0;
        ParticleVector pvec({a, a, a});

        Faunus::molecules.at(0).atoms.resize(3);

        spc.push_back(0, pvec);
        spc.push_back(0, pvec);
        spc.push_back(0, pvec);

        spc.groups.at(0).id = 0;
        spc.groups.at(1).id = 1;
        spc.groups.at(2).id = 0;

        for (size_t i = 0; i < spc.p.size(); i++)
            spc.p[i].charge = double(i);

        CHECK(spc.p.size() == 9);
        CHECK(spc.groups.size() == 3);

        CHECK(spc.numMolecules<Space::Tgroup::ANY>(0) == 2);
        CHECK(spc.numMolecules<Space::Tgroup::ANY>(1) == 1);

        spc.groups[0].deactivate(spc.p.begin(), spc.p.begin() + 1);
        spc.groups[1].deactivate(spc.groups[1].begin(), spc.groups[1].end());

        CHECK(spc.groups[0].size() == 2);
        CHECK(spc.groups[1].size() == 0);
        CHECK(spc.groups[2].size() == 3);

        CHECK(spc.numMolecules<Space::Tgroup::ACTIVE>(0) == 2);
        CHECK(spc.numMolecules<Space::Tgroup::ACTIVE | Space::Tgroup::NEUTRAL>(0) == 0);
        CHECK(spc.numMolecules<Space::Tgroup::ACTIVE>(1) == 0);
        CHECK(spc.numMolecules<Space::Tgroup::INACTIVE>(1) == 1);
        CHECK(spc.numMolecules<Space::Tgroup::INACTIVE | Space::Tgroup::NEUTRAL>(1) == 0);

        auto p = getActiveParticles(spc);
        size_t size = 0;
        std::vector<int> vals;
        for (const auto &i : p) {
            size++;
            vals.push_back(int(i.charge));
        }

        CHECK(vals == std::vector<int>({1, 2, 6, 7, 8}));
        CHECK(size == p.size());

        // now let's check the rangev3 implementation
        // in `activeParticles()`:
        auto p2 = spc.activeParticles();
        CHECK(std::distance(p2.begin(), p2.end()) == size);
        vals.clear();
        for (const auto &i : p2) {
            vals.push_back(int(i.charge));
        }
        CHECK(vals == std::vector<int>({1, 2, 6, 7, 8}));
    }

    SUBCASE("SpaceFactory") {
        Space spc;
        SpaceFactory::makeNaCl(spc, 10, R"( {"type": "cuboid", "length": 20} )"_json);
        CHECK(spc.numParticles() == 20);
    }
}

TEST_SUITE_END();
} // namespace Faunus
