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
}

TEST_CASE("[Faunus] Space") {
    Tspace spc1;
    spc1.geo = R"( {"type": "sphere", "radius": 1e9} )"_json;

    // check molecule insertion
    molecules.at(0).atomic = false;
    REQUIRE_EQ(molecules.at(0).atomic, false);
    atoms.resize(2);
    CHECK(atoms.at(0).mw == 1);
    Particle a;
    a.id = 0;
    a.pos.setZero();
    Tspace::Tpvec p(2, a);
    CHECK(p[0].traits().mw == 1);
    p[0].pos.x() = 2;
    p[1].pos.x() = 3;
    spc1.push_back(0, p);
    CHECK(spc1.p.size() == 2);
    CHECK(spc1.groups.size() == 1);
    CHECK(spc1.groups.front().id == 0);
    CHECK(spc1.groups.front().cm.x() == doctest::Approx(2.5));

    // check `positions()`
    CHECK(&spc1.positions()[0] == &spc1.p[0].pos);

    // sync groups
    Change c;
    c.all = true;
    c.dV = true;
    c.groups.resize(1);
    c.groups[0].index = 0;
    c.groups[0].all = true;
    Tspace spc2;
    spc2.sync(spc1, c);
    CHECK(spc2.p.size() == 2);
    CHECK(spc2.groups.size() == 1);
    CHECK(spc2.groups.front().id == 0);
    CHECK(spc2.groups.front().begin() != spc1.groups.front().begin());
    CHECK(spc2.p.front().pos.x() == doctest::Approx(2));

    // nothing should be synched (all==false)
    spc2.p.back().pos.z() = -0.1;
    c.all = false;
    c.groups[0].all = false;
    spc1.sync(spc2, c);
    CHECK(spc1.p.back().pos.z() != -0.1);

    // everything should be synched (all==true)
    c.groups[0].all = true;
    spc1.sync(spc2, c);
    CHECK(spc1.p.back().pos.z() == doctest::Approx(-0.1));

    SUBCASE("getActiveParticles") {
        // add three groups to space
        Tspace spc;
        spc.geo = R"( {"type": "sphere", "radius": 1e9} )"_json;
        Particle a;
        a.pos.setZero();
        a.id = 0;
        typename Tspace::Tpvec pvec({a, a, a});

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
