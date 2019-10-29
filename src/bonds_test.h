#include "bonds.h"
#include "units.h"

namespace Faunus {
namespace Potential {

using doctest::Approx;

TEST_SUITE_BEGIN("Bonds");

TEST_CASE("[Faunus] BondData") {
    ParticleVector p_4a(2, Particle());
    p_4a[0].pos = {2.0, 1.0, -2.0};
    p_4a[1].pos = {2.0, 1.0, 2.0};
    ParticleVector p_60deg_4a(3, Particle());
    p_60deg_4a[0].pos = {1.0, 1.0 + std::sqrt(3), 4.0};
    p_60deg_4a[1].pos = {1.0, 1.0, 1.0};
    p_60deg_4a[2].pos = {1.0, 5.0, 1.0};

    Geometry::DistanceFunction distance = [](const Point &a, const Point &b) -> Point { return b - a; };
    Geometry::DistanceFunction distance_3a = [](const Point &, const Point &) -> Point { return {0, 3, 0}; };
    Geometry::DistanceFunction distance_5a = [](const Point &, const Point &) -> Point { return {0, 3, 4}; };

    typedef std::shared_ptr<BondData> BondDataPtr;
    BondDataPtr bond_ptr;

    SUBCASE("HarmonicBond") {
        SUBCASE("HarmonicBond Energy") {
            HarmonicBond bond(100.0, 5.0, {0, 1});
            bond.setEnergyFunction(p_4a);
            CHECK_EQ(bond.energy(distance_5a), Approx(0));
            CHECK_EQ(bond.energy(distance_3a), Approx(200));
            CHECK_EQ(bond.energy(distance), Approx(50));
        }
        SUBCASE("HarmonicBond JSON") {
            json j = R"({"harmonic": {"index":[1,2], "k":10.0, "req":2.0}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<HarmonicBond>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energy(distance), Approx(10.0_kJmol / 2 * 4));
        }
        SUBCASE("HarmonicBond JSON Invalid") {
            CHECK_NOTHROW((R"({"harmonic": {"index":[0,9], "k":0.5, "req":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"harmoNIC": {"index":[2,3], "k":0.5, "req":2.1}})"_json).get<BondDataPtr>()); // exact match required
            CHECK_THROWS((R"({"harmonic": {"index":[2], "k":0.5, "req":2.1}})"_json).get<BondDataPtr>());   // 2 atom indices
            CHECK_THROWS((R"({"harmonic": {"index":[2,3], "req":2.1}})"_json).get<BondDataPtr>()); // k missing
            CHECK_THROWS((R"({"harmonic": {"index":[2,3], "k":2.1}})"_json).get<BondDataPtr>());   // req missing
        }
    }

    SUBCASE("FENEBond") {
        SUBCASE("FENEBond Energy") {
            FENEBond bond(100.0, 5.0, {0, 1});
            bond.setEnergyFunction(p_4a);
            CHECK_EQ(bond.energy(distance_5a), pc::infty);
            CHECK_EQ(bond.energy(distance_3a), Approx(557.86));
            CHECK_EQ(bond.energy(distance), Approx(1277.06));
        }
        SUBCASE("FENEBond JSON") {
            json j = R"({"fene": {"index":[1,2], "k":8, "rmax":6.0 }})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<FENEBond>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energy(distance), Approx(84.641_kJmol));
        }
        SUBCASE("FENEBond JSON Invalid") {
            CHECK_NOTHROW((R"({"fene": {"index":[0,9], "k":0.5, "rmax":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"FENE": {"index":[0,9], "k":0.5, "rmax":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene": {"index":[2,3,4], "k":1, "rmax":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene": {"index":[2,3], "rmax":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene": {"index":[2,3], "k":1}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("FENEWCABond") {
        SUBCASE("FENEWCABond Energy") {
            FENEWCABond bond(100.0, 5.0, 20.0, 3.2, {0, 1});
            bond.setEnergyFunction(p_4a);
            CHECK_EQ(bond.energy(distance_5a), pc::infty);
            CHECK_EQ(bond.energy(distance_3a), Approx(557.86 + 18.931));
            CHECK_EQ(bond.energy(distance), Approx(1277.06));
        }
        SUBCASE("FENEWCABond JSON") {
            json j = R"({"fene+wca": {"index":[1,2], "k":8, "rmax":6.0, "eps":3.5, "sigma":4.5}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<FENEWCABond>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energy(distance), Approx(92.805_kJmol));
        }
        SUBCASE("FENEWCABond JSON Invalid") {
            CHECK_NOTHROW((R"({"fene+wca": {"index":[0,9], "k":1, "rmax":2.1, "eps":2.48, "sigma":2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3,4], "k":1, "rmax":2.1, "eps":2.48, "sigma":2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3], "rmax":2.1, "eps":2.48, "sigma":2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3], "k":1, "eps":2.48, "sigma":2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3], "k":1, "rmax":2.1, "eps":2.48}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3], "k":1, "rmax":2.1, "sigma":2}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("HarmonicTorsion") {
        SUBCASE("HarmonicTorsion Energy") {
            HarmonicTorsion bond(100.0, 45.0_deg, {0, 1, 2});
            bond.setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond.energy(distance), Approx(100.0 / 2 * std::pow(15.0_deg, 2)));
        }
        SUBCASE("HarmonicTorsion JSON") {
            json j = R"({"harmonic_torsion": {"index":[0,1,2], "k":0.5, "aeq":65}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<HarmonicTorsion>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energy(distance), Approx(0.5_kJmol / 2 * std::pow(5.0_deg, 2)));
        }
        SUBCASE("HarmonicTorsion JSON Invalid") {
            CHECK_NOTHROW((R"({"harmonic_torsion": {"index":[0,1,9], "k":0.5, "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"harmonic_torsion": {"index":[2], "k":0.5, "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"harmonic_torsion": {"index":[0,1,2], "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"harmonic_torsion": {"index":[0,1,3], "k":0.5}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("GromosTorsion") {
        SUBCASE("GromosTorsion Energy") {
            GromosTorsion bond(100.0, cos(45.0_deg), {0, 1, 2});
            bond.setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond.energy(distance), Approx(100.0 / 2 * std::pow(cos(60.0_deg) - cos(45.0_deg), 2)));
        }
        SUBCASE("GromosTorsion JSON") {
            json j = R"({"gromos_torsion": {"index":[0,1,2], "k":0.5, "aeq":65}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<GromosTorsion>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energy(distance), Approx(0.5_kJmol / 2 * std::pow(cos(60.0_deg) - cos(65.0_deg), 2)));
        }
        SUBCASE("GromosTorsion JSON Invalid") {
            CHECK_NOTHROW((R"({"gromos_torsion": {"index":[0,1,9], "k":0.5, "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"gromos_torsion": {"index":[2], "k":0.5, "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"gromos_torsion": {"index":[0,1,2], "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"gromos_torsion": {"index":[0,1,3], "k":0.5}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("PeriodicDihedral") {
        ParticleVector p_45deg(4, Particle());
        p_45deg[1].pos = {0.0, 0.0, 0.0};
        p_45deg[2].pos = {0.0, 0.0, 2.0};
        p_45deg[0].pos = {5.0, 0.0, 0.0};
        p_45deg[3].pos = {10.0, 10.0, 2.0};

        ParticleVector p_90deg(p_45deg);
        p_90deg[3].pos[0] *= 0;
        ParticleVector p_60deg(p_45deg);
        p_60deg[3].pos[1] *= std::sqrt(3);
        ParticleVector p_120deg(p_60deg);
        p_120deg[3].pos[0] *= -1;

        SUBCASE("PeriodicDihedral Energy") {
            PeriodicDihedral bond(100.0, 0.0_deg, 3, {0, 1, 2, 3});
            bond.setEnergyFunction(p_120deg);
            CHECK_EQ(bond.energy(distance), Approx(200.0));
            bond.setEnergyFunction(p_60deg);
            CHECK_EQ(bond.energy(distance), Approx(0.0));
            bond.setEnergyFunction(p_90deg);
            CHECK_EQ(bond.energy(distance), Approx(100.0));
        }
        SUBCASE("PeriodicDihedral JSON") {
            json j = R"({"periodic_dihedral": {"index":[0,1,2,3], "k":10, "phi":0.0, "n": 3}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<PeriodicDihedral>(bond_ptr)->setEnergyFunction(p_90deg);
            CHECK_EQ(bond_ptr->energy(distance), Approx(10.0_kJmol));
        }
        SUBCASE("PeriodicDihedral JSON Invalid") {
            CHECK_NOTHROW((R"({"periodic_dihedral": {"index":[0,1,2,9], "k":0.5, "phi":2.1, "n": 2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"periodic_dihedral": {"index":[0,1,2], "k":0.5, "phi":2.1, "n": 2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"periodic_dihedral": {"index":[0,1,2,3], "phi":2.1, "n": 2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"periodic_dihedral": {"index":[0,1,2,3], "k":0.5, "n": 2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"periodic_dihedral": {"index":[0,1,2,3], "k":0.5, "phi":2.1}})"_json).get<BondDataPtr>());
        }
    }

    // test bond filter
    SUBCASE("filterBonds()") {
        std::vector<std::shared_ptr<BondData>> bonds = {
            R"({"fene":      {"index":[2,3], "k":1, "rmax":2.1, "eps":2.48}} )"_json,
            R"({"harmonic" : {"index":[2,3], "k":0.5, "req":2.1} } )"_json};
        auto filt = filterBonds(bonds, BondData::HARMONIC);
        CHECK(filt.size() == 1);
        CHECK(filt[0]->type() == BondData::HARMONIC);
        CHECK(filt[0] == bonds[1]); // filt should contain references to bonds
    }
}
TEST_SUITE_END();
} // namespace Potential
} // namespace Faunus
