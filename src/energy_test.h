#pragma once
#include "energy.h"
#include "core.h"
#include "units.h"

namespace Faunus {
namespace Energy {

using doctest::Approx;

TEST_SUITE_BEGIN("Energy");
// Mocking is unfortunately not supported by doctest yet, hence the complete Space object has to be constructed.
// Furthermore, construction of multiple independent spaces is not straightforwardly possible because of
// the global variables containing atom and molecule types.

TEST_CASE("[Faunus] Ewald - EwaldData") {
    using doctest::Approx;

    Space spc;
    EwaldData data(R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "kcutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json);

    CHECK(data.policy == EwaldData::PBC);
    CHECK(data.const_inf == 1);
    CHECK(data.alpha == 0.894427190999916);

    // Check number of wave-vectors using PBC
    PolicyIonIon ionion(spc);
    ionion.updateBox(data, Point(10, 10, 10));
    CHECK(data.kVectors.cols() == 2975);
    CHECK(data.Qion.size() == data.kVectors.cols());

    // Check number of wave-vectors using IPBC
    data.policy = EwaldData::IPBC;
    PolicyIonIonIPBC ionionIPBC(spc);
    ionionIPBC.updateBox(data, Point(10, 10, 10));
    CHECK(data.kVectors.cols() == 846);
    CHECK(data.Qion.size() == data.kVectors.cols());
}

TEST_CASE("[Faunus] Ewald - IonIonPolicy") {
    using doctest::Approx;
    Space spc;
    spc.p.resize(2);
    spc.geo = R"( {"type": "cuboid", "length": 10} )"_json;
    spc.p[0] = R"( {"pos": [0,0,0], "q": 1.0} )"_json;
    spc.p[1] = R"( {"pos": [1,0,0], "q": -1.0} )"_json;
    Group<Particle> g(spc.p.begin(), spc.p.end());
    spc.groups.push_back(g);

    EwaldData data = R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "kcutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json;
    Change c;
    c.all = true;
    data.policy = EwaldData::PBC;

    SUBCASE("PBC") {
        PolicyIonIon ionion(spc);
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data);
        CHECK(ionion.selfEnergy(data, c) == Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK(ionion.surfaceEnergy(data, c) == Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK(ionion.reciprocalEnergy(data) == Approx(0.21303063979675319 * data.bjerrum_length));
    }

    SUBCASE("PBCEigen") {
        PolicyIonIonEigen ionion(spc);
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data);
        CHECK(ionion.selfEnergy(data, c) == Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK(ionion.surfaceEnergy(data, c) == Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK(ionion.reciprocalEnergy(data) == Approx(0.21303063979675319 * data.bjerrum_length));
    }

    SUBCASE("IPBC") {
        PolicyIonIonIPBC ionion(spc);
        data.policy = EwaldData::IPBC;
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data);
        CHECK(ionion.selfEnergy(data, c) == Approx(-1.0092530088080642 * data.bjerrum_length));
        CHECK(ionion.surfaceEnergy(data, c) == Approx(0.0020943951023931952 * data.bjerrum_length));
        CHECK(ionion.reciprocalEnergy(data) == Approx(0.0865107467 * data.bjerrum_length));
    }

    // IPBCEigen is under construction
    /*SUBCASE("IPBCEigen") {
        PolicyIonIonIPBCEigen ionion(spc);
        data.type = EwaldData::IPBCEigen;
        ionion.updateBox(data, spc.geo.getLength());
        ionion.updateComplex(data);
        CHECK(ionion.selfEnergy(data, c) == Approx(-1.0092530088080642 * data.lB));
        CHECK(ionion.surfaceEnergy(data, c) == Approx(0.0020943951023931952 * data.lB));
        CHECK(ionion.reciprocalEnergy(data) == Approx(0.0865107467 * data.lB));
    }*/
}

#ifdef ENABLE_FREESASA
TEST_CASE("[Faunus] FreeSASA") {
    Change change; // change object telling that a full energy calculation
    change.all = true;
    pc::temperature = 300.0_K;
    atoms = R"([
        { "A": { "sigma": 4.0, "tfe": 1.0 } },
        { "B": { "sigma": 2.4, "tfe": 1.0 } }
    ])"_json.get<decltype(atoms)>();
    molecules = R"([
        { "M": { "atoms": ["A", "B"], "atomic": true } }
    ])"_json.get<decltype(molecules)>();
    json j = R"({
        "geometry": {"type": "sphere", "radius": 100 },
        "insertmolecules": [ { "M": { "N": 1 } } ]
    })"_json;
    Space spc = j;
    spc.p[0].pos = {0.0, 0.0, 0.0};
    spc.p[1].pos = {0.0, 0.0, 20.0};

    SUBCASE("Separated atoms") {
        SASAEnergy sasa(spc, 1.5_molar, 1.4_angstrom);
        CHECK(sasa.energy(change) == Approx(4 * pc::pi * (3.4 * 3.4 + 2.6 * 2.6) * 1.5 * 1.0_kJmol));
    }

    SUBCASE("Intersecting atoms") {
        SASAEnergy sasa(spc, 1.5_molar, 1.4_angstrom);
        std::vector<double> distance = {0.0, 2.5, 5.0, 7.5, 10.0};
        std::vector<double> sasa_energy = {87.3576, 100.4612, 127.3487, 138.4422, 138.4422};
        for (size_t i = 0; i < distance.size(); ++i) {
            spc.p[1].pos = {0.0, 0.0, distance[i]};
            CHECK(sasa.energy(change) == Approx(sasa_energy[i]).epsilon(0.02));
        }
    }

    SUBCASE("PBC") {}
}
#endif

TEST_SUITE_END();
} // namespace Energy
} // namespace Faunus
