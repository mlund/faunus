#pragma once
#include "energy.h"
#include "core.h"
#include "units.h"


namespace Faunus {

using doctest::Approx;

TEST_SUITE_BEGIN("Energy");
// Mocking is unfortunately not supported by doctest yet, hence the complete Space object has to be constructed.
// Furthermore, construction of multiple independent spaces is not straightforwardly possible because of
// the global variables containing atom and molecule types.

#ifdef ENABLE_FREESASA
TEST_CASE( "[Faunus] FreeSASA") {
    Change change;          // change object telling that a full energy calculation
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
    spc.p[0].pos = {0.0,0.0,0.0};
    spc.p[1].pos = {0.0,0.0,20.0};

    SUBCASE("Separated atoms") {
        Energy::SASAEnergy sasa(spc, 1.5_molar, 1.4_angstrom);
        CHECK(sasa.energy(change) == Approx(4*pc::pi*(3.4*3.4 + 2.6*2.6) * 1.5 * 1.0_kJmol));
    }

    SUBCASE("Intersecting atoms") {
        Energy::SASAEnergy sasa(spc, 1.5_molar, 1.4_angstrom);
        std::vector<double> distance = {0.0, 2.5, 5.0, 7.5, 10.0};
        std::vector<double> sasa_energy = {87.3576, 100.4612, 127.3487, 138.4422, 138.4422};
        for(size_t i = 0; i < distance.size(); ++i) {
            spc.p[1].pos = {0.0,0.0, distance[i]};
            CHECK(sasa.energy(change) == Approx(sasa_energy[i]).epsilon(0.02));
        }

    }

    SUBCASE("PBC") {

    }
}
#endif

TEST_SUITE_END();
} // namespace Faunus
