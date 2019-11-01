#include "externalpotential.h"

namespace Faunus {
namespace Energy {

using namespace std::string_literals;
using doctest::Approx;

TEST_CASE("[Faunus] ExternalPotential") {
    Faunus::atoms = R"([
        { "A": { "sigma": 4.0, "tfe": 1.0 } },
        { "B": { "sigma": 2.4, "tfe": 1.0 } }
    ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([
        { "M": { "atoms": ["A", "B"], "atomic": true } }
    ])"_json.get<decltype(molecules)>();

    json j = R"({
        "geometry": {"type": "sphere", "radius": 100 },
        "insertmolecules": [ { "M": { "N": 1 } } ]
    })"_json;

    SUBCASE("ParticleSelfEnergy") {
        Space spc = j;
        ParticleSelfEnergy pot(spc, [](const Particle &) { return 0.5; });
        Change change;
        change.all = true; // if both particles have changed
        CHECK(pot.energy(change) == Approx(0.5 + 0.5));
    }
}
} // namespace Energy
} // namespace Faunus
