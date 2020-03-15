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

TEST_CASE("[Faunus] Gouy-Chapman") {
    Geometry::Slit slit(50, 50, 50);
    Geometry::Chameleon geometry(slit, Geometry::SLIT);
    json j = {{"molarity", 0.1}, {"epsr", 80}, {"linearise", false}, {"rhoinv", 100.0}};
    auto phi = Energy::createGouyChapmanPotential(j, geometry);
    Particle p;
    p.charge = 1.0;
    p.pos = {0, 0, -25};                            // potential at charged surface
    CHECK(phi(p) == doctest::Approx(0.2087776151)); // = phi_0

    p.pos = {0, 0, 0}; // potential at mid-plane
    CHECK(phi(p) == doctest::Approx(0.0160227029));

    j = {{"molarity", 0.1}, {"epsr", 80}, {"linearise", false}, {"phi0", 0.2087776151}};
    phi = Energy::createGouyChapmanPotential(j, geometry);
    CHECK(phi(p) == doctest::Approx(0.0160227029));

    j = {{"molarity", 0.1}, {"epsr", 80}, {"linearise", false}, {"rho", 0.01}};
    phi = Energy::createGouyChapmanPotential(j, geometry);
    CHECK(phi(p) == doctest::Approx(0.0160227029));

    j = {{"molarity", 0.1}, {"epsr", 80}, {"linearise", true}, {"rho", 0.01}};
    phi = Energy::createGouyChapmanPotential(j, geometry);
    CHECK(phi(p) == doctest::Approx(0.0160371645));
}
} // namespace Energy
} // namespace Faunus
