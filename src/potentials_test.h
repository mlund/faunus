#include "potentials.h"

namespace Faunus {
namespace Potential {

using namespace std::string_literals;
using doctest::Approx;

/*
 * A minimal testcase for potentials derived from the MixerPairPotentialBase class includes
 * 1. evaluation with the explicit combination rule(s) and no custom pairs
 * 2. evaluation with the explicit combination rule and custom pairs
 */

TEST_CASE("[Faunus] LennardJones") {
    atoms = R"([{"A": {"sigma":2, "eps":0.9}},
                 {"B": {"sigma":8, "eps":0.1}}])"_json.get<decltype(atoms)>();
    Particle a, b;
    a = atoms[0];
    b = atoms[1];
    LennardJones lj_lb = R"({"mixing": "LB"})"_json;
    LennardJones lj_geom = R"({"mixing": "geometric"})"_json;
    LennardJones lj_custom = R"({"mixing": "LB", "custom": {"A B": {"eps": 0.5, "sigma": 8}}})"_json;

    double d = 0.9_nm;
    auto lj_func = [d](double sigma, double eps) -> double {
        return 4 * eps * (std::pow(sigma / d, 12) - std::pow(sigma / d, 6));
    };

    CHECK(lj_lb(a, a, {0, 0, d}) == Approx(lj_func(0.2_nm, 0.9_kJmol)));
    CHECK(lj_lb(a, b, {0, 0, d}) == Approx(lj_func(0.5_nm, 0.3_kJmol)));
    CHECK(lj_geom(a, b, {0, 0, d}) == Approx(lj_func(0.4_nm, 0.3_kJmol)));
    CHECK(lj_custom(a, b, {0, 0, d}) == Approx(lj_func(0.8_nm, 0.5_kJmol)));
    CHECK(lj_lb(a, a, {0, 0, d}) == lj_custom(a, a, {0, 0, d}));
    CHECK_THROWS_AS(LennardJones lj_unknown = R"({"mixing": "unknown"})"_json, std::runtime_error);
    // alternative notation for custom as an array: custom: []
    CHECK_NOTHROW(LennardJones lj_custom_alt =
                      R"({"mixing": "LB", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json);
}

TEST_CASE("[Faunus] HardSphere") {
    atoms = R"([{"A": {"sigma": 2}}, {"B": {"sigma": 8}}])"_json.get<decltype(atoms)>();
    Particle a, b;
    a = atoms[0];
    b = atoms[1];
    HardSphere hs = R"({"mixing": "arithmetic"})"_json;
    HardSphere hs_custom = R"({"custom": {"A B": {"sigma": 6}}})"_json;

    CHECK(hs(a, a, {0, 0, 2.1_angstrom}) == 0);
    CHECK(hs(a, a, {0, 0, 1.9_angstrom}) == pc::infty);
    CHECK(hs(a, b, {0, 0, 5.1_angstrom}) == 0);
    CHECK(hs(a, b, {0, 0, 4.9_angstrom}) == pc::infty);
    CHECK(hs_custom(a, b, {0, 0, 6.1_angstrom}) == 0);
    CHECK(hs_custom(a, b, {0, 0, 5.9_angstrom}) == pc::infty);

    CHECK_NOTHROW(HardSphere hs_default = R"({})"_json);
    CHECK_THROWS_AS(HardSphere hs_unknown = R"({"mixing": "unknown"})"_json, std::runtime_error);
}

TEST_CASE("[Faunus] SquareWell") {
    atoms = R"([{"A": { "r": 5, "sigma_sw":4, "eps_sw":0.2 }},
                 {"B": { "r": 10, "sigma_sw":2, "eps_sw":0.1 }} ])"_json.get<decltype(atoms)>();
    Particle a, b;
    a = atoms[0];
    b = atoms[1];
    SquareWell pot = R"({"mixing": "LB"})"_json;

    CHECK(pot(a, b, {0, 0, 5 + 10 + 5.99}) == Approx(-std::sqrt(0.2_kJmol * 0.1_kJmol)));
    CHECK(pot(a, b, {0, 0, 5 + 10 + 6.01}) == Approx(0));
}

TEST_CASE("[Faunus] CustomPairPotential") {
    json j = R"({ "atomlist" : [
                 {"A": { "q":1.0,  "r":3, "eps":0.1 }},
                 {"B": { "q":-1.0, "r":4, "eps":0.05 }} ]})"_json;

    atoms = j["atomlist"].get<decltype(atoms)>();

    Particle a, b;
    a = atoms[0];
    b = atoms[1];

    CustomPairPotential pot = R"({
                "constants": { "kappa": 30, "lB": 7},
                "function": "lB * q1 * q2 / (s1+s2) * exp(-kappa/r) * kT + pi"})"_json;

    CHECK(pot(a, b, {0, 0, 2}) == Approx(-7 / (3.0 + 4.0) * std::exp(-30 / 2) * pc::kT() + pc::pi));
}

TEST_CASE("[Faunus] FunctorPotential") {
    json j = R"({ "atomlist" : [
                 {"A": { "q":1.0,  "r":1.1, "eps":0.1 }},
                 {"B": { "q":-1.0, "r":2.0, "eps":0.05 }},
                 {"C": { "r":1.0 }} ]})"_json;

    atoms = j["atomlist"].get<decltype(atoms)>();

    FunctorPotential u = R"(
                {
                  "default": [
                    { "coulomb" : {"epsr": 80.0, "type": "plain", "cutoff":20} }
                  ],
                  "A B" : [
                    { "coulomb" : {"epsr": 80.0, "type": "plain", "cutoff":20} },
                    { "wca" : {"mixing": "LB"} }
                  ],
                  "C C" : [
                    { "hardsphere" : {} }
                  ]
                 }
                )"_json;

    Coulomb coulomb = R"({ "coulomb": {"epsr": 80.0, "type": "plain", "cutoff":20} } )"_json;
    WeeksChandlerAndersen wca = R"({ "wca" : {"mixing": "LB"} })"_json;

    Particle a = atoms[0];
    Particle b = atoms[1];
    Particle c = atoms[2];
    Point r = {2, 0, 0};
    CHECK(u(a, a, r) == Approx(coulomb(a, a, r)));
    CHECK(u(b, b, r) == Approx(coulomb(b, b, r)));
    CHECK(u(a, b, r) == Approx(coulomb(a, b, r) + wca(a, b, r)));
    CHECK(u(c, c, r * 1.01) == 0);
    CHECK(u(c, c, r * 0.99) == pc::infty);
}

TEST_CASE("[Faunus] Dipole-dipole interactions") {
    json j = R"({ "atomlist" : [
                 {"A": { "mu":[1.0,0.0,0.0], "mulen":3.0 }},
                 {"B": { "mu":[0.0,1.0,0.0], "mulen":3.0 }},
                 {"C": { "mu":[1.0,1.0,0.0] }} ]})"_json;

    atoms = j["atomlist"].get<decltype(atoms)>();

    FunctorPotential u = R"(
                {
                  "default": [
                    { "dipoledipole" : {"epsr": 1.0, "type": "plain", "cutoff":20} }
                  ]
                 }
                )"_json;

    DipoleDipole dipoledipole = R"({ "dipoledipole": {"epsr": 1.0, "type": "plain", "cutoff":20} } )"_json;

    Particle a = atoms[0];
    Particle b = atoms[1];
    Particle c = atoms[2];
    Point r = {2, 0, 0};
    CHECK(u(a, a, r) ==
          Approx(dipoledipole(a, a,
                              r))); // interaction between two parallell dipoles, directed parallell to their seperation
    CHECK(u(b, b, r) ==
          Approx(dipoledipole(
              b, b, r))); // interaction between two parallell dipoles, directed perpendicular to their seperation
    CHECK(u(a, b, r) == Approx(dipoledipole(a, b, r))); // interaction between two perpendicular dipoles
    CHECK(u(a, a, r) == -2.25 * dipoledipole.lB);
    CHECK(u(b, b, r) == 1.125 * dipoledipole.lB);
    CHECK(u(a, c, r) == -0.75 * dipoledipole.lB);
    CHECK(u(b, c, r) == 0.375 * dipoledipole.lB);
    CHECK(u(a, b, r) == 0);

    r = {3, 0, 0};
    CHECK(u(a, a, r) ==
          Approx(dipoledipole(a, a,
                              r))); // interaction between two parallell dipoles, directed parallell to their seperation
    CHECK(u(b, b, r) ==
          Approx(dipoledipole(
              b, b, r))); // interaction between two parallell dipoles, directed perpendicular to their seperation
    CHECK(u(a, b, r) == Approx(dipoledipole(a, b, r))); // interaction between two perpendicular dipoles
    CHECK(u(a, a, r) == -(2.0 / 3.0) * dipoledipole.lB);
    CHECK(u(b, b, r) == (1.0 / 3.0) * dipoledipole.lB);
    CHECK(u(a, c, r) == -2.0 / 9.0 * dipoledipole.lB);
    CHECK(u(b, c, r) == 1.0 / 9.0 * dipoledipole.lB);
    CHECK(u(a, b, r) == 0);
}

TEST_CASE("[Faunus] Pair Potentials") {
    json j = R"({ "atomlist" : [
                 { "A": { "r": 1.5, "tension": 0.023} },
                 { "B": { "r": 2.1, "tfe": 0.98 } }]})"_json;

    atoms = j["atomlist"].get<decltype(atoms)>();
    Particle a, b;
    a.id = 0;
    b.id = 1;

    SUBCASE("SASApotential") {
        SASApotential pot;
        json in = R"({ "sasa": {"molarity": 1.0, "radius": 0.0, "shift":false}})"_json;
        pot = in["sasa"];
        double conc = 1.0 * 1.0_molar;
        double tension = atoms[a.id].tension / 2;
        double tfe = atoms[b.id].tfe / 2;
        double f = tension + conc * tfe;
        CHECK(tension > 0.0);
        CHECK(conc > 0.0);
        CHECK(tfe > 0.0);
        CHECK(f > 0.0);
        CHECK(in == json(pot));
        CHECK(pot(a, b, {0, 0, 0}) == Approx(f * 4 * pc::pi * 2.1 * 2.1));                // complete overlap
        CHECK(pot(a, b, {10, 0, 0}) == Approx(f * 4 * pc::pi * (2.1 * 2.1 + 1.5 * 1.5))); // far apart
        CHECK(pot(a, b, {2.5, 0, 0}) == Approx(f * 71.74894965974514));                   // partial overlap
    }
}

} // namespace Potential
} // namespace Faunus
