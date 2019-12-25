#include "potentials.h"

namespace Faunus {
namespace Potential {

using namespace std::string_literals;
using doctest::Approx;

TEST_SUITE_BEGIN("MixerPairPotentials");

/*
 * A minimal testcase for potentials derived from the MixerPairPotentialBase class includes
 * 1. evaluation with the explicit combination rule(s) and no custom pairs
 * 2. evaluation with the explicit combination rule and custom pairs
 */

TEST_CASE("[Faunus] PairMixer") {
    SUBCASE("Enumerated potential") {
        REQUIRE(PairMixer::combArithmetic(2.0, 8.0) == Approx(5.0));
        REQUIRE(PairMixer::combGeometric(2.0, 8.0) == Approx(4.0));
        CHECK(PairMixer::getCombinator(COMB_LORENTZ_BERTHELOT, PairMixer::COEF_SIGMA)(2.0, 8.0) ==
              PairMixer::combArithmetic(2.0, 8.0));
        CHECK(PairMixer::getCombinator(COMB_LORENTZ_BERTHELOT, PairMixer::COEF_EPSILON)(2.0, 8.0) ==
              PairMixer::combGeometric(2.0, 8.0));
        CHECK_THROWS_AS(PairMixer::getCombinator(COMB_LORENTZ_BERTHELOT), std::logic_error);

        SUBCASE("") {
            atoms =
                R"([{"A": {"sigma":2.0}}, {"B": {"sigma":8.0}}, {"C": {"sigma":18.0}}])"_json.get<decltype(atoms)>();
            REQUIRE(atoms.front().interaction.get("sigma") == Approx(2.0));
            std::vector<CustomInteractionData> pairs = R"([{"A C": {"sigma": 9.5}}, {"C B": {"sigma": 12.5}}])"_json;
            TExtractorFunc sigma = [](InteractionData a) -> double { return a.get("sigma"); };

            SUBCASE("") {
                PairMixer mixer(sigma, &PairMixer::combArithmetic);
                SUBCASE("Atom pairs") {
                    auto matrix = mixer.createPairMatrix(atoms);
                    CHECK(matrix->isApprox(matrix->transpose())); // symmetric
                    CHECK((*matrix)(0, 0) == Approx(2.0));
                    CHECK((*matrix)(0, 1) == Approx(5.0));
                }
                SUBCASE("Custom pairs") {
                    auto matrix = mixer.createPairMatrix(atoms, pairs);
                    CHECK(matrix->isApprox(matrix->transpose())); // symmetric
                    CHECK((*matrix)(0, 0) == Approx(2.0));
                    CHECK((*matrix)(0, 1) == Approx(5.0));
                    CHECK((*matrix)(2, 0) == Approx(9.5));
                    CHECK((*matrix)(2, 1) == Approx(12.5));
                }
            }
            SUBCASE("Modifier") {
                PairMixer mixer(sigma, &PairMixer::combArithmetic, [](double x) { return 10 * x; });
                auto matrix = mixer.createPairMatrix(atoms, pairs);
                CHECK(matrix->isApprox(matrix->transpose())); // symmetric
                CHECK((*matrix)(0, 0) == Approx(20.0));
                CHECK((*matrix)(0, 1) == Approx(50.0));
                CHECK((*matrix)(2, 0) == Approx(95.0));
                CHECK((*matrix)(2, 1) == Approx(125.0));
            }
            SUBCASE("Alternative JSON") {
                CHECK_NOTHROW(R"({"A C": {"sigma": 9.5}})"_json.get<std::vector<CustomInteractionData>>());
                std::vector<CustomInteractionData> alt_pairs = R"({"A C": {"sigma": 9.5}, "C B": {"sigma": 12.5}})"_json;
                CHECK_EQ(alt_pairs.size(), pairs.size());
            }
        }
    }
}

TEST_CASE("[Faunus] LennardJones") {
    atoms = R"([{"A": {"sigma":2.0, "eps":0.9}},
                {"B": {"sigma":8.0, "eps":0.1}},
                {"C": {"sigma":5.0, "eps":1.1}}])"_json.get<decltype(atoms)>();
    Particle a = atoms[0], b = atoms[1];

    double d = 0.9_nm;
    auto lj_func = [d](double sigma, double eps) -> double {
        return 4 * eps * (std::pow(sigma / d, 12) - std::pow(sigma / d, 6));
    };

    SUBCASE("JSON initilization") {
        CHECK_THROWS_AS(LennardJones lj = R"({"mixing": "unknown"})"_json, std::runtime_error);
        CHECK_NOTHROW(LennardJones lj = R"({})"_json);
        // alternative notation for custom as an array: custom: []
        CHECK_NOTHROW(LennardJones lj = R"({"mixing": "LB", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json);
    }
    SUBCASE("JSON output custom") {
        json j_in = R"({"mixing": "LB", "sigma": "sigma", "eps": "eps",
            "custom": [{"A B": {"eps": 0.5, "sigma": 8}}, {"A C": {"eps": 1.0, "sigma": 5}}]})"_json;
        LennardJones lj = j_in;
        json j_out = lj;
        // The custom pair data in the output contain a lot of ballast among the original data.
        // If the original data match the output data, patching of the output shall not change it.
        json j_custom = j_out["lennardjones"]["custom"];
        j_custom[0].merge_patch(j_in["custom"][0]);
        j_custom[1].merge_patch(j_in["custom"][1]);
        CHECK_EQ(j_out["lennardjones"]["custom"], j_custom);
    }
    SUBCASE("Lorentz-Berthelot mixing") {
        LennardJones lj = R"({"mixing": "LB"})"_json;
        CHECK(lj(a, a, d * d, {0, 0, d}) == Approx(lj_func(0.2_nm, 0.9_kJmol)));
        CHECK(lj(a, b, d * d, {0, 0, d}) == Approx(lj_func(0.5_nm, 0.3_kJmol)));
    }
    SUBCASE("Geometric mixing") {
        LennardJones lj = R"({"mixing": "geometric"})"_json;
        CHECK(lj(a, a, d * d, {0, 0, d}) == Approx(lj_func(0.2_nm, 0.9_kJmol)));
        CHECK(lj(a, b, d * d, {0, 0, d}) == Approx(lj_func(0.4_nm, 0.3_kJmol)));
    }
    SUBCASE("Custom pairs") {
        LennardJones lj = R"({"mixing": "LB", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json;
        CHECK(lj(a, b, d * d, {0, 0, d}) == Approx(lj_func(0.8_nm, 0.5_kJmol)));
        CHECK(lj(a, a, d * d, {0, 0, d}) == Approx(lj_func(0.2_nm, 0.9_kJmol)));
    }
}

TEST_CASE("[Faunus] WeeksChandlerAndersen") {
    SUBCASE("JSON initilization") {
        atoms = R"([{"A": {"sigma":2.0, "eps":0.9}},
                    {"B": {"sigma":8.0, "eps":0.1}}])"_json.get<decltype(atoms)>();
        Particle a = atoms[0], b = atoms[1];

        CHECK_THROWS_AS(WeeksChandlerAndersen wca = R"({"mixing": "unknown"})"_json, std::runtime_error);
        CHECK_NOTHROW(WeeksChandlerAndersen wca = R"({})"_json);
        CHECK_NOTHROW(WeeksChandlerAndersen wca =
                          R"({"mixing": "LB", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json);

        SUBCASE("Missing coefficient") {
            WeeksChandlerAndersen wca = R"({"mixing": "LB", "sigma": "sigma_wca"})"_json;
            CHECK_EQ(std::isnan(wca(a, a, 10.0 * 10.0, {0, 0, 10.0})), true);
            CHECK_EQ(std::isnan(wca(a, a, 1.0 * 1.0, {0, 0, 1.0})), true);
        }
    }
    SUBCASE("JSON initilization custom coefficient names") {
        atoms = R"([{"A": {"sigma_wca":2.0, "eps":0.9}},
                    {"B": {"sigma_wca":8.0, "eps":0.1}}])"_json.get<decltype(atoms)>();
        Particle a = atoms[0], b = atoms[1];

        CHECK_NOTHROW(WeeksChandlerAndersen wca = R"({"mixing": "LB", "sigma": "sigma_wca"})"_json);
        // Shall throw after non-default potentials are created properly,
        // i.e., not needed pairs are not evaluated at all for the matrices
        // CHECK_THROWS_AS_MESSAGE(WeeksChandlerAndersen wca = R"({"mixing": "LB"})"_json, std::runtime_error,
        //                         "unknown atom property");
        // CHECK_THROWS_AS_MESSAGE(WeeksChandlerAndersen wca = R"({"mixing": "LB", "sigma": "unknown"})"_json,
        //                         std::runtime_error, "unknown atom property");
        // different atom and custom coefficient names are not allowed
        // CHECK_THROWS_AS_MESSAGE(
        //     WeeksChandlerAndersen wca =
        //         R"({"mixing": "LB", "sigma": "sigma_wca", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json,
        //     std::runtime_error, "unknown atom property");
    }
    SUBCASE("JSON serialization") {
        WeeksChandlerAndersen wca;
        wca = R"({"mixing": "LB", "sigma": "sigma_wca"})"_json;
        json j = wca;
        json &j_wca(j["wca"]);
        CHECK_EQ(j_wca["mixing"], "lorentz_berthelot");
        CHECK_EQ(j_wca["sigma"], "sigma_wca");
        CHECK_EQ(j_wca["eps"], "eps");
    }
}

TEST_CASE("[Faunus] HardSphere") {
    atoms = R"([{"A": {"sigma": 2}}, {"B": {"sigma": 8}}])"_json.get<decltype(atoms)>();
    Particle a = atoms[0], b = atoms[1];

    SUBCASE("JSON initialization") {
        CHECK_NOTHROW(HardSphere hs_default = R"({})"_json);
        CHECK_THROWS_AS(HardSphere hs_unknown = R"({"mixing": "unknown"})"_json, std::runtime_error);
    }
    SUBCASE("Undefined mixing") {
        HardSphere hs = R"({"mixing": "undefined"})"_json;
        CHECK(hs(a, a, 1.99 * 1.99, {0, 0, 1.99}) == pc::infty);
        CHECK(hs(a, a, 2.01 * 2.01, {0, 0, 2.01}) == 0.0);
        // CHECK(std::isnan(hs(a, b, {0, 0, 4.99}))); // fails
        // CHECK(std::isnan(hs(a, b, {0, 0, 5.01}))); // fails
    }
    SUBCASE("Arithmetic mixing") {
        HardSphere hs = R"({"mixing": "arithmetic"})"_json;
        CHECK(hs(a, a, 2.01_angstrom * 2.01_angstrom, {0, 0, 2.01_angstrom}) == 0);
        CHECK(hs(a, a, 1.99_angstrom * 1.99_angstrom, {0, 0, 1.99_angstrom}) == pc::infty);
        CHECK(hs(a, b, 5.01_angstrom * 5.01_angstrom, {0, 0, 5.01_angstrom}) == 0);
        CHECK(hs(a, b, 4.99_angstrom * 4.99_angstrom, {0, 0, 4.99_angstrom}) == pc::infty);
    }
    SUBCASE("Custom pairs with implicit mixing") {
        HardSphere hs = R"({"custom": [{"A B": {"sigma": 6}}]})"_json;
        CHECK(hs(a, a, 2.01_angstrom * 2.01_angstrom, {0, 0, 2.01_angstrom}) == 0);
        CHECK(hs(a, a, 1.99_angstrom * 1.99_angstrom, {0, 0, 1.99_angstrom}) == pc::infty);
        CHECK(hs(a, b, 6.01_angstrom * 6.01_angstrom, {0, 0, 6.01_angstrom}) == 0);
        CHECK(hs(a, b, 5.99_angstrom * 5.99_angstrom, {0, 0, 5.99_angstrom}) == pc::infty);
    }
}

TEST_CASE("[Faunus] SquareWell") {
    atoms = R"([{"A": { "r":5,  "sigma":4, "eps":0.2 }},
                {"B": { "r":10, "sigma":2, "eps":0.1 }} ])"_json.get<decltype(atoms)>();
    Particle a = atoms[0], b = atoms[1];

    SUBCASE("JSON initilization") {
        CHECK_THROWS_AS(SquareWell sw = R"({"mixing": "unknown"})"_json, std::runtime_error);
        CHECK_NOTHROW(SquareWell sw = R"({})"_json);
    }
    SUBCASE("Undefined mixing") {
        SquareWell sw = R"({"mixing": "undefined"})"_json;
        CHECK(sw(a, a, 3.99 * 3.99, {0, 0, 0}) == Approx(-0.2_kJmol));
        CHECK(sw(a, a, 4.01 * 4.01, {0, 0, 0}) == Approx(0.0));
        // CHECK(std::isnan(sw(a, b, {0, 0, 5.99}))); // fails
        // CHECK(std::isnan(sw(a, b, {0, 0, 6.01}))); // fails
    }
    SUBCASE("Lorentz-Berthelot mixing") {
        SquareWell sw = R"({"mixing": "LB"})"_json;
        CHECK(sw(a, b, 2.99 * 2.99, {0, 0, 0}) == Approx(-std::sqrt(0.2_kJmol * 0.1_kJmol)));
        CHECK(sw(a, b, 3.01 * 3.01, {0, 0, 0}) == Approx(0));
    }
}

TEST_CASE("[Faunus] Hertz") {

    json j = R"({ "atomlist" : [
                 { "A": { "eps": 1.0, "sigma": 1.3} },
                 { "B": { "eps": 2.0, "sigma": 1.0 } }]})"_json;

    atoms = j["atomlist"].get<decltype(atoms)>();
    Particle a = atoms[0], b = atoms[1];

    SUBCASE("JSON serialization") {
        json hertz_json = R"({ "hertz": {"mixing": "lorentz_berthelot", "eps": "eps", "sigma": "sigma"}})"_json;
        Hertz hertz = hertz_json;
        CHECK(hertz_json == json(hertz));
    }
    SUBCASE("Lorentz-Berthelot mixing") {
        pc::temperature = 298.15_K;
        Hertz hertz = R"({"mixing": "lorentz_berthelot"})"_json;
        CHECK(hertz(a, b, 0.7 * 0.7, {0.7, 0, 0}) == Approx(0.0546424449)); // within cut-off
        CHECK(hertz(a, b, 1.15 * 1.15, {1.15, 0, 0}) == Approx(0.0));       // at cut-off
        CHECK(hertz(a, b, 2.0 * 2.0, {2.0, 0, 0}) == Approx(0.0));          // outside of cut-off
    }
}

TEST_SUITE_END();

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

    CHECK(pot(a, b, 2 * 2, {0, 0, 2}) == Approx(-7 / (3.0 + 4.0) * std::exp(-30 / 2) * pc::kT() + pc::pi));
}

TEST_CASE("[Faunus] FunctorPotential") {
    json j = R"({ "atomlist" : [
                 {"A": { "q":1.0,  "r":1.1, "eps":0.1 }},
                 {"B": { "q":-1.0, "r":2.0, "eps":0.05 }},
                 {"C": { "r":1.0, "mu":[2,0,0] }} ]})"_json;

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
    double r2 = r.squaredNorm();
    CHECK(u(a, a, r2, r) == Approx(coulomb(a, a, r2, r)));
    CHECK(u(b, b, r2, r) == Approx(coulomb(b, b, r2, r)));
    CHECK(u(a, b, r2, r) == Approx(coulomb(a, b, r2, r) + wca(a, b, r2, r)));
    CHECK(u(c, c, (r * 1.01).squaredNorm(), r * 1.01) == 0);
    CHECK(u(c, c, (r * 0.99).squaredNorm(), r * 0.99) == pc::infty);

    SUBCASE("selfEnergy()") {
        // let's check that the self energy gets properly transferred to the functor potential
        json j = R"(
                {"default": [{ "coulomb" : {"epsr": 80.0, "type": "qpotential", "cutoff":20, "order":4} }]})"_json;

        FunctorPotential functor = j;

        NewCoulombGalore galore = j["default"][0];
        CHECK(functor.selfEnergy(a) == Approx(galore.selfEnergy(a)));

        // now test w. dipolar particles
        j = R"(
                {"default": [{ "multipole" : {"epsr": 80.0, "type": "qpotential", "cutoff":20, "order":4} }]})"_json;

        functor = j;
        Multipole multipole = j["default"][0];
        CHECK(functor.selfEnergy(a) == Approx(multipole.selfEnergy(a))); // q=1, mu=0
        CHECK(functor.selfEnergy(c) == Approx(multipole.selfEnergy(c))); // q=0, mu=2
    }
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
                    { "multipole" : {"epsr": 1.0, "type": "plain", "cutoff":20} }
                  ]
                 }
                )"_json;

    Multipole dipoledipole = R"({"epsr": 1.0, "type": "plain", "cutoff":20})"_json;

    Particle a = atoms[0];
    Particle b = atoms[1];
    Particle c = atoms[2];
    Point r = {2, 0, 0};
    double r2 = r.squaredNorm();
    CHECK(u(a, a, r2, r) ==
          Approx(dipoledipole(a, a, r2,
                              r))); // interaction between two parallell dipoles, directed parallell to their seperation
    CHECK(u(b, b, r2, r) ==
          Approx(dipoledipole(
              b, b, r2, r))); // interaction between two parallell dipoles, directed perpendicular to their seperation
    CHECK(u(a, b, r2, r) == Approx(dipoledipole(a, b, r2, r))); // interaction between two perpendicular dipoles
    CHECK(u(a, a, r2, r) == -2.25 * dipoledipole.lB);
    CHECK(u(b, b, r2, r) == 1.125 * dipoledipole.lB);
    CHECK(u(a, c, r2, r) == -0.75 * dipoledipole.lB);
    CHECK(u(b, c, r2, r) == 0.375 * dipoledipole.lB);
    CHECK(u(a, b, r2, r) == 0);

    r = {3, 0, 0};
    r2 = 3 * 3;
    CHECK(u(a, a, r2, r) ==
          Approx(dipoledipole(a, a, r2,
                              r))); // interaction between two parallell dipoles, directed parallell to their seperation
    CHECK(u(b, b, r2, r) ==
          Approx(dipoledipole(
              b, b, r2, r))); // interaction between two parallell dipoles, directed perpendicular to their seperation
    CHECK(u(a, b, r2, r) == Approx(dipoledipole(a, b, r2, r))); // interaction between two perpendicular dipoles
    CHECK(u(a, a, r2, r) == -(2.0 / 3.0) * dipoledipole.lB);
    CHECK(u(b, b, r2, r) == (1.0 / 3.0) * dipoledipole.lB);
    CHECK(u(a, c, r2, r) == -2.0 / 9.0 * dipoledipole.lB);
    CHECK(u(b, c, r2, r) == 1.0 / 9.0 * dipoledipole.lB);
    CHECK(u(a, b, r2, r) == 0);
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
        CHECK(pot(a, b, 0, {0, 0, 0}) == Approx(f * 4 * pc::pi * 2.1 * 2.1));                      // complete overlap
        CHECK(pot(a, b, 10 * 10, {10, 0, 0}) == Approx(f * 4 * pc::pi * (2.1 * 2.1 + 1.5 * 1.5))); // far apart
        CHECK(pot(a, b, 2.5 * 2.5, {2.5, 0, 0}) == Approx(f * 71.74894965974514));                 // partial overlap
    }
}

} // namespace Potential
} // namespace Faunus
