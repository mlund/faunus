#include <doctest/doctest.h>
#include "units.h"
#include <numeric>
#include <cmath>
#include <spdlog/fmt/fmt.h>

double Faunus::PhysicalConstants::temperature = 298.15;

std::string Faunus::u8::bracket(std::string_view sv) {
    return fmt::format("\u27e8{}\u27e9", sv);
}

TEST_CASE("[Faunus] infinite math") {
    using namespace Faunus;
    CHECK(std::isfinite(1.0 + pc::infty) == false);
    CHECK(std::isinf(1.0 + pc::infty));
    CHECK(std::isinf(1.0 + pc::neg_infty));
    CHECK(std::isinf(1.0 / 0.0));
    CHECK(std::isnan(pc::infty / pc::infty));
    CHECK(1.0 / pc::infty == doctest::Approx(0.0));
    CHECK(std::isfinite(0.0 / 0.0) == false);
    CHECK(std::isnan(0.0 / 0.0));
    CHECK(std::isinf(std::exp(pc::infty)));
    CHECK(std::exp(pc::neg_infty) == doctest::Approx(0.0));
    CHECK(std::isinf(std::log(pc::infty)));
    CHECK(std::log(0.0) == pc::neg_infty);
    CHECK(std::log(-0.0) == pc::neg_infty);
    CHECK(std::isnan(std::log(-1.0)));
    CHECK(std::isnan(std::sqrt(-1.0))); // is this needed?
}

/**
 * @param concentration Salt concentration (arbitrary units)
 * @param valency Absolute valencies, i.e. 1:1, 2:1, 2:2 salt
 * @return Ionic strength (same unit as `conc`)
 * @note For STL gcd algorithm for >=2 numbers,
 *       see https://github.com/sol-prog/cpp-gcd-example/blob/master/gcd_4.cpp
 */
double Faunus::ionicStrength(double concentration, const std::array<int, 2> &valency) {
    double gcd = std::gcd(valency[0], valency[1]); // greatest common divider
    double mu0 = valency[1] / gcd;                 // stoichiometric coefficient
    double mu1 = valency[0] / gcd;                 // stoichiometric coefficient
    return 0.5 * concentration * (mu0 * valency[0] * valency[0] + mu1 * valency[1] * valency[1]);
}
TEST_CASE("[Faunus] ionicStrength") {
    using namespace Faunus;
    using doctest::Approx;
    CHECK(ionicStrength(0.1, {1, 1}) == Approx(0.1));
    CHECK(ionicStrength(0.1, {2, 2}) == Approx(0.5 * (0.1 * 4 + 0.1 * 4)));
    CHECK(ionicStrength(0.1, {2, 1}) == Approx(0.5 * (0.1 * 4 + 0.2)));
    CHECK(ionicStrength(0.1, {1, 2}) == Approx(0.5 * (0.2 + 0.1 * 4)));
    CHECK(ionicStrength(0.1, {1, 3}) == Approx(0.5 * (0.3 + 0.1 * 9)));
    CHECK(ionicStrength(0.1, {3, 1}) == Approx(0.5 * (0.3 + 0.1 * 9)));
    CHECK(ionicStrength(0.1, {1, 3}) == Approx(0.5 * (0.3 + 0.1 * 9)));
    CHECK(ionicStrength(0.1, {2, 3}) == Approx(0.5 * (0.3 * 4 + 0.2 * 9)));
    SUBCASE("debyeLength") { CHECK(debyeLength(0.03, {1, 1}, 7) == Approx(17.7376102214)); }
}

/**
 * @param molarity Molar salt concentration
 * @param valency Absolute valencies, i.e. 1:1, 2:1, 2:2 salt
 * @param bjerrum_length Bjerrum length in Angstrom
 * @return Debye screening length in Angstrom
 */
double Faunus::debyeLength(double molarity, const std::array<int, 2> &valency, double bjerrum_length) {
    double ionic_strength = ionicStrength(molarity, valency) * 1.0_molar;
    return 1.0 / std::sqrt(8.0 * pc::pi * bjerrum_length * ionic_strength);
}

TEST_CASE("[Faunus] Units and string literals") {
    using doctest::Approx;
    using namespace Faunus;
    pc::temperature = 298.15_K;
    CHECK(1.0e-10_m == 1);
    CHECK((1 / 1.0_debye) == Approx(4.8032));
    CHECK(1.0_debye == Approx(3.33564e-30_Cm));
    CHECK(1.0_debye == Approx(0.20819434_eA));
    CHECK(360.0_deg == Approx(2 * std::acos(-1)));
    CHECK((1.0_mol / 1.0_liter) == Approx(1.0_molar));
    CHECK(1.0_bar == Approx(0.987_atm));
    CHECK(1.0_atm == Approx(101325._Pa));
    CHECK(1.0_kT == Approx(2.47897_kJmol));
    CHECK(1.0_hartree == Approx(2625.499_kJmol));
}
