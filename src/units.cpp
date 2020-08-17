#include <doctest/doctest.h>
#include "units.h"

double Faunus::PhysicalConstants::temperature = 298.15;

std::string Faunus::u8::bracket(const std::string &s) {
    return "\u27e8" + s + "\u27e9";
}

/**
 * @param concentration Salt concentration (arbitrary units)
 * @param valency Absolute valencies, i.e. 1:1, 2:1, 2:2 salt
 * @return Ionic strength (same unit as `conc`)
 * @note For STL gcd algorithm for >=2 numbers,
 *       see https://github.com/sol-prog/cpp-gcd-example/blob/master/gcd_4.cpp
 */
double Faunus::ionicStrength(const double concentration, const std::array<unsigned int, 2> &valency) {
    if (valency[0] == 0 or valency[1] == 0) {
        throw std::runtime_error("valencies cannot be zero");
    }
    const double gcd = std::gcd(valency[0], valency[1]); // greatest common divider
    const double nu0 = valency[1] / gcd;                 // stoichiometric coefficient
    const double nu1 = valency[0] / gcd;                 // stoichiometric coefficient
    return 0.5 * concentration * (nu0 * valency[0] * valency[0] + nu1 * valency[1] * valency[1]);
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
    CHECK_THROWS(ionicStrength(0.1, {1, 0}));
    CHECK_THROWS(ionicStrength(0.1, {0, 1}));
    CHECK_THROWS(ionicStrength(0.1, {0, 0}));
    SUBCASE("debyeLength") { CHECK(debyeLength(0.03, {1, 1}, 7) == Approx(17.7376102214)); }
}

/**
 * @param molar_ionic_strength Molar ionic strength
 * @param bjerrum_length Bjerrum length in Angstrom
 * @return Debye screening length in Angstrom
 */
double Faunus::debyeLength(const double molar_ionic_strength, const double bjerrum_length) {
    return 1.0 / std::sqrt(8.0 * pc::pi * bjerrum_length * molar_ionic_strength * 1.0_molar);
}

/**
 * @param molarity Molar salt concentration
 * @param valency Absolute valencies, i.e. 1:1, 2:1, 2:2 salt
 * @param bjerrum_length Bjerrum length in Angstrom
 * @return Debye screening length in Angstrom
 */
double Faunus::debyeLength(const double molarity, const std::array<unsigned int, 2> &valency,
                           const double bjerrum_length) {
    const double molar_ionic_strength = ionicStrength(molarity, valency);
    return debyeLength(molar_ionic_strength, bjerrum_length);
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
