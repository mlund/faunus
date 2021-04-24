#include <doctest/doctest.h>
#include "units.h"
#include <cmath>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/all.hpp>
#include <nlohmann/json.hpp>

double Faunus::PhysicalConstants::temperature = 298.15;

std::string Faunus::u8::bracket(const std::string &s) {
    return "\u27e8" + s + "\u27e9";
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
 * @param molarity Molar salt concentration
 * @param valencies valencies for participating ions {1,-1} ~ NaCl, {2,-1} ~ MgCl2, {1,3,-2} ~ KAl(SO4)2
 * @throw If stoichiometry cannot be resolved
 */
Faunus::IonicSalt::IonicSalt(const double molarity, const std::vector<int>& valencies)
    : molarity(molarity), valencies(valencies) {
    using namespace ranges; // @todo migrate to cpp20 ranges
    const auto sum_positive = accumulate(valencies | cpp20::views::filter([](auto v) { return v > 0; }), 0);
    const auto sum_negative = accumulate(valencies | cpp20::views::filter([](auto v) { return v < 0; }), 0);
    const auto gcd = std::gcd(sum_positive, sum_negative);
    if (sum_positive == 0 || sum_negative == 0 || gcd == 0) {
        throw std::runtime_error("cannot resolve stoichiometry; did you provide both + and - ions?");
    }
    auto nu_times_squared_valency = valencies | cpp20::views::transform([&](const auto valency) {
                                        const auto nu = (valency > 0 ? -sum_negative : sum_positive) / gcd;
                                        return nu * valency * valency;
                                    });
    molar_ionic_strength = 0.5 * molarity * static_cast<double>(accumulate(nu_times_squared_valency, 0));
}

/**
 * The salt composition is automatically resolved, and the ionic strength, I, is calculated according to:
 *
 *     I = 0.5 * concentration * sum( nu_i * valency_i^2 )
 *
 * where `nu` are the minimum stoichiometric coefficients, deduced by assuming a net-neutral salt.
 */
double Faunus::IonicSalt::ionicStrength() const { return molar_ionic_strength; }

/**
 * @param bjerrum_length Bjerrum length in Angstrom
 * @return Debye screening length in Angstrom
 */
double Faunus::IonicSalt::debyeLength(const double bjerrum_length) const {
    return 1.0 / std::sqrt(8.0 * pc::pi * bjerrum_length * 1.0_angstrom * molar_ionic_strength * 1.0_molar);
}

/**
 * This will search the json object for:
 * 1. `ionicstrength` and if found assume 1:1 salt
 * 2. `molarity` and assume 1:1 salt unless `valencies` are found
 *
 * @throw If `molarity` cannot be found
 */
Faunus::IonicSalt Faunus::makeIonicSalt(const nlohmann::json& j) {
    const auto molarity = j.at("molarity").get<double>();
    std::vector<int> valencies = {1, -1}; // default 1:1 salt
    valencies = j.value("valencies", valencies);
    return IonicSalt(molarity, valencies);
}

TEST_CASE("[Faunus] IonicSalt") {
    using namespace Faunus;
    using doctest::Approx;
    CHECK(IonicSalt(0.1, {1, -1}).ionicStrength() == Approx(0.1));                                        // NaCl
    CHECK(IonicSalt(0.1, {2, -2}).ionicStrength() == Approx(0.5 * (0.1 * 4 + 0.1 * 4)));                  // CaSO4
    CHECK(IonicSalt(0.1, {2, -1}).ionicStrength() == Approx(0.5 * (0.1 * 4 + 0.2)));                      // CaCl2
    CHECK(IonicSalt(0.1, {1, -2}).ionicStrength() == Approx(0.5 * (0.2 + 0.1 * 4)));                      // K2SO4
    CHECK(IonicSalt(0.1, {1, -3}).ionicStrength() == Approx(0.5 * (0.3 + 0.1 * 9)));                      // Na3Cit
    CHECK(IonicSalt(0.1, {3, -1}).ionicStrength() == Approx(0.5 * (0.3 + 0.1 * 9)));                      // LaCl3
    CHECK(IonicSalt(0.1, {2, -3}).ionicStrength() == Approx(0.5 * (0.3 * 4 + 0.2 * 9)));                  // Ca3(PO4)2
    CHECK(IonicSalt(0.1, {1, 3, -2}).ionicStrength() == Approx(0.5 * (0.1 * 1 + 0.1 * 9 + 0.1 * 2 * 4))); // KAl(SO4)2
    CHECK(IonicSalt(1.0, {2, 3, -2}).ionicStrength() == Approx(0.5 * (2 * 4 + 2 * 9 + 5 * 4))); // Ca2Al2(SO4)5
    CHECK_THROWS(IonicSalt(0.1, {1, 1}));
    CHECK_THROWS(IonicSalt(0.1, {-1, -1}));
    CHECK_THROWS(IonicSalt(0.1, {0, 0}));
    SUBCASE("debyeLength") {
        CHECK(IonicSalt(0.03, {1, -1}).debyeLength(7.0) == Approx(17.7376102214));
        CHECK(getDebyeLength(R"({"debyelength": 30.0})"_json, 7.0) == Approx(30.0));
        CHECK(getDebyeLength(R"({"molarity": 0.03, "epsr": 80})"_json, 7.0) == Approx(17.7376102214));
        CHECK_THROWS(getDebyeLength(R"({"molarity": 0.03, "valencies": [1,1]})"_json, 7.0));
    }
}

/**
 * The json object if searched in this order:
 *
 * 1. `debyelength`; if found, the given bjerrum_length is ignored.
 * 2. `salt/molarity` and optionally `salt/valencies` (default 1:-1 ~ NaCl).
 */
double Faunus::getDebyeLength(const nlohmann::json& j, const double bjerrum_length) {
    using namespace std::string_literals;
    try {
        if (auto it = j.find("salt"); it != j.end()) {
            return getDebyeLength(*it, bjerrum_length);
        }
        if (auto it = j.find("debyelength"); it != j.end()) {
            return it->get<double>();
        }
        return makeIonicSalt(j).debyeLength(bjerrum_length);
    } catch (std::exception& e) {
        throw std::runtime_error("ionic strength: "s + e.what());
    }
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