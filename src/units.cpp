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
