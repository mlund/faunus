#include "smart_montecarlo.h"
#include <doctest/doctest.h>

namespace Faunus::SmarterMonteCarlo {

/**
 * @param outside_acceptance Probability to accept if outside region ("p")
 * @param n_total Total number of elements
 * @param n_inside Total number of elements inside region
 * @param direction Did element exit or enter the region? Returns zero if no boundary crossing
 *
 * - See Allen and Tildesley p. 318 (2017 ed.).
 * - Original reference: [doi:10/frvx8j](https://doi.org/frvx8j)
 */
double bias(double outside_acceptance, const int n_total, const int n_inside,
            BiasDirection direction)
{
    const auto p = outside_acceptance;
    const auto n_prime = p * n_total + (1.0 - p) * n_inside;
    switch (direction) {
    case BiasDirection::EXIT_REGION: // in --> out
        return -std::log(p / ((1.0 - (1.0 - p) / n_prime)));
    case BiasDirection::ENTER_REGION: // out --> in
        return std::log(p * (1.0 + (1.0 - p) / n_prime));
    default:
        return 0.0;
    }
}

TEST_CASE("[Faunus] SmartMonteCarlo::bias")
{
    using doctest::Approx;
    CHECK_EQ(bias(1.0, 20, 5, BiasDirection::NO_CROSSING), Approx(0.0));
    CHECK_EQ(bias(1.0, 20, 5, BiasDirection::EXIT_REGION), Approx(0.0));
    CHECK_EQ(bias(1.0, 20, 5, BiasDirection::ENTER_REGION), Approx(0.0));

    CHECK_EQ(bias(0.1, 20, 5, BiasDirection::NO_CROSSING), Approx(0.0));
    CHECK_EQ(bias(0.1, 20, 5, BiasDirection::EXIT_REGION), Approx(2.1535495138));
    CHECK_EQ(bias(0.1, 20, 5, BiasDirection::ENTER_REGION), Approx(-2.1729072697));
}

/**
 * @param outside_acceptance Probability to accept element outside the region
 * @param region Region to preferentially pick from
 */
RegionSampler::RegionSampler(const double outside_acceptance,
                             std::unique_ptr<Region::RegionBase> region)
    : outside_acceptance(outside_acceptance)
    , region(std::move(region))
{
    if (outside_acceptance <= pc::epsilon_dbl || outside_acceptance > 1.0) {
        throw ConfigurationError("outside_acceptance (p), must be in the range (0,1]");
    }
}

/** Determines the direction of a transition */
BiasDirection RegionSampler::getDirection(const bool inside_before, const bool inside_after)
{
    if (inside_before && (not inside_after)) {
        return BiasDirection::EXIT_REGION;
    }
    else if (inside_after && (not inside_before)) {
        return BiasDirection::ENTER_REGION;
    }
    return BiasDirection::NO_CROSSING;
}

void RegionSampler::to_json(json& j) const
{
    j["region"] = static_cast<json>(*region);
    j["p"] = outside_acceptance;
    if (fixed_number_inside) {
        j["fixed <N_inside>"] = fixed_number_inside.value();
    }
}

} // namespace Faunus::SmarterMonteCarlo