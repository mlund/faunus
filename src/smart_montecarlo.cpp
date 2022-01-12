#include "smart_montecarlo.h"

namespace Faunus::SmarterMonteCarlo {

/**
 * @param outside_rejection_probability Probability to reject if outside region
 * @param n_total Total number of elements
 * @param n_inside Total number of elements inside region
 * @param direction Did element exit or enter the region? Returns zero if no boundary crossing
 *
 * - See Allen and Tildesley p. 318 (2017 ed.).
 * - Original reference: [doi:10/frvx8j](https://doi.org/frvx8j)
 */
double bias(double outside_rejection_probability, const int n_total, const int n_inside, BiasDirection direction) {
    const auto p = outside_rejection_probability;
    const auto n_prime = p * n_total + (1.0 - p) * n_inside;
    switch (direction) {
    case BiasDirection::EXIT_REGION: // in --> out
        return -std::log(p / ((1.0 - (1.0 - p) / n_prime)));
    case BiasDirection::ENTER_REGION: // out --> in
        return std::log(p * (1.0 + (1.0 - p) / n_prime));
    case BiasDirection::NO_CROSSING:
        return 0.0;
    }
}

TEST_CASE("[Faunus] SmartMonteCarlo::bias") {
    using doctest::Approx;
    CHECK(bias(1.0, 20, 5, BiasDirection::NO_CROSSING) == Approx(0.0));
    CHECK(bias(1.0, 20, 5, BiasDirection::EXIT_REGION) == Approx(0.0));
    CHECK(bias(1.0, 20, 5, BiasDirection::ENTER_REGION) == Approx(0.0));

    CHECK(bias(0.1, 20, 5, BiasDirection::NO_CROSSING) == Approx(0.0));
    CHECK(bias(0.1, 20, 5, BiasDirection::EXIT_REGION) == Approx(2.1535495138));
    CHECK(bias(0.1, 20, 5, BiasDirection::ENTER_REGION) == Approx(-2.1729072697));
}

/**
 * @param outside_rejection_probability Probability that any element outside the region will be rejected
 * @param region Region to preferentially pick from
 */
RegionSampler::RegionSampler(double symmetry, std::unique_ptr<Region::RegionBase> region)
    : symmetry(symmetry)
    , region(std::move(region)) {
    if (symmetry <= pc::epsilon_dbl || symmetry > 1.0) {
        throw ConfigurationError("'symmetry' must be in the range (0,1]");
    }
}

/** Determines the direction of a transition */
BiasDirection RegionSampler::getDirection(const bool inside_before, const bool inside_after) const {
    if (inside_before && (not inside_after)) {
        return BiasDirection::EXIT_REGION;
    } else if (inside_after && (not inside_before)) {
        return BiasDirection::ENTER_REGION;
    }
    return BiasDirection::NO_CROSSING;
}

void RegionSampler::to_json(json& j) const {
    j["region"] = static_cast<json>(*region);
    j["symmetry"] = symmetry;
    if (fixed_number_inside) {
        j["Fixed <N_inside>"] = fixed_number_inside.value();
    }
}

} // namespace Faunus::SmarterMonteCarlo