#include "smart_montecarlo.h"

namespace Faunus::SmartMonteCarlo {

/**
 * @param outside_rejection_probability Probability to reject if outside region
 * @param number_total Total number of elements
 * @param number_inside Total number of elements inside region
 * @param direction Did element exit or enter the region? Returns zero if no boundary crossing
 */
double bias(double outside_rejection_probability, const int number_total, const int number_inside,
            BiasDirection direction) {
    const auto p = outside_rejection_probability;
    switch (direction) {
    case BiasDirection::EXIT_REGION:
        return -std::log(p / (1.0 - (1.0 - p) / (p * number_total + (1.0 - p) * number_inside)));
    case BiasDirection::ENTER_REGION:
        return -std::log(1.0 / (1.0 + (1.0 - p) / (p * number_total + (1.0 - p) * number_inside)));
    case BiasDirection::NO_CROSSING:
        return 0.0;
    }
}

TEST_CASE("[Faunus] SmartMonteCarlo::bias") {
    using doctest::Approx;
    CHECK(bias(1.0, 20, 5, BiasDirection::NO_CROSSING) == Approx(0.0));
    CHECK(bias(1.0, 20, 5, BiasDirection::EXIT_REGION) == Approx(0.0));
    CHECK(bias(1.0, 20, 5, BiasDirection::ENTER_REGION) == Approx(0.0));

    CHECK(bias(0.0, 20, 5, BiasDirection::NO_CROSSING) == Approx(0.0));
    CHECK(std::isinf(bias(0.0, 20, 5, BiasDirection::EXIT_REGION)));
    CHECK(bias(0.0, 20, 5, BiasDirection::ENTER_REGION) == Approx(0.1823215568));

    CHECK(bias(0.1, 20, 5, BiasDirection::NO_CROSSING) == Approx(0.0));
    CHECK(bias(0.1, 20, 5, BiasDirection::EXIT_REGION) == Approx(2.1535495138));
    CHECK(bias(0.1, 20, 5, BiasDirection::ENTER_REGION) == Approx(0.1296778233));
}

/**
 * @param outside_rejection_probability Probability that any element outside the region will be rejected
 * @param region Region to preferentially pick from
 */
RegionSampler::RegionSampler(double outside_rejection_probability, std::unique_ptr<Region::RegionBase> region)
    : outside_rejection_probability(outside_rejection_probability), region(std::move(region)) {}

void to_json(json& j, const RegionSampler& smc) { smc.to_json(j); }
void RegionSampler::to_json(json& j) const {
    j["region"] = static_cast<json>(*region);
    j["reject_outside"] = outside_rejection_probability;
}

} // namespace Faunus::SmartMonteCarlo