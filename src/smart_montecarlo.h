#pragma once

#include "group.h"
#include "regions.h"
#include <optional>
#include <range/v3/algorithm/count_if.hpp>

namespace Faunus::SmartMonteCarlo {

/**
 * @brief Particle or groups selection for smart Monte Carlo (SMC) moves
 *
 * This contains data useful for performing a SMC move. This object should be created
 * only with `SmartMonteCarlo::select` and hence has a protected constructor.
 */
template <typename T> class Selection {
  public:
    T* item = nullptr; //!< Selected group, particle etc.
    int number_total = 0;
    int number_inside = 0;
    bool item_is_inside = false;
    Selection(T& item, int number_total, int number_inside, bool item_is_inside);
};
/**
 * @tparam T Particle or group (`item')
 * @param item Reference to item
 * @param number_total Total number of items
 * @param number_inside Number of items inside region
 * @param item_is_inside Is current `item' inside region?
 */
template <typename T>
Selection<T>::Selection(T& item, int number_total, int number_inside, bool item_is_inside)
    : item(&item), number_total(number_total), number_inside(number_inside), item_is_inside(item_is_inside) {}

enum class BiasDirection { ENTER_REGION, EXIT_REGION, NO_CROSSING }; //!< Controls the bias direction

double bias(double outside_rejection_probability, int number_total, int number_inside,
            BiasDirection direction); //!< Bias energy (kT)

/**
 * Helper class for smart MC moves that preferentially samples within Regions
 *
 * Responsibilities:
 * - Randomly select groups or particles with respect to an arbitrary `Region`
 * - Calculate the corresponding energy bias due to the non-uniform sampling
 */
class RegionSampler {
  public:
    friend void to_json(json&, const RegionSampler&);

  private:
    double outside_rejection_probability = 1.0; //!< Probability to reject a particle outside region
    std::unique_ptr<Region::RegionBase> region; //!< This defines the smart MC region
  protected:
    virtual void to_json(json& j) const; //!< Serialise to json

  public:
    std::optional<int> fixed_count_inside; //!< Set this to speed up (skip) region check

    RegionSampler(double outside_rejection_probability, std::unique_ptr<Region::RegionBase> region);
    virtual ~RegionSampler() = default;

    /**
     * @brief Select a random element from range and perform region check. Typically *before* a move.
     * @tparam Range range of particles or groups
     * @tparam T typically `Particle` or `Group`
     * @param range Range of particles or groups
     * @param random Random number object used to perform a random selection
     * @return Optional selection object. May be empty if outside and `outside_rejection_probability < 1`
     */
    template <typename T, typename Range> std::optional<Selection<T>> select(Range& range, Random& random) {
        auto item = random.sample(range.begin(), range.end()); // random particle or group
        if (item == range.end()) {
            return std::nullopt;
        }
        const auto item_is_inside = region->inside(*item);
        if (!item_is_inside && random() < outside_rejection_probability) {
            return std::nullopt; // reject outside items w. `outside_rejection_probability`
        }

        const auto number_total = ranges::distance(range.begin(), range.end());
        const auto number_inside = getNumberInside(range); // number of elements inside region

        return Selection<T>(*item, number_total, number_inside, item_is_inside);
    }

    /**
     * Determines the number of items inside the region. If `fixed_count_inside` is set,
     * this will be used and thus skip the potentially expensive conditional count.
     */
    template <typename Range> int getNumberInside(Range& range) const {
        if (fixed_count_inside) { // is it already fixed?
            return fixed_count_inside.value();
        }
        return ranges::count_if(range, [&](const auto& i) { return region->inside(i); }); // expensive
    }

    /**
     * @tparam T Typically `Particle` or `Group`
     * @returns Bias energy (kT) always zero if no `selection` has been set
     *
     * This function is typically called *after* a MC move and will determine the
     * bias due to non-uniform sampling.
     */
    template <typename T> double bias(const Selection<T>& selection) {
        BiasDirection direction;
        const auto moved_item_is_inside = region->inside(*(selection.item)); // may have changed due to move
        if (selection.item_is_inside && !moved_item_is_inside) {
            direction = BiasDirection::EXIT_REGION;
        } else if (!selection.item_is_inside && moved_item_is_inside) {
            direction = BiasDirection::ENTER_REGION;
        } else {
            direction = BiasDirection::NO_CROSSING;
        }
        return SmartMonteCarlo::bias(outside_rejection_probability, selection.number_total, selection.number_inside,
                                     direction);
    }
};

void to_json(json& j, const RegionSampler& smc);

} // namespace Faunus::SmartMonteCarlo
