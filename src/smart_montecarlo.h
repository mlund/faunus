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
    int n_total = 0;
    int n_inside = 0;
    bool is_inside = false;
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
    : item(&item), n_total(number_total), n_inside(number_inside), is_inside(item_is_inside) {}

enum class BiasDirection { ENTER_REGION, EXIT_REGION, NO_CROSSING }; //!< Controls the bias direction

double bias(double outside_rejection_probability, int n_total, int n_inside,
            BiasDirection direction); //!< Bias energy (kT)

/**
 * Helper class for smart MC moves that preferentially samples within Regions
 *
 * Responsibilities:
 *
 * - Randomly select groups or particles with respect to an arbitrary `Region`
 * - Calculate the corresponding energy bias due to the non-uniform sampling
 */
class RegionSampler {
  private:
    double symmetry = 1.0; //!< Or "p" between ]0:1]; 1 --> uniform sampling (no regional preference)
  protected:
    std::unique_ptr<Region::RegionBase> region; //!< This defines the smart MC region
    std::optional<int> fixed_number_inside;     //!< Optionally give a fixed number inside
    template <typename Range> int getNumberInside(Range& range) const;

  public:
    RegionSampler(double symmetry, std::unique_ptr<Region::RegionBase> region);
    virtual ~RegionSampler() = default;
    void to_json(json& j) const; //!< Serialise to json
    template <typename T, typename Range> std::optional<Selection<T>> select(Range& range, Random& random);
    template <typename T> double bias(const Selection<T>& selection);
};

/**
 * Determines the number of items inside the region. If `fixed_number_inside` is set,
 * this will be used and thus skip the potentially expensive conditional count.
 */
template <typename Range> int RegionSampler::getNumberInside(Range& range) const {
    if (fixed_number_inside) {
        return fixed_number_inside.value();
    }
    return ranges::count_if(range, [&](const auto& i) { return region->inside(i); }); // expensive
}

/**
 * @brief Select a random element from range and perform region check. Typically *before* a move.
 * @tparam Range range of particles or groups
 * @tparam T typically `Particle` or `Group`
 * @param range Range of particles or groups
 * @param random Random number object used to perform a random selection
 * @return Optional selection object
 *
 * This will loop until a perturbable element is found. If a max number of attempts is reached
 * without any selection, an empty object will be returned and a warning issued. The maximum
 * number of attempts is currently set to ten times the size of the given range.
 */
template <typename T, typename Range> std::optional<Selection<T>> RegionSampler::select(Range& range, Random& random) {
    const auto n_total = ranges::distance(range.begin(), range.end());
    int selection_attemps = 10 * n_total;
    do {
        auto it = random.sample(range.begin(), range.end()); // random particle or group
        if (it == range.end()) {
            return std::nullopt;
        }
        const auto inside = region->inside(*it); // is element inside or outside region?
        if (not inside && symmetry < random()) {
            continue;
        }
        return Selection<T>(*it, n_total, getNumberInside(range), inside);
    } while (selection_attemps-- > 0);
    faunus_logger->warn("Max number selection attempts reached. Increase 'symmetry'?");
    return std::nullopt;
}

/**
 * @tparam T Typically `Particle` or `Group`
 * @returns Bias energy (kT) always zero if no `selection` has been set
 *
 * This function is typically called *after* a MC move and will determine the
 * bias due to non-uniform sampling.
 */
template <typename T> double RegionSampler::bias(const Selection<T>& selection) {
    BiasDirection direction;
    const auto moved_item_is_inside = region->inside(*(selection.item)); // may have changed due to move
    if (selection.is_inside && (!moved_item_is_inside)) {
        direction = BiasDirection::EXIT_REGION;
    } else if (moved_item_is_inside && (!selection.is_inside)) {
        direction = BiasDirection::ENTER_REGION;
    } else {
        direction = BiasDirection::NO_CROSSING;
    }
    return SmartMonteCarlo::bias(symmetry, selection.n_total, selection.n_inside, direction);
}

/**
 * Helper class for constructing smart monte carlo moves
 */
template <typename T> class MoveSupport {
  private:
    using OptionalElement = std::optional<std::reference_wrapper<T>>; //!< Reference to selected element
    std::optional<SmartMonteCarlo::Selection<T>> selection; //!< Contains data on currently selected group (if any)
    SmartMonteCarlo::RegionSampler region_sampler;
    int bias_update_interval = 100;
    Average<double> mean_bias_energy;
    Average<double> mean_selected_inside;    //!< Fraction of selections found inside
    AverageStdev<double> mean_number_inside; //!< Average number of groups found inside region
    void analyseSelection(int count_inside); //!< Track and analyze inside count

  public:
    MoveSupport(const Space& spc, const json& j);
    double bias();
    void to_json(json& j) const;
    template <typename Range> OptionalElement select(Range& mollist, Random& random);
};

template <typename T>
MoveSupport<T>::MoveSupport(const Space& spc, const json& j)
    : region_sampler(j.at("symmetry").get<double>(), Region::createRegion(spc, j)) {}

/**
 * This is used to sample the average number of groups inside the region. Once converged, the
 * expensive region search can be disabled by giving `fixed_count_inside` a mean value which will
 * be used for all further bias evaluations.
 */
template <typename T> void MoveSupport<T>::analyseSelection(int count_inside) {
    assert(selection);
    mean_number_inside += static_cast<double>(count_inside);
    mean_selected_inside += static_cast<double>(selection->is_inside);

    if (mean_number_inside.size() % bias_update_interval == 0) {
        if (mean_number_inside.stdev() / mean_number_inside.avg() < 0.05) {
            // fixed_count_inside = static_cast<int>(mean_count_inside.avg());
            // faunus_logger->info("Stopping bias update since threshold reached.");
        }
    }
}
template <typename T> double MoveSupport<T>::bias() {
    auto bias_energy(0.0);
    if (selection) {
        analyseSelection(selection->n_inside);
        bias_energy = region_sampler.bias(*selection);
    }
    mean_bias_energy += bias_energy;
    return bias_energy;
}
template <typename T> void MoveSupport<T>::to_json(json& j) const {
    region_sampler.to_json(j);
    if (!mean_number_inside.empty()) {
        j["mean number inside"] = mean_number_inside.avg();
    }
    if (!mean_bias_energy.empty()) {
        j["mean bias energy (kT)"] = mean_bias_energy.avg();
    }
    if (!mean_selected_inside.empty()) {
        j["selected inside fraction"] = mean_selected_inside.avg();
    }
}

template <typename T>
template <typename Range>
typename MoveSupport<T>::OptionalElement MoveSupport<T>::select(Range& mollist, Random& random) {
    selection = region_sampler.select<T>(mollist, random);
    if (selection) {
        return *(selection->item);
    }
    return std::nullopt;
}

} // namespace Faunus::SmartMonteCarlo
