#pragma once

#include "regions.h"
#include <optional>
#include <range/v3/algorithm/count_if.hpp>

namespace Faunus::SmarterMonteCarlo {

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
    : item(&item)
    , n_total(number_total)
    , n_inside(number_inside)
    , is_inside(item_is_inside) {}

enum class BiasDirection { ENTER_REGION, EXIT_REGION, NO_CROSSING }; //!< Controls the bias direction

double bias(double outside_acceptance, int n_total, int n_inside,
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
    const double outside_acceptance = 1.0; //!< Or "p" between ]0:1]; 1 --> uniform sampling (no regional preference)
    BiasDirection getDirection(bool inside_before, bool inside_after) const;
    template <typename Range> double getNumberInside(Range& range) const;

  protected:
    const std::unique_ptr<Region::RegionBase> region; //!< This defines the smart MC region

  public:
    RegionSampler(double outside_acceptance, std::unique_ptr<Region::RegionBase> region);
    virtual ~RegionSampler() = default;
    void to_json(json& j) const; //!< Serialise to json
    template <typename T, typename Range> std::optional<Selection<T>> select(Range& range, Random& random);
    template <typename T> double bias(const Selection<T>& selection);
    std::optional<double> fixed_number_inside; //!< Optionally give a fixed number inside
};

/**
 * Determines the number of items inside the region. If `fixed_number_inside` is set,
 * this will be used and thus skip the potentially expensive conditional count.
 *
 * @returns (mean) number of elements inside region
 */
template <typename Range> double RegionSampler::getNumberInside(Range& range) const {
    if (fixed_number_inside) {
        return fixed_number_inside.value();
    }
    return static_cast<double>(ranges::count_if(range, [&](const auto& i) { return region->inside(i); })); // expensive
}

/**
 * @brief Select a random element from range and perform region check. Typically *before* a move.
 * @tparam Range range of particles or groups
 * @tparam T typically `Particle` or `Group`
 * @param range Range of particles or groups
 * @param random Random number object used to perform a random selection
 * @return Optional selection object
 *
 * Loops until a perturbable element is found. If a max number of attempts is reached
 * without any selection, an empty object is returned and a warning issued. The maximum
 * number of attempts is currently set to 10x the size of the given range.
 */
template <typename T, typename Range> std::optional<Selection<T>> RegionSampler::select(Range& range, Random& random) {
    const auto n_total = ranges::distance(range.begin(), range.end());
    int max_selection_attempts = 10 * n_total;
    do {
        auto it = random.sample(range.begin(), range.end()); // random particle or group
        if (it == range.end()) {
            return std::nullopt;
        }
        const auto inside = region->inside(*it); // is element inside or outside region?
        if (not inside && outside_acceptance < random()) {
            continue;
        }
        return Selection<T>(*it, n_total, getNumberInside(range), inside);
    } while (max_selection_attempts-- > 0);
    faunus_logger->warn("Max selection attempts reached. Increase outside acceptance (p) ?");
    return std::nullopt;
}

/**
 * @tparam T Typically `Particle` or `Group`
 * @returns Bias energy (kT) always zero if no `selection` has been set
 *
 * This function is typically called *after* a MC move and will determine the
 * bias due to non-uniform sampling which depends on the transition directions,
 * i.e. if an element is moved from _inside_ of the region to the _outside_ etc.
 */
template <typename T> double RegionSampler::bias(const Selection<T>& selection) {
    const auto is_inside_after = region->inside(*(selection.item)); // may have changed due to move
    const auto direction = getDirection(selection.is_inside, is_inside_after);
    return SmarterMonteCarlo::bias(outside_acceptance, selection.n_total, selection.n_inside, direction);
}

/**
 * Helper class for constructing smart monte carlo moves
 */
template <typename T> class MoveSupport {
  private:
    using OptionalElement = std::optional<std::reference_wrapper<T>>; //!< Reference to selected element
    std::optional<SmarterMonteCarlo::Selection<T>> selection; //!< Contains data on currently selected group (if any)
    SmarterMonteCarlo::RegionSampler region_sampler;          //!< Biased sampling of elements
    Average<double> mean_bias_energy;                         //!< Statistics of bias energy
    Average<double> mean_selected_inside;                     //!< Fraction of selections found inside
    AverageStdev<double> mean_number_inside;                  //!< Average number of groups found inside region
    AverageStdev<double> blocks_inside;                       //!< Used for block averaging of `mean_number_inside`
    void analyseSelection(int count_inside);                  //!< Track and analyze inside count

  public:
    MoveSupport(const Space& spc, const json& j);
    double bias();
    void to_json(json& j) const;
    template <typename Range> OptionalElement select(Range& mollist, Random& random);
};

template <typename T>
MoveSupport<T>::MoveSupport(const Space& spc, const json& j)
    : region_sampler(j.at("p").get<double>(), Region::createRegion(spc, j)) {}

/**
 * This is used to sample and analyse the average number of groups inside the region.
 * Once converged, the expensive region search can be disabled by giving
 * `fixed_count_inside` a mean value which will be used for further bias evaluations.
 *
 * @todo convergence thresholds currently hardcoded; bootstrapping would be better...
 */
template <typename T> void MoveSupport<T>::analyseSelection(int count_inside) {
    assert(selection);
    mean_selected_inside += static_cast<double>(selection->is_inside);

    if (region_sampler.fixed_number_inside) {
        return; // skip the rest as we have already fixed n_inside
    }

    mean_number_inside += static_cast<double>(count_inside);
    if (mean_number_inside.size() % 100 == 0) {                            // build up a block average
        blocks_inside += mean_number_inside.avg();                         // add a new block
        if (blocks_inside.size() % 2 == 0 && blocks_inside.rsd() < 0.04) { // check if block has converged
            region_sampler.fixed_number_inside = blocks_inside.avg();      // disable all further sampling!
            faunus_logger->info("Regional <N_inside> = {:.1f} converged and settled after {} iterations",
                                region_sampler.fixed_number_inside.value(), mean_number_inside.size());
        }
    }
}

/**
 * Returns the bias energy. Requires the total number
 * of elements inside the region, Nin, which is explicitly sampled
 * until Nin is converged whereafter a constant, mean value is assigned.
 * (speed optimization; approximation)
 *
 * @return Bias energy (kT)
 */
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
    if (mean_number_inside) {
        j["mean number inside"] = mean_number_inside.avg();
    }
    if (mean_bias_energy) {
        j["mean bias energy (kT)"] = mean_bias_energy.avg();
    }
    if (mean_selected_inside) {
        j["selected inside fraction"] = mean_selected_inside.avg();
    }
}

/**
 * @tparam T Element type, typically Group or Particle
 * @param mollist Range of elements
 * @param random Random number generator
 * @return Optional reference to selected element
 */
template <typename T>
template <typename Range>
typename MoveSupport<T>::OptionalElement MoveSupport<T>::select(Range& mollist, Random& random) {
    selection = region_sampler.select<T>(mollist, random);
    if (selection) {
        return *(selection->item);
    }
    return std::nullopt;
}

} // namespace Faunus::SmarterMonteCarlo
