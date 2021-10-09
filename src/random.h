#pragma once
#include <random>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <nlohmann/json.hpp>

#ifdef ENABLE_PCG
#include <pcg_random.hpp>
#endif

namespace Faunus {

#ifdef ENABLE_PCG
using RandomNumberEngine = pcg32;
#else
using RandomNumberEngine = std::mt19937;
#endif

/**
 * @brief Random number generator
 *
 * Example code:
 *
 * ```{.cpp}
 *     Random r1;                                     // default deterministic seed
 *     Random r2 = R"( {"seed" : "hardware"} )"_json; // non-deterministic seed
 *     r1.seed();                                     // non-deterministic seed
 * ```
 */
class Random {
  private:
    std::uniform_real_distribution<double> dist01; //!< Uniform real distribution [0,1)
  public:
    RandomNumberEngine engine; //!< Random number engine used for all operations
    Random();            //!< Constructor with deterministic seed
    void seed();         //!< Set a non-deterministic ("hardware") seed
    double operator()(); //!< Random double in uniform range [0,1)

    /**
     * @brief Integer in uniform range [min:max]
     * @param min minimum value
     * @param max maximum value
     * @return random integer in [min:max] range
     */
    template <typename IntType = int> IntType range(IntType min, IntType max) {
        static_assert(std::is_integral<IntType>::value, "Integral required");
        return std::uniform_int_distribution<IntType>(min, max)(engine);
    }

    /**
     * @brief Pick random element in iterable range
     * @param begin Begin iterator
     * @param end End iterator
     * @return Iterator to random element
     */
    template <typename Iterator> Iterator sample(Iterator begin, Iterator end) {
        std::advance(begin, range<size_t>(0, std::distance(begin, end) - 1));
        return begin;
    }
};

void to_json(nlohmann::json &, const Random &);   //!< Random to json conversion
void from_json(const nlohmann::json &, Random &); //!< json to Random conversion

extern Random random; //!< global instance of Random

/**
 * @brief Stores a series of elements with given weight
 *
 * Elements is accessed with `get()` that will
 * randomly pick from the weighted distribution.
 * Add elements with `addGroup()`
 * where the default weight is _unity_.
 *
 * @tparam T Data type to store
 */
template <typename T> class WeightedDistribution {
  private:
    std::discrete_distribution<> distribution;
    std::vector<double> weights; //!< weights for each data point
    size_t latest_index;         //!< index from latest element access via addGroup or get
  public:
    std::vector<T> data;                        //!< raw vector of T
    auto size() const { return data.size(); }   //!< Number of data points
    bool empty() const { return data.empty(); } //!< True if no data points
    size_t getLastIndex() const {
        assert(!data.empty());
        return latest_index;
    } //!< index of last `get()` or `addGroup()` element

    void clear() {
        data.clear();
        weights.clear();
    } //!< Clear all data

    /**
     * @brief Set weights for all data points in `vec`
     * @param begin Begin iterator for weights
     * @param end End iterator for weights
     *
     * The iterable range must match the `size()` of the stored data,
     * otherwise an exception is thrown.
     */
    template <typename Iterator> void setWeight(Iterator begin, Iterator end) {
        if (auto size = std::distance(begin, end); size == data.size()) {
            weights.resize(size);
            std::copy(begin, end, weights.begin());
            distribution = std::discrete_distribution(weights.begin(), weights.end());
            assert(size_t(distribution.max()) == data.size() - 1);
        } else {
            throw std::runtime_error("number of weights must match data");
        }
    }

    /**
     * @brief Add data and it's weight
     * @param value data
     * @param weight weight (default: 1.0)
     */
    void push_back(const T &value, double weight = 1.0) {
        data.push_back(value);
        weights.push_back(weight);
        setWeight(weights.begin(), weights.end());
        latest_index = data.size() - 1;
    }

    /**
     * @brief Get random data point respecting the weighted distribution
     * @param engine Random number engine
     * @return Reference to data point
     */
    template <typename RandomGenerator> const T &sample(RandomGenerator &engine) {
        assert(not empty());
        latest_index = distribution(engine);
        return data.at(latest_index);
    }
};
} // namespace