#pragma once

#include <random>
#include <vector>
#include <cassert>
#include <nlohmann/json_fwd.hpp>

namespace Faunus {

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
    std::mt19937 engine; //!< Random number engine used for all operations
    Random();            //!< Constructor with deterministic seed
    void seed();         //!< Set a non-deterministic ("hardware") seed
    double operator()(); //!< Random double in uniform range [0,1)
    int range(int, int); //!< Integer in uniform range [min:max]

    /**
     * @brief Pick random element in iterable range
     * @param begin Begin iterator
     * @param end End iterator
     * @return Iterator to random element
     */
    template <typename Iterator> Iterator sample(Iterator begin, Iterator end) {
        std::advance(begin, range(0, std::distance(begin, end) - 1));
        return begin;
    }
};

void to_json(nlohmann::json &, const Random &);   //!< Random to json conversion
void from_json(const nlohmann::json &, Random &); //!< json to Random conversion

extern Random random; //!< global instance of Random

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] Random") {
    Random slump, slump2; // local instances

    CHECK(slump() == slump2()); // deterministic initialization by default; the global random variable cannot
                                // be used for comparison as its state is not reset at the beginning of
                                // each test case

    int min = 10, max = 0, N = 1e6;
    double x = 0;
    for (int i = 0; i < N; i++) {
        int j = slump.range(0, 9);
        if (j < min)
            min = j;
        if (j > max)
            max = j;
        x += j;
    }
    CHECK(min == 0);
    CHECK(max == 9);
    CHECK(std::fabs(x / N) == doctest::Approx(4.5).epsilon(0.01));

    Random r1 = R"( {"seed" : "hardware"} )"_json; // non-deterministic seed
    Random r2;                                     // default is a deterministic seed
    CHECK(r1() != r2());
    Random r3 = nlohmann::json(r1); // r1 --> json --> r3
    CHECK(r1() == r3());

    // check if random_device works
    Random a, b;
    CHECK(a() == b());
    a.seed();
    b.seed();
    CHECK(a() != b());
}
#endif

/**
 * @brief Stores a series of elements with given weight
 *
 * Elements is accessed with `get()` that will
 * randomly pick from the weighted distribution.
 * Add elements with `push_back()`
 * where the default weight is _unity_.
 */
template <typename T> class WeightedDistribution {
  private:
    std::discrete_distribution<> dist;
    std::vector<double> weights;
    size_t index; //!< index from latest element access via push_back or get
  public:
    std::vector<T> vec;                           //!< raw vector of T
    auto size() const { return vec.size(); }      //!< Number of data points
    bool empty() const { return vec.empty(); }    //!< True if no data points
    size_t getLastIndex() const { return index; } //!< index of last `get()` or `push_back()` call

    void clear() {
        vec.clear();
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
        if (auto size = std::distance(begin, end); size == vec.size()) {
            weights.resize(size);
            std::copy(begin, end, weights.begin());
            dist = std::discrete_distribution(weights.begin(), weights.end());
            assert(size_t(dist.max()) == vec.size() - 1);
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
        vec.push_back(value);
        weights.push_back(weight);
        setWeight(weights.begin(), weights.end()); // recalc. weight distribution
        index = vec.size() - 1;
    }

    /**
     * @brief Get random data point respecting the weighted distribution
     * @param engine Random number engine
     * @return Reference to data point
     */
    template <typename RandomGenerator> const T &sample(RandomGenerator &engine) {
        assert(not empty());
        index = dist(engine);
        return vec.at(index);
    }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] WeightedDistribution") {
    WeightedDistribution<double> v;

    v.push_back(0.5);
    CHECK(v.getLastIndex() == 0);
    CHECK(v.size() == 1);

    v.push_back(0.1, 4);
    CHECK(v.getLastIndex() == 1);
    CHECK(v.size() == 2);
    CHECK(not v.empty());

    int N = 1e4;
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += v.sample(random.engine);
    }
    CHECK(sum / N == doctest::Approx((0.5 * 1 + 0.1 * 4) / (1 + 4)).epsilon(0.05));

    std::vector<int> weights = {2, 1};
    v.setWeight(weights.begin(), weights.end());
    sum = 0;
    for (int i = 0; i < N; i++) {
        sum += v.sample(random.engine);
    }
    CHECK(sum / N == doctest::Approx((0.5 * 2 + 0.1 * 1) / (2 + 1)).epsilon(0.05));

    weights = {2, 1, 1};
    CHECK_THROWS(v.setWeight(weights.begin(), weights.end()));

    v.clear();
    CHECK(v.empty());
}
#endif
} // namespace Faunus
