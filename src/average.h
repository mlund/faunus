#pragma once
#include <doctest/doctest.h>
#include <limits>
#include <ostream>
#include <istream>
#include <cmath>
#include <fstream>
#include <spdlog/fmt/fmt.h>
#include <nlohmann/json.hpp>

namespace Faunus {

/**
 * @brief Class to collect averages
 * @todo replace static assert w. concept in c++20
 */
template <class value_type = double, class counter_type = unsigned long int> class Average {
    static_assert(std::is_floating_point<value_type>::value, "floating point type required");
    static_assert(std::is_unsigned<counter_type>::value, "unsigned integer required");

  protected:
    counter_type number_of_samples = 0; //!< number of values in average
    value_type value_sum = 0.0;         //!< Sum of all recorded values
  public:
    void clear() { *this = Average(); }                                                 //!< Clear all data
    bool empty() const { return number_of_samples == 0; }                               //!< True if empty
    auto size() const { return number_of_samples; }                                     //!< Number of samples
    auto avg() const { return value_sum / static_cast<value_type>(number_of_samples); } //!< Average
    explicit operator value_type() const { return avg(); }                              //!< Static cast operator
    bool operator<(const Average& other) const { return avg() < other.avg(); }          //!< Compare means

    /**
     * @brief Add value to average
     * @param value Value to add
     * @throw If overflow in either the counter, or sum of squared values
     */
    void add(const value_type value) {
        if (number_of_samples == std::numeric_limits<counter_type>::max()) {
            throw std::overflow_error("max. number of samples reached");
        }
        value_sum += value;
        number_of_samples++;
    }

    /**
     * @brief Add value to average
     * @param value Value to add
     * @throw If overflow in either the counter, or sum of squared values
     */
    Average& operator+=(const value_type value) {
        add(value);
        return *this;
    }

    bool operator==(const Average& other) const { return (size() == other.size()) && (value_sum == other.value_sum); }

    /** Clear and assign a new value */
    Average& operator=(const value_type value) {
        clear();
        add(value);
        return *this;
    }

    /**
     * @brief Merge two averages with correct weights
     * @param other Other average
     * @return Merged average
     * @throw if numeric overflow
     */
    auto operator+(const Average& other) const {
        if (std::numeric_limits<counter_type>::max() - other.number_of_samples < number_of_samples) {
            throw std::overflow_error("maximum number of samples reached");
        }
        if (std::numeric_limits<value_type>::max() - other.value_sum < value_sum) {
            throw std::overflow_error("value too large");
        }
        Average summed_average = *this;
        summed_average.number_of_samples += other.number_of_samples;
        summed_average.value_sum += other.value_sum;
        return summed_average;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Average& average) {
        stream << average.number_of_samples << " " << average.value_sum;
        return stream;
    } // serialize to stream

    Average& operator<<(std::istream& stream) {
        stream >> number_of_samples >> value_sum;
        return *this;
    } // de-serialize from stream

    template <class Archive> void serialize(Archive& archive) { archive(value_sum, number_of_samples); }
};

/**
 * @brief Class to collect averages and standard deviation
 * @todo inherit from `Average`
 */
template <class value_type = double, class counter_type = unsigned long int>
class AverageStdev : public Average<value_type, counter_type> {
  private:
    value_type squared_value_sum = 0.0; //!< Square sum of all recorded values
    using base = Average<value_type, counter_type>;
    using base::number_of_samples;
    using base::value_sum;

  public:
    using base::avg;
    using base::empty;
    void clear() { *this = AverageStdev(); } //!< Clear all data
    auto rms() const {
        return std::sqrt(squared_value_sum / static_cast<value_type>(number_of_samples));
    } //!< Root-mean-square

    value_type stdev() const {
        if (empty()) {
            return 0.0;
        }
        const auto mean = avg();
        return std::sqrt(
            (squared_value_sum + static_cast<value_type>(number_of_samples) * mean * mean - 2.0 * value_sum * mean) /
            static_cast<value_type>(number_of_samples - 1));
    } //!< Standard deviation

    /**
     * @brief Add value to average
     * @param value Value to add
     * @throw If counter overflow
     */
    void add(const value_type value) {
        base::add(value);
        squared_value_sum += value * value;
    }

    /**
     * @brief Add value to average
     * @param value Value to add
     * @throw If overflow in either the counter, or sum of squared values
     */
    AverageStdev& operator+=(const value_type value) {
        add(value);
        return *this;
    }

    bool operator==(const AverageStdev& other) const {
        return (this->size() == other.size()) && (value_sum == other.value_sum) &&
               (squared_value_sum == other.squared_value_sum);
    }

    /** Clear and assign a new value */
    AverageStdev& operator=(const value_type value) {
        clear();
        add(value);
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& stream, const AverageStdev& average) {
        stream << average.number_of_samples << " " << average.value_sum << " " << average.squared_value_sum;
        return stream;
    } // serialize to stream

    AverageStdev& operator<<(std::istream& stream) {
        stream >> number_of_samples >> value_sum >> squared_value_sum;
        return *this;
    } // de-serialize from stream

    template <class Archive> void serialize(AverageStdev& archive) {
        archive(value_sum, squared_value_sum, number_of_samples);
    }
};

/**
 * @brief "Decorrelation" class from https://dx.doi.org/10.1002/jcc.20746
 *
 * See also Kent's dissertation, page 66:
 * https://thesis.library.caltech.edu/748/1/david_randall_kent_iv-dissertation.pdf
 *
 * Pecularities of the pseudo-code in the article:
 *
 * 1. vector index begins at 1 (see e.g. comparisons with `size`)
 * 2. is it really needed to initialize vectors with one element?
 * 3. The indentation of the following code should likely be moved one level to the left:
 * ~~~
 *    carry = new_sample
 *    i = 1
 *    done = false
 * ~~~
 *
 * @todo Bundle current vectors into single vector of struct; unittests; check correctness
 * @warning Under construction and untested
 */
template <class value_type = double, class counter_type = unsigned long int> class Decorrelation {
  private:
    using Statistic = AverageStdev<value_type, counter_type>;
    unsigned int nsamples = 0;
    std::vector<Statistic> blocked_statistics = {Statistic()};
    std::vector<value_type> waiting_sample = {0.0};
    std::vector<bool> waiting_sample_exists = {false};

  public:
    Decorrelation& add(value_type new_sample) {
        nsamples++;
        if (nsamples >= static_cast<unsigned int>(std::pow(2, blocked_statistics.size()))) {
            blocked_statistics.template emplace_back();
            waiting_sample.push_back(0.0);
            waiting_sample_exists.push_back(false);
        }
        auto carry = new_sample;
        blocked_statistics.at(0) += new_sample;

        for (size_t i = 1; i < blocked_statistics.size(); i++) {
            if (waiting_sample_exists.at(i)) {
                new_sample = 0.5 * (waiting_sample.at(i) + carry);
                carry = new_sample;
                blocked_statistics.at(i) += new_sample;
                waiting_sample_exists.at(i) = false;
            } else {
                waiting_sample_exists.at(i) = true;
                waiting_sample.at(i) = carry;
                break;
            }
        }
        return *this;
    }

    auto size() const { return nsamples; }
    bool empty() const { return nsamples == 0; }

    void to_disk(const std::string& filename) const {
        if (auto stream = std::ofstream(filename); stream) {
            for (size_t i = 0; i < blocked_statistics.size(); i++) {
                stream << fmt::format("{} {} {}\n", i, blocked_statistics[i].avg(), blocked_statistics[i].stdev());
            }
        }
    }

    void to_json(nlohmann::json& j) const {
        j = nlohmann::json::array();
        for (size_t i = 0; i < blocked_statistics.size(); i++) {
            j.push_back(nlohmann::json::array({blocked_statistics[i].avg(), blocked_statistics[i].stdev()}));
        }
    }
};

#if __cplusplus > 201703L
/** Requirements for types used with `AverageObj` */
template <class T> concept Averageable = requires(T a) {
    a * 1.0; // please implement `T operator*(double) const`
    a* a;    // please implement `T operator*(const T&) const`
    a += a;  // please implement `T& operator+=(const T&)`
};
#endif
    /**
     * @brief Simple class to average data contained in objects
     * @tparam T Type to average
     * @tparam int_t Unsigned interger type
     *
     * It is required that `T` has the following operator overloads:
     * - `T operator*(double) const`
     * - `T operator*(const T&) const`
     * - `T& operator+=(const &T)`
     */
#if __cplusplus > 201703L
template <Averageable T, typename int_t = unsigned long int> class AverageObj {
#else
template <typename T, typename counter_type = unsigned long int> class AverageObj {
#endif
  protected:
    counter_type number_of_samples = 0;
    T sum; // make sure constructors zero this!
  public:
    AverageObj() : sum(T()){}; //!< Construct from empty object

    AverageObj(const T& value) : number_of_samples(1), sum(value){};

    //! Add to average
    AverageObj& operator+=(const T& value) {
        if (number_of_samples == std::numeric_limits<counter_type>::max()) {
            throw std::overflow_error("maximum samples reached");
        }
        sum += value;
        ++number_of_samples;
        return *this;
    }
    //! Calculate average
    T avg() const {
        if (number_of_samples > 0) {
            return sum * (1.0 / static_cast<double>(number_of_samples));
        }
        return sum;
    }

    //! Convert to T
    operator T() const { return avg(); }

    //! Compare operator
    bool operator<(const AverageObj& other) const { return avg() < other.avg(); }

    //! True if empty
    bool empty() const { return number_of_samples == 0; }

    //!< Clear all data
    void clear() { *this = AverageObj<T>(); }

    //! Number of samples
    auto size() const { return number_of_samples; }
};

/*
 * Experimental extension of AverageObj that includes stdev() and rms(). This
 * requires a pow(T, double) function for squaring and taking the square-root.
 * How could this best be implemented and maintain compatibility when T=double?
 */
template <typename T, typename int_t = unsigned long int> class AverageObjStdev : public AverageObj<T> {
  private:
    T sum_squared; // make sure constructors zero this!
    using AverageObj<T>::avg;
    using AverageObj<T>::sum;
    using AverageObj<T>::number_of_samples;

  public:
    AverageObjStdev() : sum_squared(T()){}; //!< Construct from empty object

    explicit AverageObjStdev(const T& value) : AverageObj<T>(value), sum_squared(value * value){};

    //! Add to average
    AverageObjStdev& operator+=(const T& value) {
        AverageObj<T>::operator+=()(value);
        sum_squared += std::pow(value, 2);
        return *this;
    }

    //! Root-mean-square
    T rms() const {
        if (number_of_samples == 0) {
            return T();
        }
        return std::pow(sum_squared * (1.0 / static_cast<double>(number_of_samples)), 0.5);
    }

    //! Standard deviation
    T stdev() const {
        if (number_of_samples == 0) {
            return T();
        }
        const auto N = static_cast<double>(number_of_samples);
        const auto mean = avg();
        return std::pow((sum_squared + N * mean * mean - 2.0 * sum * mean) * (1.0 / (N - 1.0)), 0.5);
    }

    //!< Clear all data
    void clear() { *this = AverageObjStdev<T>(); }
};

} // namespace Faunus
