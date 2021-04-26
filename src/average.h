#pragma once
#include <doctest/doctest.h>
#include <limits>
#include <ostream>
#include <istream>
#include <cmath>

namespace Faunus {

/**
 * @brief Class to collect averages
 *
 * @todo Should be split to a simpler class that stores only `cnt` and `sum` while
 *       rms etc should be implemented in derived classes. This could be done with
 *       the current `sqsum` method, or with automatic block averaging using an
 *       auto-correlation function.
 */
template <class T = double> class Average {
  public:
    unsigned long int cnt = 0; //!< number of values in average
    T sum = 0.0;               //!< Sum of all recorded values
    T sqsum = 0.0;             //!< Square sum of all recorded values

    void clear() { *this = Average<T>(); }                                     //!< Clear all data
    bool empty() const { return cnt == 0; }                                    //!< True if empty
    auto size() const { return cnt; }                                          //!< Number of samples
    auto avg() const { return sum / static_cast<T>(cnt); }                     //!< Average
    auto rms() const { return std::sqrt(sqsum / static_cast<T>(cnt)); }        //!< Root-mean-square
    operator double() const { return avg(); }                                  //!< Static cast operator
    bool operator<(const Average& other) const { return avg() < other.avg(); } //!< Compare operator

    auto stdev() const {
        if (cnt > 1) {
            const auto mean = avg();
            return std::sqrt((sqsum + static_cast<T>(cnt) * mean * mean - 2.0 * sum * mean) / static_cast<T>(cnt - 1));
        } else {
            return 0.0;
        }
    } //!< Standard deviation

    void add(const T value) {
        const auto value_squared = value * value;
        if (cnt == std::numeric_limits<decltype(cnt)>::max()) {
            throw std::overflow_error("max. number of samples reached");
        }
        if (std::numeric_limits<T>::max() - value_squared < sqsum) {
            throw std::overflow_error("value too large");
        }
        cnt++;
        sum += value;
        sqsum += value_squared;
    } //!< Add value

    Average& operator+=(const T value) {
        add(value);
        return *this;
    } //!< Add value

    Average& operator=(const T value) {
        clear();
        add(value);
        return *this;
    } //!< Assign to value

    /**
     * @brief Merge two averages with correct weights
     * @param other Other average
     * @return Merged average
     * @throw if numeric overflow
     */
    const Average operator+(const Average& other) const {
        if (std::numeric_limits<decltype(cnt)>::max() - other.cnt < cnt) {
            throw std::overflow_error("maximum number of samples reached");
        }
        if (std::numeric_limits<T>::max() - other.sqsum < sqsum) {
            throw std::overflow_error("value too large");
        }
        Average<T> summed_average = *this;
        summed_average.cnt += other.cnt;
        summed_average.sum += other.sum;
        summed_average.sqsum += other.sqsum;
        return summed_average;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Average<T>& average) {
        stream << average.cnt << " " << average.sum << " " << average.sqsum;
        return stream;
    } // serialize to stream

    Average<T>& operator<<(std::istream& stream) {
        stream >> cnt >> sum >> sqsum;
        return *this;
    } // de-serialize from stream

    template <class Archive> void serialize(Archive& archive) { archive(sum, sqsum, cnt); }
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
template <typename T, typename int_t = unsigned long int> class AverageObj {
#endif
  protected:
    int_t number_of_samples = 0;
    T sum; // make sure constructors zero this!
  public:
    AverageObj() : sum(T()){}; //!< Construct from empty object

    AverageObj(const T& value) : number_of_samples(1), sum(value){};

    //! Add to average
    AverageObj& operator+=(const T& value) {
        if (number_of_samples == std::numeric_limits<int_t>::max()) {
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
        } else {
            return sum;
        }
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
    int_t size() const { return number_of_samples; }
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

    AverageObjStdev(const T& value) : AverageObj<T>(value), sum_squared(value * value){};

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
        } else {
            return std::pow(sum_squared * (1.0 / static_cast<double>(number_of_samples)), 0.5);
        }
    }

    //! Standard deviation
    T stdev() const {
        if (number_of_samples == 0) {
            return T();
        } else {
            const auto N = static_cast<double>(number_of_samples);
            const auto mean = avg();
            return std::pow((sum_squared + N * mean * mean - 2.0 * sum * mean) * (1.0 / (N - 1.0)), 0.5);
        }
    }

    //!< Clear all data
    void clear() { *this = AverageObjStdev<T>(); }
};

} // namespace Faunus
