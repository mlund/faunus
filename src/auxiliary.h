#pragma once

#include "average.h"
#include <functional>
#include <fstream>
#include <vector>
#include <string>

/**
 * @file auxiliary.h
 *
 * This file contains auxiliary functionality that
 * have no dependencies other than STL and can hence
 * be copied to other projects.
 */
namespace Faunus {

/**
 * @brief Convert floating point number to integral number. Perform range check and rounding.
 * @tparam TOut integral type
 * @tparam TIn floating point type
 * @param number (floating point type)
 * @return number (integral type)
 * @throw std::overflow_error
 */
template <typename TOut, typename TIn> inline TOut numeric_cast(const TIn number) {
    static_assert(std::is_floating_point<TIn>::value, "TIn must be floating point.");
    static_assert(std::is_integral<TOut>::value, "TOut must be integer.");
    if (std::isfinite(number)) {
        // The number is finite ...
        if (number < std::nextafter(static_cast<TIn>(std::numeric_limits<TOut>::max()), 0) &&
            number > std::nextafter(static_cast<TIn>(std::numeric_limits<TOut>::min()), 0)) {
            // ... and fits into the integral type range.
            // The nextafter function is used to mitigate possible rounding up or down.
            return static_cast<TOut>(number >= 0 ? number + 0.5 : number - 0.5); // round before cast
        }
    }
    // not-a-number, infinite, or outside the representable range
    throw std::overflow_error("numeric cast overflow");
}

/**
 * @brief Iterate over pairs in container, call a function on the elements, and sum results
 * @tparam T Floating point type. Default: `double)`
 * @param begin Begin iterator
 * @param end End iterator
 * @param f Function to apply to the pair
 * @param aggregator Function to aggregate the result from each pair. Default: `std::plus<T>`
 */
template <typename Titer, typename Tfunction, typename T = double, typename Taggregate_function>
T for_each_unique_pair(Titer begin, Titer end, Tfunction f, Taggregate_function aggregator = std::plus<T>()) {
    T x = T();
    for (auto i = begin; i != end; ++i) {
        for (auto j = i; ++j != end;) {
            x = aggregator(x, f(*i, *j));
        }
    }
    return x;
}
#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] for_each_pair") {
    int x;
    std::vector<int> a = {1, 2, 3};
    x = for_each_unique_pair(
        a.begin(), a.end(), [](int i, int j) { return i * j; }, std::plus<int>());
    CHECK(x == 2 + 3 + 6);
    a.resize(1);
    x = for_each_unique_pair(
        a.begin(), a.end(), [](int i, int j) { return i * j; }, std::plus<int>());
    CHECK(x == 0);
}
#endif

/** @brief Erase from container `a` all values found in container `b` */
template <typename T> T erase_range(T a, const T &b) {
    a.erase(std::remove_if(a.begin(), a.end(),
                           [b](typename T::value_type i) { return std::find(b.begin(), b.end(), i) != b.end(); }),
            a.end());
    return a;
}

/**
 * @brief Ordered pair where `first<=second`
 *
 * Upon construction, the smallest element is placed in `first`
 * so that `ordered_pair<int>(i,j)==ordered_pair<int>(j,i)` is always true.
 *
 * @todo Add std::pair copy operator
 */
template <class T> struct ordered_pair : public std::pair<T, T> {
    using base = std::pair<T, T>;
    ordered_pair() = default;
    ordered_pair(const T &a, const T &b) : base(std::minmax(a, b)) {}
    bool contains(const T &value) const {
        if (value != base::first && value != base::second) {
            return false;
        } else {
            return true;
        }
    }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] ordered_pair") {
    ordered_pair<int> a = {1, 2}, b = {2, 1};
    CHECK((a.first == 1 && a.second == 2));
    CHECK((b.first == 1 && b.second == 2));
    CHECK(a == b);
    CHECK(a.contains(1));
    CHECK(a.contains(2));
    CHECK(!a.contains(3));
}
#endif

/**
 * @brief Round a floating point to integer for use with histogram binning.
 *
 * This will round to nearest integer, assuming a certain resolution
 * (default 1). This can be useful for binning floating point data
 * into a histogram or table of resolution `dx` (second argument).
 *
 * Example:
 *
 * ~~~~
 *     double x=21.3;
 *     to_bin(x);     // -> 21
 *     to_bin(x,2);   // -> 11
 *     to_bin(x,0.5); // -> 43
 * ~~~~
 *
 */
template <class T, class Tint = int> Tint to_bin(T x, T dx = 1) {
    return (x < 0) ? Tint(x / dx - 0.5) : Tint(x / dx + 0.5);
}

/**
 * @brief Round and bin numbers for use with tables and histograms
 *
 * This will round a float to nearest value divisible by `dx` or
 * convert to an integer corresponding to a binning value. In the
 * latter case, a minimum value must be specified upon construction.
 *
 * Example:
 *
 * ~~~{.cpp}
 * Quantize<double> Q(0.5, -0.5);
 * Q = 0.61;
 * double x = Q;      // --> 0.5
 * int bin = Q;       // --> 2
 * bin = Q(-0.5);     // --> 0
 * x = Q.frombin(2);  // --> 0.5
 *
 * std::vector<bool> v(10);
 * v[ Q(0.61) ] = false; // 2nd element set to false
 * ~~~
 *
 */
template <class Tfloat = double> class Quantize {
  private:
    Tfloat xmin, dx, x;

  public:
    /**
     * @brief Constructor
     * @param dx resolution
     * @param xmin minimum value if converting to integral type (binning)
     */
    Quantize(Tfloat dx, Tfloat xmin = 0) : xmin(xmin), dx(dx) {}

    /** @brief Assigment operator */
    Quantize &operator=(Tfloat val) {
        assert(val >= xmin);
        if (val >= 0) {
            x = int(val / dx + 0.5) * dx;
        } else {
            x = int(val / dx - 0.5) * dx;
        }
        assert(x >= xmin);
        return *this;
    }

    Quantize &frombin(unsigned int i) {
        x = i * dx + xmin;
        return *this;
    }

    /** @brief Assignment with function operator */
    Quantize &operator()(Tfloat val) {
        *this = val;
        return *this;
    }

    /** @brief Implicit convertion to integral (bin) or float (rounded) */
    template <typename T> operator T() {
        if (std::is_integral<T>::value) {
            return T((x - xmin) / dx + 0.5);
        }
        return x;
    }
};

/**
 * @brief Convert whitespace separated words into vector of given type
 *
 * Example:
 *
 * ~~~~
 * auto v = textio::words2vec<double>( "0.2 1 100" );
 * for (auto i : v)
 *   cout << 2*i << " "; // -> 0.4 2 200
 * ~~~~
 */
template <class T> std::vector<T> words2vec(const std::string &string_of_words) {
    auto number_of_words =
        std::distance(std::istream_iterator<std::string>(std::istringstream(string_of_words) >> std::ws),
                      std::istream_iterator<std::string>());
    std::vector<T> vector_of_T(number_of_words);
    std::stringstream stream(string_of_words);
    size_t i = 0;
    while (i < vector_of_T.size()) {
        stream >> vector_of_T[i++];
    }
    return vector_of_T;
} // space separated string to vector

template <class T> std::string vec2words(const std::vector<T> &v) {
    std::ostringstream o;
    if (!v.empty()) {
        o << v.front();
        for (size_t i = 1; i < v.size(); i++) {
            o << " " << v[i];
        }
    }
    return o.str();
} // vector to space separated string w values

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] Text manipulation") {
    CHECK(vec2words<double>({1.0, -1.2, 0}) == "1 -1.2 0");
    CHECK(words2vec<double>("1 -1.2 0") == std::vector<double>({1.0, -1.2, 0}));
}
#endif

template <typename T> struct BasePointerVector {
    std::vector<std::shared_ptr<T>> vec; //!< Vector of shared pointers to base class

    auto begin() noexcept { return vec.begin(); }
    auto begin() const noexcept { return vec.begin(); }
    auto end() noexcept { return vec.end(); }
    auto end() const noexcept { return vec.end(); }
    auto empty() const noexcept { return vec.empty(); }
    auto size() const noexcept { return vec.size(); }
    auto &back() noexcept { return vec.back(); }
    auto &back() const noexcept { return vec.back(); }
    auto &front() noexcept { return vec.front(); }
    auto &front() const noexcept { return vec.front(); }
    auto &at(size_t n) { return vec.at(n); }
    auto &at(size_t n) const { return vec.at(n); }

    template <typename Tderived, class... Args, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    void emplace_back(Args &... args) {
        vec.push_back(std::make_shared<Tderived>(args...));
    } //!< Create an (derived) instance and append a pointer to it to the vector

    template <typename Tderived, class... Args, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    void emplace_back(const Args &... args) {
        vec.push_back(std::make_shared<Tderived>(args...));
    } //!< Create an (derived) instance and append a pointer to it to the vector

    template <typename Tderived, class Arg, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    void push_back(std::shared_ptr<Arg> arg) {
        vec.push_back(arg);
    } //!< Append a pointer to a (derived) instance to the vector

    template <typename Tderived, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    auto &operator=(const BasePointerVector<Tderived> &d) {
        vec.assign(d.vec.begin(), d.vec.end());
        return *this;
    } //!< Allow assignment to a vector of ancestors

    template <typename Tderived, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>> auto find() const {
        BasePointerVector<Tderived> _v;
        for (auto base : vec) {
            if (auto derived = std::dynamic_pointer_cast<Tderived>(base); derived) {
                _v.template push_back<Tderived>(derived);
            }
        }
        return _v;
    } //!< Pointer list to all matching type

    template <typename U> friend void to_json(json &, const BasePointerVector<U> &); //!< Allow serialization to JSON
}; //!< Helper class for storing vectors of base pointers

template <typename T> void to_json(json &j, const BasePointerVector<T> &b) {
    try {
        for (auto shared_ptr : b.vec) {
            j.push_back(*shared_ptr);
        }
    } catch (const std::exception &e) {
        throw std::runtime_error("error converting to json: "s + e.what());
    }
}

template <typename T> void from_json(const json &j, BasePointerVector<T> &b) {
    try {
        for (auto it : j) {
            std::shared_ptr<T> ptr = it;
            b.template push_back<T>(ptr);
        }
    } catch (const std::exception &e) {
        throw std::runtime_error("error converting from json: "s + e.what());
    }
}
} // namespace Faunus
