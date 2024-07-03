#pragma once
#include "average.h"
#include <nlohmann/json.hpp>
#include <ranges>
#include <functional>
#include <iterator>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <concepts>

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
template <std::integral TOut, std::floating_point TIn> inline TOut numeric_cast(const TIn number)
{
    if (std::isfinite(number)) {
        // The number is finite ...
        if (number < std::nextafter(static_cast<TIn>(std::numeric_limits<TOut>::max()), 0) &&
            number > std::nextafter(static_cast<TIn>(std::numeric_limits<TOut>::min()), 0)) {
            // ... and fits into the integral type range.
            // The nextafter function is used to mitigate possible rounding up or down.
            return static_cast<TOut>(number >= 0 ? number + TIn(0.5)
                                                 : number - TIn(0.5)); // round before cast
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
template <std::forward_iterator Titer, typename Tfunction, typename T = double,
          typename Taggregate_function>
T for_each_unique_pair(Titer begin, Titer end, Tfunction f,
                       Taggregate_function aggregator = std::plus<T>())
{
    T x = T();
    for (auto i = begin; i != end; ++i) {
        for (auto j = i; ++j != end;) {
            x = aggregator(x, f(*i, *j));
        }
    }
    return x;
}

/** @brief Erase from `target` range all values found in `values` range */
template <std::ranges::range T> T erase_range(T target, const T& values)
{
    target.erase(std::remove_if(target.begin(), target.end(),
                                [&](auto i) {
                                    return std::find(values.begin(), values.end(), i) !=
                                           values.end();
                                }),
                 target.end());
    return target;
}

/**
 * @brief Ordered pair where `first<=second`
 *
 * Upon construction, the smallest element is placed in `first`
 * so that `ordered_pair<int>(i,j)==ordered_pair<int>(j,i)` is always true.
 *
 * @todo Add std::pair copy operator
 */
template <class T> struct ordered_pair : public std::pair<T, T>
{
    using base = std::pair<T, T>;
    ordered_pair() = default;

    ordered_pair(const T& a, const T& b)
        : base(std::minmax(a, b))
    {
    }

    bool contains(const T& value) const { return value == base::first || value == base::second; }
};

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
template <std::floating_point T, std::integral Tint = int> Tint to_bin(T x, T dx = 1)
{
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
template <std::floating_point Tfloat = double> class Quantize
{
  private:
    Tfloat xmin, dx, x;

  public:
    /**
     * @brief Constructor
     * @param dx resolution
     * @param xmin minimum value if converting to integral type (binning)
     */
    Quantize(Tfloat dx, Tfloat xmin = 0)
        : xmin(xmin)
        , dx(dx)
    {
    }

    /** @brief Assigment operator */
    Quantize& operator=(Tfloat val)
    {
        assert(val >= xmin);
        if (val >= 0) {
            x = int(val / dx + 0.5) * dx;
        }
        else {
            x = int(val / dx - 0.5) * dx;
        }
        assert(x >= xmin);
        return *this;
    }

    Quantize& frombin(unsigned int i)
    {
        x = i * dx + xmin;
        return *this;
    }

    /** @brief Assignment with function operator */
    Quantize& operator()(Tfloat val)
    {
        *this = val;
        return *this;
    }

    /** @brief Implicit convertion to integral (bin) or float (rounded) */
    template <typename T> operator T()
    {
        if (std::is_integral<T>::value) {
            return T((x - xmin) / dx + 0.5);
        }
        return x;
    }
};

template <class T>
concept StringStreamable = requires(T a) {
    { std::istringstream() >> a };
    { std::ostringstream() << a };
};

/**
 * @brief Convert space separated words into vector of type T
 *
 * Example:
 *
 *     textio::splitConvert<double>("0.2 1 100"); // -> {0.2, 1.0, 100.0}
 *
 * @returns std::vector of type T
 */
template <StringStreamable T> auto splitConvert(const std::string& words)
{
    auto stream = std::istringstream(words);
    return std::vector<T>(std::istream_iterator<T>(stream), std::istream_iterator<T>());
} // space separated string to vector of values

/**
 * @brief Convert range of values into space separated string
 * @param values Range (vector, set, ...) of values to convert
 * @return String with space sepatated values
 */
template <std::ranges::range Range>
std::string joinToString(const Range& values)
    requires StringStreamable<std::ranges::range_value_t<Range>>
{
    std::ostringstream o;
    if (!values.empty()) {
        o << *values.begin();
        std::for_each(std::next(values.begin()), values.end(), [&](auto& val) { o << " " << val; });
    }
    return o.str();
}

template <typename T> struct BasePointerVector
{
    using value_type = std::shared_ptr<T>;
    std::vector<value_type> vec; //!< Vector of shared pointers to base class

    auto begin() noexcept { return vec.begin(); }

    auto begin() const noexcept { return vec.begin(); }

    auto end() noexcept { return vec.end(); }

    auto end() const noexcept { return vec.end(); }

    auto empty() const noexcept { return vec.empty(); }

    auto size() const noexcept { return vec.size(); }

    auto& back() noexcept { return vec.back(); }

    auto& back() const noexcept { return vec.back(); }

    auto& front() noexcept { return vec.front(); }

    auto& front() const noexcept { return vec.front(); }

    auto& at(size_t n) { return vec.at(n); }

    auto& at(size_t n) const { return vec.at(n); }

    template <typename Tderived, class... Args,
              class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    void emplace_back(Args&... args)
    {
        vec.push_back(std::make_shared<Tderived>(args...));
    } //!< Create an (derived) instance and append a pointer to it to the vector

    template <typename Tderived, class... Args,
              class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    void emplace_back(const Args&... args)
    {
        vec.push_back(std::make_shared<Tderived>(args...));
    } //!< Create an (derived) instance and append a pointer to it to the vector

    //    template <typename Tderived, class Arg, class = std::enable_if_t<std::is_base_of<T,
    //    Tderived>::value>> auto& push_back(std::shared_ptr<Arg> arg) {
    //        vec.push_back(arg);
    //        return vec.back(); // reference to element just added
    //    }                      //!< Append a pointer to a (derived) instance to the vector

    auto& push_back(value_type arg)
    {
        vec.push_back(arg);
        return vec.back(); // reference to element just added
    } //!< Append a pointer to a (derived) instance to the vector

    template <typename Tderived, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    auto& operator=(const BasePointerVector<Tderived>& d)
    {
        vec.assign(d.vec.begin(), d.vec.end());
        return *this;
    } //!< Allow assignment to a vector of ancestors

    template <typename Tderived, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    auto find() const
    {
        BasePointerVector<Tderived> _v;
        for (auto& base : vec) {
            if (auto derived = std::dynamic_pointer_cast<Tderived>(base); derived) {
                _v.push_back(derived);
            }
        }
        return _v;
    } //!< Pointer list to all matching type

    template <typename Tderived, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
    std::shared_ptr<Tderived> findFirstOf() const
    {
        for (auto& base : vec) {
            if (auto derived = std::dynamic_pointer_cast<Tderived>(base); derived) {
                return derived;
            }
        }
        return nullptr;
    }

    template <typename U>
    friend void to_json(nlohmann::json&,
                        const BasePointerVector<U>&); //!< Allow serialization to JSON
}; //!< Helper class for storing vectors of base pointers

template <typename T> void to_json(nlohmann::json& j, const BasePointerVector<T>& b)
{
    using namespace std::string_literals;
    try {
        for (auto& shared_ptr : b.vec) {
            j.push_back(*shared_ptr);
        }
    }
    catch (const std::exception& e) {
        throw std::runtime_error("error converting to json: "s + e.what());
    }
}

template <typename T> void from_json(const nlohmann::json& j, BasePointerVector<T>& b)
{
    using namespace std::string_literals;
    try {
        for (const auto& it : j) {
            std::shared_ptr<T> ptr = it;
            b.push_back(ptr);
        }
    }
    catch (const std::exception& e) {
        throw std::runtime_error("error converting from json: "s + e.what());
    }
}
} // namespace Faunus
