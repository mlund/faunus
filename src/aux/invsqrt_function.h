#pragma once
#include <doctest/doctest.h>
#include <vector>
#include <type_traits>
#include <cmath>
#include <cstdint>

namespace Faunus {

/**
 * @brief Fast inverse square-root approximation
 *
 * This is Carmack et al's "magic" inverse square root routine, modified to
 * work with both double and float and with one (less precise) or two (more precise)
 * iterations. Conditionals and sanity checks are evaluated at compile time.
 *
 * @tparam float_t Floating point type (float or double)
 * @tparam iterations Number of iterations (1 or 2)
 * @param x value to operate on
 * @returns Approximate inverse square root of x
 * @remarks Code comments supposedly from the original Quake III Arena code
 * @see https://en.wikipedia.org/wiki/Fast_inverse_square_root
 */
template <typename float_t = double, char iterations = 1> inline float_t inv_sqrt(float_t x)
{
    static_assert(iterations == 1 or iterations == 2);
    static_assert(std::is_floating_point<float_t>::value);
    typedef typename std::conditional<sizeof(float_t) == 8, std::int64_t, std::int32_t>::type int_t;
    static_assert(sizeof(float_t) == sizeof(int_t));
    float_t x_half = x * 0.5;
    float_t y = x;
    int_t i = *reinterpret_cast<int_t*>(&y); // evil floating point bit level hacking
    i = (sizeof(float_t) == 8 ? 0x5fe6eb50c7b537a9 : 0x5f3759df) - (i >> 1); // what the fuck?
    y = *reinterpret_cast<float_t*>(&i);
    y = y * (1.5 - x_half * y * y); // 1st iteration
    if constexpr (iterations == 2) {
        y = y * (1.5 - x_half * y * y); // 2nd iteration, this can be removed
    }
    return y;
}
} // namespace Faunus

TEST_CASE_TEMPLATE("inv_sqrt", T, double, float)
{
    std::vector<T> vals = {0.23, 3.3, 10.2, 100.45, 512.06};
    for (auto x : vals) {
        CHECK_EQ(Faunus::inv_sqrt<T, 2>(x), doctest::Approx(1.0 / std::sqrt(x)));
    }
}
