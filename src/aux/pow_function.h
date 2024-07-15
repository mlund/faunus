#pragma once
#include <doctest/doctest.h>

namespace Faunus {
/**
 * @brief n'th integer power of float
 *
 * On GCC/Clang this will use the fast `__builtin_powi` function.
 *
 * See also:
 * - https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp
 * - https://martin.ankerl.com/2007/10/04/optimized-pow-approximation-for-java-and-c-c/
 */
template <class T> inline constexpr T powi(T x, unsigned int n)
{
#if defined(__GNUG__)
    return __builtin_powi(x, n);
#else
    return n > 0 ? x * powi(x, n - 1) : 1;
#endif
}

TEST_CASE("[Faunus] powi")
{
    using doctest::Approx;
    double x = 3.1;
    CHECK_EQ(powi(x, 0), Approx(1));
    CHECK_EQ(powi(x, 1), Approx(x));
    CHECK_EQ(powi(x, 2), Approx(x * x));
    CHECK_EQ(powi(x, 4), Approx(x * x * x * x));
}
} // namespace Faunus
