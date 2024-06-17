#pragma once
#include <doctest/doctest.h>
#include <type_traits>
#include <cmath>

namespace Faunus {

/**
 * @brief Approximation of erfc-function
 * @param x Value for which erfc should be calculated
 * @details Reference for this approximation is found in Abramowitz and Stegun,
 *          Handbook of mathematical functions, eq. 7.1.26
 *
 * @f[
 *     \erf(x) = 1 - (a_1t + a_2t^2 + a_3t^3 + a_4t^4 + a_5t^5)e^{-x^2} + \epsilon(x)
 * @f]
 * @f[
 *     t = \frac{1}{1 + px}
 * @f]
 * @f[
 *     |\epsilon(x)| \le 1.5\times 10^{-7}
 * @f]
 *
 * @warning Needs modification if x < 0
 */
template <std::floating_point T> inline T erfc_x(T x) {
    T t = 1.0 / (1.0 + 0.3275911 * x);
    const T a1 = 0.254829592;
    const T a2 = -0.284496736;
    const T a3 = 1.421413741;
    const T a4 = -1.453152027;
    const T a5 = 1.061405429;
    return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))) * std::exp(-x * x);
}

TEST_CASE("[Faunus] erfc_x") {
    double infty = std::numeric_limits<double>::infinity();
    using doctest::Approx;
    CHECK_EQ(erfc_x(infty), Approx(0));
    CHECK_EQ(std::erfc(-infty), Approx(2 - erfc_x(infty)));
    CHECK_EQ(erfc_x(0.0), Approx(1.0));
    CHECK_EQ(2 - erfc_x(0.2), Approx(std::erfc(-0.2)));
    CHECK_EQ(erfc_x(0.2), Approx(std::erfc(0.2)));
}

/**
 * @brief Approximate 1 - erfc_x
 * @param x Value for which erf should be calculated
 */
template <typename T> T inline erf_x(T x) { return (1 - erfc_x(x)); }

} // namespace Faunus