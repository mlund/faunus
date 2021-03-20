#pragma once
#include <array>
namespace Faunus {
/**
 * @brief Evaluate n'th degree Legendre polynomial
 *
 * Example:
 * @code
 * Legendre<float, 10> l;
 * auto P = l.eval(1.3)
 * std::cout << P[3]; --> third order value
 * @endcode
 *
 * @author Mikael Lund
 * @date Canberra 2005-2020
 * @note Since C++17 there's `std::legendre` but this seems more efficient
 *       if a range of degrees are needed
 * @tparam T floating point type (double, float, ...)
 * @tparam max_order Maximum order to evaluate
 * @tparam use_table Use lookup table for 1+1/i? Default = false
 * @todo Benchmark `use_table`
 */
template <typename T, std::size_t max_order, bool use_table = false> class Legendre {
  private:
    std::array<T, max_order + 1> y; //!< Lookup table for 1+1/i (overkill?)
    std::array<T, max_order + 1> P; //!< Legendre terms stored here
  public:
    Legendre() {
        P[0] = 1.0;
        if constexpr (use_table) {
            for (std::size_t i = 1; i < max_order; ++i) {
                y[i] = 1.0 + 1.0 / T(i);
            }
        }
    }

    /** @brief Evaluate polynomials at x */
    const auto &eval(T x) {
        if constexpr (max_order > 0) {
            P[1] = x;
            for (std::size_t i = 1; i < max_order; ++i) {
                if constexpr (use_table) {
                    P[i + 1] = ((y[i] + 1.0) * x * P[i] - P[i - 1]) / y[i];
                } else {
                    P[i + 1] = ((2.0 + 1.0 / T(i)) * x * P[i] - P[i - 1]) / (1.0 + 1.0 / T(i));
                }
            }
        }
        return P;
    }
};
#ifdef DOCTEST_LIBRARY_INCLUDED__
TEST_CASE_TEMPLATE("[Faunus] Legendre", LegendreType, Legendre<double, 3, false>, Legendre<double, 3, true>) {
    using doctest::Approx;
    LegendreType l;
    double x = 2.2;
    auto P = l.eval(x);
    CHECK(P[0] == Approx(1));
    CHECK(P[1] == Approx(x));
    CHECK(P[2] == Approx(0.5 * (3 * x * x - 1)));
    CHECK(P[3] == Approx(0.5 * (5 * x * x * x - 3 * x)));
}
#endif

} // namespace Faunus