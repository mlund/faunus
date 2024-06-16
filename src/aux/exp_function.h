#pragma once
#include <doctest/doctest.h>
#include <limits>
#include <cstdint>

namespace Faunus {
/**
 * @brief Approximate exp() function
 * @note see [Cawley 2000](http://dx.doi.org/10.1162/089976600300015033)
 * @warning Does not work in big endian systems, nor on gcc
 *
 * Update 2019: http://www.federicoperini.info/wp-content/uploads/FastExp-Final.pdf
 */
template <class Tint = std::int32_t> double exp_cawley(double y)
{
    static_assert(2 * sizeof(Tint) == sizeof(double), "Approximate exp() requires 4-byte integer");
    union {
        double d;
        struct {
            Tint j, i;
        } n; // little endian
        // struct { int i, j; } n;  // bin endian
    } eco;
    eco.n.i = 1072632447 + (Tint)(y * 1512775.39519519);
    eco.n.j = 0;
    return eco.d;
}

inline double exp_untested(double y)
{
    typedef std::int32_t Tint;
    static_assert(2 * sizeof(Tint) == sizeof(double), "Approximate exp() requires 4-byte integer");
    double d(0);
    *((Tint*)(&d) + 0) = 0;
    *((Tint*)(&d) + 1) = (Tint)(1512775 * y + 1072632447);
    return d;
}

TEST_CASE("[Faunus] exp_cawley")
{
    double infty = std::numeric_limits<double>::infinity();
    using doctest::Approx;
    WARN(exp_cawley(-infty) == Approx(0)); // clang=OK; GCC=not OK
    CHECK_EQ(exp_cawley(0), Approx(0.9710078239));
    CHECK_EQ(exp_cawley(2), Approx(7.3096199036));
    CHECK_EQ(exp_cawley(-2), Approx(0.13207829));
}

TEST_CASE("[Faunus] exp_untested")
{
    double infty = std::numeric_limits<double>::infinity();
    using doctest::Approx;
    CHECK_EQ(exp_untested(-infty), Approx(0)); // clang=OK; GCC=not OK
    CHECK_EQ(exp_untested(0), Approx(0.9710078239));
    CHECK_EQ(exp_untested(2), Approx(7.3096199036));
    CHECK_EQ(exp_untested(-2), Approx(0.13207829));
}

} // namespace Faunus
