#include "auxiliary.h"
#include <cmath>

namespace Faunus {

using doctest::Approx;

TEST_SUITE_BEGIN("Auxiliary");

TEST_CASE("numeric_cast") {
    CHECK_EQ(numeric_cast<int>(2.2), 2);
    CHECK_EQ(numeric_cast<int>(-2.2), -2);
    CHECK_EQ(numeric_cast<int>(2.8), 3);
    CHECK_EQ(numeric_cast<int>(-2.8), -3);

    double short_overflow = std::numeric_limits<short>::max() + 1;
    double short_underflow = std::numeric_limits<short>::min() - 1;
    CHECK_NOTHROW(numeric_cast<int>(short_overflow));
    CHECK_THROWS_AS(numeric_cast<short>(short_overflow), std::overflow_error);
    CHECK_NOTHROW(numeric_cast<int>(short_underflow));
    CHECK_THROWS_AS(numeric_cast<short>(short_underflow), std::overflow_error);
    CHECK_NOTHROW(numeric_cast<int>(-1.0));
    CHECK_THROWS_AS(numeric_cast<unsigned>(-1.0), std::overflow_error);
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdiv-by-zero"
    CHECK_THROWS_AS(numeric_cast<int>(1. / 0), std::overflow_error);
    #pragma GCC diagnostic pop
}

TEST_SUITE_END();
} // namespace Faunus
