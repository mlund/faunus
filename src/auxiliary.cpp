#include <doctest/doctest.h>
#include "auxiliary.h"
#include <string>
#include <cmath>

namespace Faunus {

using doctest::Approx;

TEST_SUITE_BEGIN("Auxiliary");

TEST_CASE("[Faunus] for_each_pair") {
    int x;
    std::vector<int> a = {1, 2, 3};
    x = for_each_unique_pair(
        a.begin(), a.end(), [](int i, int j) { return i * j; }, std::plus<>());
    CHECK_EQ(x, 2 + 3 + 6);
    a.resize(1);
    x = for_each_unique_pair(
        a.begin(), a.end(), [](int i, int j) { return i * j; }, std::plus<>());
    CHECK_EQ(x, 0);
}

TEST_CASE("[Faunus] ordered_pair") {
    ordered_pair<int> a = {1, 2}, b = {2, 1};
    CHECK(((a.first == 1 && a.second == 2)));
    CHECK(((b.first == 1 && b.second == 2)));
    CHECK_EQ(a, b);
    CHECK(a.contains(1));
    CHECK(a.contains(2));
    CHECK(!a.contains(3));
}

TEST_CASE("[Faunus] Text manipulation") {
    CHECK_EQ(joinToString(std::vector<double>{1.0, -1.2, 0}), "1 -1.2 0");
    CHECK_EQ(splitConvert<double>("1 -1.2 0"), std::vector<double>{1.0, -1.2, 0});
}

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
