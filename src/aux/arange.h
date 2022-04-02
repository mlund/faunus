/*
The MIT License (MIT)
Copyright © 2021 Mikael Lund

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the “Software”), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.
THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once
#include <cmath>
#include <concepts>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

namespace Faunus {

/**
 * @brief Return evenly spaced values within a given interval
 *
 * Values are lazily generated within the half-open interval [start, stop)
 * In other words, the interval including start but excluding stop.
 * This function has the same behavior as Numpy's `arange()`
 *
 * Examples:
 *
 *     arange(4.0, 10.0, 1.0); // --> 4 5 6 7 8 9
 *     arange(4.0, 20.0, 3.0); // --> 4 7 10 13 16 19
 *     arange(-1.0, 1.0, 0.5); // --> -1 -0.5 0 0.5
 *
 * @tparam T Value type; must be floating point or integral
 * @param start Start of interval (included)
 * @param stop End of interval (excluded)
 * @param step Spacing between values
 * @return Range of lazily generated values
 */
template <typename T> constexpr auto arange(const T start, const T stop, const T step) {
    static_assert(std::is_floating_point_v<T> || std::is_integral_v<T>, "floating point or integral type required");
    using float_type = typename std::conditional<std::is_floating_point_v<T>, T, long double>::type;
    using int_type = typename std::conditional<std::is_integral_v<T>, T, int>::type;
    const auto length = static_cast<int_type>(std::ceil((stop - start) / static_cast<float_type>(step)));
    return ranges::cpp20::views::iota(int_type(0), length) |
           ranges::cpp20::views::transform([start, step](auto i) -> T { return start + static_cast<T>(i) * step; });
}

TEST_CASE("[Faunus] arange") {
    SUBCASE("Step = 1 (float)") {
        auto r = arange(4.0, 10.0, 1.0); // --> 4 5 6 7 8 9
        auto pos = r.begin();
        CHECK(ranges::size(r) == 6);
        CHECK(*(pos++) == doctest::Approx(4.0));
        CHECK(*(pos++) == doctest::Approx(5.0));
        CHECK(*(pos++) == doctest::Approx(6.0));
        CHECK(*(pos++) == doctest::Approx(7.0));
        CHECK(*(pos++) == doctest::Approx(8.0));
        CHECK(*(pos++) == doctest::Approx(9.0));
    }
    SUBCASE("Step = 1 (int)") {
        auto r = arange(4, 10, 1); // --> 4 5 6 7 8 9
        auto pos = r.begin();
        CHECK(ranges::size(r) == 6);
        CHECK(*(pos++) == 4);
        CHECK(*(pos++) == 5);
        CHECK(*(pos++) == 6);
        CHECK(*(pos++) == 7);
        CHECK(*(pos++) == 8);
        CHECK(*(pos++) == 9);
    }
    SUBCASE("Step > 1 (float)") {
        auto r = arange(4.0, 20.0, 3.0); // --> 4 7 10 13 16 19
        auto pos = r.begin();
        CHECK(ranges::size(r) == 6);
        CHECK(*(pos++) == doctest::Approx(4.0));
        CHECK(*(pos++) == doctest::Approx(7.0));
        CHECK(*(pos++) == doctest::Approx(10.0));
        CHECK(*(pos++) == doctest::Approx(13.0));
        CHECK(*(pos++) == doctest::Approx(16.0));
        CHECK(*(pos++) == doctest::Approx(19.0));
    }
    SUBCASE("Step > 1 (int)") {
        auto r = arange(4, 20, 3); // --> 4 7 10 13 16 19
        auto pos = r.begin();
        CHECK(ranges::size(r) == 6);
        CHECK(*(pos++) == 4);
        CHECK(*(pos++) == 7);
        CHECK(*(pos++) == 10);
        CHECK(*(pos++) == 13);
        CHECK(*(pos++) == 16);
        CHECK(*(pos++) == 19);
    }

    SUBCASE("Step < 1") {
        auto r = arange(-1.0, 1.0, 0.5); // --> -1 -0.5 0 0.5
        auto pos = r.begin();
        CHECK(ranges::size(r) == 4);
        CHECK(*(pos++) == doctest::Approx(-1.0));
        CHECK(*(pos++) == doctest::Approx(-0.5));
        CHECK(*(pos++) == doctest::Approx(0.0));
        CHECK(*(pos++) == doctest::Approx(0.5));
    }
}

} // namespace Faunus
