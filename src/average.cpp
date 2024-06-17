#include <doctest/doctest.h>
#include "average.h"

namespace Faunus {

TEST_CASE("[Faunus] Average")
{
    Faunus::Average<double> a;
    a += 1.0;
    a += 2.0;
    CHECK_EQ(a.size(), 2);
    CHECK_EQ(a.avg(), 1.5);
    CHECK_EQ(static_cast<double>(a), 1.5); // conversion to double

    auto b = a;        // copy
    CHECK(!b.empty()); // check not empty()
    CHECK_EQ(a, b);
    CHECK_EQ(a, b);
    b.clear();        // reset all data
    CHECK(b.empty()); // check empty()
    b += 2.0;
    b += 3.0;
    CHECK((a < b)); // a.avg() < b.avg()
    CHECK_EQ((a + b).avg(), doctest::Approx(2.0));

    b = 1.0; // assign from double
    CHECK_EQ(b.size(), 1);

    SUBCASE("overflow")
    {
        Average<double, std::uint8_t> a;
        const auto max_number_of_samples = std::numeric_limits<std::uint8_t>::max();
        for (std::uint8_t i = 0; i < max_number_of_samples; i++) {
            a += 1.0;
        }
        CHECK_EQ(a.size(), max_number_of_samples);
        CHECK_THROWS(a += 1.0);
    }
}

TEST_CASE("[Faunus] AverageStd")
{
    Faunus::AverageStdev<double> a;
    a += 1.0;
    a += 2.0;
    CHECK_EQ(a.size(), 2);
    CHECK_EQ(a.rms(), doctest::Approx(1.5811388301));
    CHECK_EQ(a.avg(), 1.5);
    CHECK_EQ(static_cast<double>(a), 1.5); // conversion to double

    auto b = a;        // copy
    CHECK(!b.empty()); // check not empty()
    CHECK_EQ(a, b);
    CHECK_EQ(a, b);
    b.clear();        // reset all data
    CHECK(b.empty()); // check empty()
    b += 2.0;
    b += 3.0;
    CHECK((a < b)); // a.avg() < b.avg()

    b = 1.0; // assign from double
    CHECK_EQ(b.size(), 1);
}

TEST_CASE("[Faunus] AverageObj")
{
    using doctest::Approx;
    Faunus::AverageObj<double> a;
    a += 1.0;
    a += 2.0;
    CHECK_EQ(a.avg(), Approx(1.5));
    CHECK_EQ(a, Approx(1.5)); // implicit conversion to double

    struct MyClass
    {
        double x;

        MyClass& operator+=(const MyClass& other)
        {
            x += other.x;
            return *this;
        } // required

        MyClass operator*(double value) const
        {
            MyClass scaled;
            scaled.x = x * value;
            return scaled;
        } // required
    };

    Faunus::AverageObj<MyClass> b;
    MyClass k1, k2;
    k1.x = 10.0;
    k2.x = 20.0;
    b += k1;
    b += k2;
    CHECK_EQ(b.avg().x, Approx(15.0));
}

} // namespace Faunus