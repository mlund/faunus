#include <doctest/doctest.h>
#include "average.h"

namespace Faunus {

TEST_CASE("[Faunus] Average") {
    Faunus::Average<double> a;
    a += 1.0;
    a += 2.0;
    CHECK(a.cnt == 2);
    CHECK(a.sum == 3);
    CHECK(a.sqsum == 5);
    CHECK(a.avg() == 1.5);
    CHECK(a == 1.5); // implicit conversion to double

    auto b = a;        // copy
    CHECK(!b.empty()); // check not empty()
    CHECK(a == b);
    b.clear();        // reset all data
    CHECK(b.empty()); // check empty()
    b += 2.0;
    b += 3.0;
    CHECK(a < b); // a.avg() < b.avg()
    CHECK((a + b).avg() == doctest::Approx(2));

    b = 1.0; // assign from double
    CHECK(b.size() == 1);
}

TEST_CASE("[Faunus] AverageObj") {
    using doctest::Approx;
    Faunus::AverageObj<double> a;
    a += 1.0;
    a += 2.0;
    CHECK(a.avg() == Approx(1.5));
    CHECK(a == Approx(1.5)); // implicit conversion to double

    struct MyClass {
        double x;
        MyClass &operator+=(const MyClass &other) {
            x += other.x;
            return *this;
        } // required
        MyClass operator*(double value) const {
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
    CHECK(b.avg().x == Approx(15.0));
}

} // namespace Faunus