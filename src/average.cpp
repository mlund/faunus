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

} // namespace Faunus