#pragma once
#include <nlohmann/json.hpp>
#include <Eigen/Geometry>
#include "units.h"
#include "tensor.h"

namespace Faunus {

using doctest::Approx;

TEST_SUITE_BEGIN("Tensor");

TEST_CASE("[Faunus] Tensor") {
    using doctest::Approx;
    Tensor Q1 = Tensor(1, 2, 3, 4, 5, 6);
    Tensor Q2 = json(Q1);        // Q1 --> json --> Q2
    CHECK(json(Q1) == json(Q2)); // Q1 --> json == json <-- Q2 ?
    CHECK(Q2 == Tensor(1, 2, 3, 4, 5, 6));

    auto m = Eigen::AngleAxisd(pc::pi / 2, Point(0, 1, 0)).toRotationMatrix();
    Q1.rotate(m);
    CHECK(Q1(0, 0) == Approx(6));
    CHECK(Q1(0, 1) == Approx(5));
    CHECK(Q1(0, 2) == Approx(-3));
    CHECK(Q1(1, 1) == Approx(4));
    CHECK(Q1(1, 2) == Approx(-2));
    CHECK(Q1(2, 2) == Approx(1));
}
TEST_SUITE_END();
} // namespace Faunus
