#include <doctest/doctest.h>
#include "tensor.h"
#include "units.h"
#include "core.h"
#include <Eigen/Geometry>
#include <nlohmann/json.hpp>

namespace Faunus {

Tensor::Tensor() { base::setZero(); }

Tensor::Tensor(double xx, double xy, double xz, double yy, double yz, double zz) {
    (*this) << xx, xy, xz, xy, yy, yz, xz, yz, zz;
}

void Tensor::rotate(const Tensor::base &rotation_matrix) {
    (*this) = rotation_matrix * (*this) * rotation_matrix.transpose();
}

void Tensor::eye() { *this = base::Identity(3, 3); }

void to_json(nlohmann::json &j, const Tensor &t) { j = {t(0, 0), t(0, 1), t(0, 2), t(1, 1), t(1, 2), t(2, 2)}; }

void from_json(const nlohmann::json &j, Tensor &t) {
    if (j.size() == 6 and j.is_array())
        t = Tensor(j[0], j[1], j[2], j[3], j[4], j[5]);
    else
        throw std::runtime_error("tensor expects array with exactly six coefficients");
}

TEST_SUITE_BEGIN("Tensor");

TEST_CASE("[Faunus] Tensor") {
    using doctest::Approx;
    using namespace nlohmann;
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
