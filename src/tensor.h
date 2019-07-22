#pragma once

#include <nlohmann/json.hpp>
#include <Eigen/Geometry>
#ifdef DOCTEST_LIBRARY_INCLUDED
#include "units.h"
#endif

namespace Faunus {
/**
 * @brief Tensor class
 * @todo add documentation
 */
struct Tensor : public Eigen::Matrix3d {
    typedef Eigen::Matrix3d base;

    Tensor(); //!< Constructor, clear data

    /**
     * @brief Constructor
     * @todo explain input
     */
    Tensor(double xx, double xy, double xz, double yy, double yz, double zz); //!< Construct from input

    void rotate(const base &m); //!< Rotate using rotation matrix. Remove?

    void eye();

    template <typename T> Tensor(const Eigen::MatrixBase<T> &other) : base(other) {}

    template <typename T> Tensor &operator=(const Eigen::MatrixBase<T> &other) {
        base::operator=(other);
        return *this;
    }
}; //!< Tensor class

void to_json(nlohmann::json &j, const Tensor &t);   //!< Tensor -> Json
void from_json(const nlohmann::json &j, Tensor &t); //!< Json -> Tensor

#ifdef DOCTEST_LIBRARY_INCLUDED
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
#endif

} // namespace Faunus
