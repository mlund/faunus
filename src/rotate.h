#pragma once

#include <Eigen/Geometry>

namespace Faunus {

/* @brief Quaternion rotation routine using the Eigen library
 * */
struct QuaternionRotate : public std::pair<Eigen::Quaterniond, Eigen::Matrix3d> {

    typedef Eigen::Vector3d Point;
    typedef std::pair<Eigen::Quaterniond, Eigen::Matrix3d> base;
    using base::first;
    using base::second;

    double angle = 0; //!< Rotation angle

    QuaternionRotate() = default;
    QuaternionRotate(double angle, Point u);
    void set(double angle, Point u);
    Point operator()(Point a, std::function<void(Point &)> boundary = [](Point &) {},
                     const Point &shift = {0, 0, 0}) const; //!< Rotate point w. optional PBC boundaries
    auto operator()(const Eigen::Matrix3d &a) const;        //!< Rotate matrix/tensor
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] QuaternionRotate") {
    using doctest::Approx;
    QuaternionRotate qrot;
    Point a = {1, 0, 0};
    qrot.set(pc::pi / 2, {0, 1, 0}); // rotate around y-axis
    CHECK(qrot.angle == Approx(pc::pi / 2));

    SUBCASE("rotate using quaternion") {
        a = qrot(a); // rot. 90 deg.
        CHECK(a.x() == Approx(0));
        a = qrot(a); // rot. 90 deg.
        CHECK(a.x() == Approx(-1));
    }

    SUBCASE("rotate using rotation matrix") {
        a = qrot.second * a;
        CHECK(a.x() == Approx(0));
        a = qrot.second * a;
        CHECK(a.x() == Approx(-1));
    }
}
#endif

} // namespace Faunus
