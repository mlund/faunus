#pragma once

#include <Eigen/Geometry>

namespace Faunus {

/**
 * @brief Rotation routine using the Eigen library
 * @todo Get rotate tensors using quaternion
 *
 * This can rotate vectors and tensors using quaternion and rotation matrix (for matrix rotation).
 */
class QuaternionRotate {
  private:
    using Point = Eigen::Vector3d;
    Eigen::Quaterniond quaternion;
    Eigen::Matrix3d rotation_matrix;

  public:
    double angle = 0.0; //!< Current rotation angle
    QuaternionRotate() = default;
    QuaternionRotate(double angle, Point axis);           //!< Set angle and rotation axis
    void set(double angle, Point axis);                   //!< Set angle and rotation axis
    auto operator()(const Eigen::Matrix3d &matrix) const; //!< Rotate matrix/tensor
    const Eigen::Quaterniond &getQuaternion() const;      //!< Get current quaternion
    const Eigen::Matrix3d &getRotationMatrix() const;     //!< Get current rotation matrix
    Point operator()(
        Point vector, std::function<void(Point &)> boundary = [](Point &) {},
        const Point &shift = {0, 0, 0}) const; //!< Rotate point w. optional PBC boundaries
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
        a = qrot.getRotationMatrix() * a;
        CHECK(a.x() == Approx(0));
        a = qrot.getRotationMatrix() * a;
        CHECK(a.x() == Approx(-1));
    }
}
#endif

} // namespace Faunus
