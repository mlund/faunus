#include <doctest/doctest.h>
#include "rotate.h"
#include "core.h"
#include "units.h"

namespace Faunus {

/**
 * @param angle rotation angle (radians)
 * @param axis rotation axis
 */
QuaternionRotate::QuaternionRotate(double angle, Point axis) { set(angle, axis); }

/**
 * @param angle rotation angle (radians)
 * @param axis rotation axis
 */
void QuaternionRotate::set(double angle, Point axis) {
    this->angle = angle;
    axis.normalize(); // make unit vector
    quaternion = Eigen::AngleAxisd(angle, axis);
    rotation_matrix << 0.0, -axis.z(), axis.y(), axis.z(), 0.0, -axis.x(), -axis.y(), axis.x(), 0.0;
    rotation_matrix = Eigen::Matrix3d::Identity() + rotation_matrix * std::sin(angle) +
                      rotation_matrix * rotation_matrix * (1.0 - std::cos(angle));
    // Quaternion can be converted to rotation matrix:
    // second = first.toRotationMatrix()
}

/**
 * @param vector Point to rotate
 * @param boundary Boundary function to handle PBC
 * @param shift Shift used to aid in boundary removal
 * @return Rotated vector
 */
QuaternionRotate::Point QuaternionRotate::operator()(Point vector, std::function<void(Point &)> boundary,
                                                     const Point &shift) const {
    vector = vector - shift;
    boundary(vector);
    vector = quaternion * vector + shift;
    boundary(vector);
    return vector;
    // https://www.cc.gatech.edu/classes/AY2015/cs4496_spring/Eigen.html
}

/**
 * @param matrix Matrix or tensor to rotate
 * @return Rotated matrix or tensor
 */
auto QuaternionRotate::operator()(const Eigen::Matrix3d &matrix) const {
    return rotation_matrix * matrix * rotation_matrix.transpose();
}

const Eigen::Quaterniond &QuaternionRotate::getQuaternion() const { return quaternion; }
const Eigen::Matrix3d &QuaternionRotate::getRotationMatrix() const { return rotation_matrix; }

} // namespace Faunus

using namespace Faunus;

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
