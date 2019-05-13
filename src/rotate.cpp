#include "rotate.h"
#include <Eigen/Geometry>

namespace Faunus {
QuaternionRotate::QuaternionRotate(double angle, Point u) { set(angle, u); }

void QuaternionRotate::set(double angle, Point u) {
    this->angle = angle;
    u.normalize();
    first = Eigen::AngleAxisd(angle, u);
    second << 0, -u.z(), u.y(), u.z(), 0, -u.x(), -u.y(), u.x(), 0;
    second = Eigen::Matrix3d::Identity() + second * std::sin(angle) + second * second * (1 - std::cos(angle));

    // Quaternion can be converted to rotation matrix:
    // second = first.toRotationMatrix()
}

Point QuaternionRotate::operator()(Point a, std::function<void(Point &)> boundary, const Point &shift) const {
    a = a - shift;
    boundary(a);
    a = first * a + shift;
    boundary(a);
    return a;
    // https://www.cc.gatech.edu/classes/AY2015/cs4496_spring/Eigen.html
}

auto QuaternionRotate::operator()(const Eigen::Matrix3d &a) const { return second * a * second.transpose(); }
} // namespace Faunus
