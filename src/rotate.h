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
    QuaternionRotate(double angle, Point axis);                     //!< Set angle and rotation axis
    void set(double angle, Point axis);                             //!< Set angle and rotation axis
    auto operator()(const Eigen::Matrix3d& matrix) const;           //!< Rotate matrix/tensor
    [[nodiscard]] const Eigen::Quaterniond& getQuaternion() const;  //!< Get current quaternion
    [[nodiscard]] const Eigen::Matrix3d& getRotationMatrix() const; //!< Get current rotation matrix
    Point operator()(
        Point vector, std::function<void(Point&)> boundary = [](Point&) {},
        const Point& shift = {0, 0, 0}) const; //!< Rotate point w. optional PBC boundaries
};

} // namespace Faunus
