#pragma once
#include <Eigen/Core>
#include <nlohmann/json_fwd.hpp>

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
} // namespace Faunus
