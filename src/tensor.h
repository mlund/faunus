#pragma once
#include <Eigen/Core>
#include <nlohmann/json.hpp>

namespace Faunus {
/**
 * @brief Tensor class
 * @todo add documentation
 */
struct Tensor : public Eigen::Matrix3d
{
    typedef Eigen::Matrix3d base;

    Tensor(); //!< Constructor, clear data

    /**
     * @brief Constructor
     * @todo explain input
     */
    Tensor(double, double, double, double, double, double); //!< Constructor

    void rotate(const base&); //!< Rotate using rotation matrix

    [[maybe_unused]] [[maybe_unused]] void eye();

    template <typename T>
    Tensor(const Eigen::MatrixBase<T>& other)
        : base(other)
    {
    }

    template <typename T> Tensor& operator=(const Eigen::MatrixBase<T>& other)
    {
        base::operator=(other);
        return *this;
    }
}; //!< Tensor class

void to_json(nlohmann::json&, const Tensor&);   //!< Tensor -> Json
void from_json(const nlohmann::json&, Tensor&); //!< Json -> Tensor
} // namespace Faunus
