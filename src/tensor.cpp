#include "tensor.h"

namespace Faunus {

Tensor::Tensor() { base::setZero(); }

Tensor::Tensor(double xx, double xy, double xz, double yy, double yz, double zz) {
    (*this) << xx, xy, xz, xy, yy, yz, xz, yz, zz;
}

void Tensor::rotate(const Tensor::base &m) { (*this) = m * (*this) * m.transpose(); }

void Tensor::eye() { *this = base::Identity(3, 3); }

void to_json(nlohmann::json &j, const Tensor &t) { j = {t(0, 0), t(0, 1), t(0, 2), t(1, 1), t(1, 2), t(2, 2)}; }

void from_json(const nlohmann::json &j, Tensor &t) {
    if (j.size() != 6 || !j.is_array())
        throw std::runtime_error("Json->Tensor: array w. exactly six coefficients expected.");
    t = Tensor(j[0], j[1], j[2], j[3], j[4], j[5]);
}

} // namespace Faunus
