//
// Created by Mikael Lund on 2022-04-07.
//

#ifndef FAUNUS_MATRIXMARKET_H
#define FAUNUS_MATRIXMARKET_H

#include <iostream>
#include <ostream>
#include <spdlog/fmt/fmt.h>
#include <Eigen/SparseCore>

namespace Faunus {

/**
 * @brief Save Matrix Market sparse matrix
 * @note File format description: https://math.nist.gov/MatrixMarket/formats.html
 * @param mat Sparse matrix
 * @param out Output stream
 * @param symmetric Set to true if symmetric matrix (half storage)
 *
 * Mainly copied from Eigen's unsupported `SparseExtra/MarketIO.h`
 *
 * @todo Support symmetric matrices
 * @warning Untested for symmetric matrices with non-zero diagonal
 */
template <typename SparseMatrixType>
bool streamMarket(const SparseMatrixType& mat, std::ostream& out, bool symmetric = false) {
    using RealScalar = typename SparseMatrixType::RealScalar;
    using Scalar = Eigen::SparseMatrix<double>::Scalar;
    if (!out) {
        return false;
    }
    out.flags(std::ios_base::scientific);
    out.precision(std::numeric_limits<RealScalar>::digits10 + 2);
    out << fmt::format("%%MatrixMarket matrix coordinate  {} {}\n", "real", symmetric ? "symmetric" : "general");
    out << mat.rows() << " " << mat.cols() << " " << mat.nonZeros() / (symmetric ? 2 : 1) << "\n";
    for (int col = 0; col < mat.outerSize(); ++col) {
        for (Eigen::SparseMatrix<Scalar>::InnerIterator it(mat, col); it; ++it) {
            if (symmetric && col > it.row()) {
                continue;
            }
            out << it.row() + 1 << " " << it.col() + 1 << " " << it.value() << "\n";
        }
    }
    return true;
}
} // namespace Faunus
#endif // FAUNUS_MATRIXMARKET_H
