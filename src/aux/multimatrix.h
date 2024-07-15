#pragma once
#include <doctest/doctest.h>
#include <array>
#include <vector>
#include <Eigen/Core>

namespace Faunus {
/**
 * @brief memory offset for arbitrary dimensions, row-major layout
 * @todo this can be used to construct an N-dimensional matrix class
 *
 * https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays
 */
template <int dim, typename Tindices = std::array<int, dim>> struct RowMajorOffset
{
    /**
     * @brief memory offset
     * @param N matrix dimentions
     * @param n matrix indices
     * @return memory offset
     */
    inline int operator()(const Tindices& N, const Tindices& n)
    {
        int offset = 0;
        for (int i = 0; i < dim; i++) {
            int prod = 1;
            for (int j = i + 1; j < dim; j++) {
                prod *= N[j];
            }
            offset += prod * n[i];
        }
        return offset;
    }
};

// dynamic 2d matrix with contiguous memory
template <typename T, typename Tindices = std::array<int, 2>> class DynamicArray2D
{
  public:
    std::vector<T> data; // contiguous data block; row- layout
    Tindices dim;        // matrix size

    void resize(const Tindices& dimensions)
    {
        dim = dimensions;
        data.resize(dim[0] * dim[1], T());
    }

    inline size_t index(const Tindices& i) { return i[1] + dim[1] * i[0]; }

    inline T& operator()(const Tindices& i) { return i[1] + dim[1] * i[0]; }

    inline const T& operator()(const Tindices& i) const { return i[1] + dim[1] * i[0]; }
};

// dynamic 3d matrix with contiguous memory
template <typename T, typename Tindices = std::array<int, 3>> class DynamicArray3D
{
  public:
    std::vector<T> data; // contiguous data block; row-column-depth layout
    Tindices dim;        // matrix size

    void resize(const Tindices& dimensions)
    {
        dim = dimensions;
        data.resize(dim[0] * dim[1] * dim[2], T());
    }

    inline size_t index(const Tindices& i) { return i[2] + dim[2] * (i[1] + dim[1] * i[0]); }

    inline T& operator()(const Tindices& i) { return data[i[2] + dim[2] * (i[1] + dim[1] * i[0])]; }

    inline const T& operator()(const Tindices& i) const
    {
        return data[i[2] + dim[2] * (i[1] + dim[1] * i[0])];
    }
};

TEST_CASE("[Faunus] RowMajor3DMatrix")
{
    DynamicArray3D<double, Eigen::Vector3i> matrix;
    Eigen::Vector3i dim = {4, 10, 2};
    Eigen::Vector3i one = {1, 1, 1};
    matrix.resize(dim);
    matrix({2, 3, 1}) = 0.1;
    CHECK_EQ(matrix.data[matrix.index({2, 3, 1})], 0.1); // access element via index
    CHECK_EQ(&matrix({0, 0, 0}), &matrix.data.front());  // access first element
    CHECK_EQ(&matrix(dim - one), &matrix.data.back());   // access last element

    RowMajorOffset<3, Eigen::Vector3i> offset;
    CHECK_EQ(matrix.index(dim - one), matrix.data.size() - 1); // index of last element
    CHECK_EQ(matrix.index(dim - one), offset(dim, dim - one)); // index of last element
    CHECK_EQ(matrix.index(dim), offset(dim, dim));             // index beyond last element
    int cnt = 0;
    for (int k = 0; k < dim[0]; k++) {
        for (int l = 0; l < dim[1]; l++) {
            for (int m = 0; m < dim[2]; m++) {
                cnt++;
            }
        }
    }
    CHECK_EQ(cnt, matrix.data.size()); // count number of elements
}
} // namespace Faunus
