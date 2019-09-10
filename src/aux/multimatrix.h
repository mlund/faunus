#pragma once
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
template <int dim, typename Tindices = std::array<int, dim>> struct RowMajorOffset {
    /**
     * @brief memory offset
     * @param N matrix dimentions
     * @param n matrix indices
     * @return memory offset
     */
    inline int operator()(const Tindices &N, const Tindices &n) {
        int offset = 0;
        for (int i = 0; i < dim; i++) {
            int prod = 1;
            for (int j = i + 1; j < dim; j++)
                prod *= N[j];
            offset += prod * n[i];
        }
        return offset;
    }
};

// dynamic 2d matrix with contiguous memory
template <typename T, typename Tindices = std::array<int, 2>> class DynamicArray2D {
  public:
    std::vector<T> data; // contiguous data block; row- layout
    Tindices dim;        // matrix size
    void resize(const Tindices &dimensions) {
        dim = dimensions;
        data.resize(dim[0] * dim[1], T());
    }
    inline size_t index(const Tindices &i) { return i[1] + dim[1] * i[0]; }
    inline T &operator()(const Tindices &i) { return i[1] + dim[1] * i[0]; }
    inline const T &operator()(const Tindices &i) const { return i[1] + dim[1] * i[0]; }
};

// dynamic 3d matrix with contiguous memory
template <typename T, typename Tindices = std::array<int, 3>> class DynamicArray3D {
  public:
    std::vector<T> data; // contiguous data block; row-column-depth layout
    Tindices dim;        // matrix size
    void resize(const Tindices &dimensions) {
        dim = dimensions;
        data.resize(dim[0] * dim[1] * dim[2], T());
    }
    inline size_t index(const Tindices &i) { return i[2] + dim[2] * (i[1] + dim[1] * i[0]); }
    inline T &operator()(const Tindices &i) { return data[i[2] + dim[2] * (i[1] + dim[1] * i[0])]; }
    inline const T &operator()(const Tindices &i) const { return data[i[2] + dim[2] * (i[1] + dim[1] * i[0])]; }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] RowMajor3DMatrix") {
    DynamicArray3D<double, Eigen::Vector3i> m;
    Eigen::Vector3i dim = {4, 10, 2};
    Eigen::Vector3i one = {1, 1, 1};
    m.resize(dim);
    m({2, 3, 1}) = 0.1;
    CHECK(m.data[m.index({2, 3, 1})] == 0.1); // access element via index
    CHECK(&m({0, 0, 0}) == &m.data.front());  // access first element
    CHECK(&m(dim - one) == &m.data.back());   // access last element

    RowMajorOffset<3, Eigen::Vector3i> offset;
    CHECK(m.index(dim - one) == m.data.size() - 1);      // index of last element
    CHECK(m.index(dim - one) == offset(dim, dim - one)); // index of last element
    CHECK(m.index(dim) == offset(dim, dim));             // index beyond last element
    int cnt = 0;
    for (int k = 0; k < dim[0]; k++)
        for (int l = 0; l < dim[1]; l++)
            for (int m = 0; m < dim[2]; m++)
                cnt++;
    CHECK(cnt == m.data.size()); // count number of elements
}
#endif
} // namespace Faunus
