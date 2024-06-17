#pragma once
#include <doctest/doctest.h>
#include <vector>
namespace Faunus {
/**
 * @brief Container for data between pairs
 *
 * Symmetric, dynamic NxN matrix for storing data
 * about pairs. Set values with `set()`. If `triangular==true`
 * the memory usage is reduced but introduces an `if`-statement
 * upon access.
 *
 * ~~~ cpp
 *     int i=2,j=3; // particle type, for example
 *     PairMatrix<double> m;
 *     m.set(i,j,12.0);
 *     cout << m(i,j);         // -> 12.0
 *     cout << m(i,j)==m(j,i); // -> true
 * ~~~
 *
 * @note Vector of vector does not allocate contiguous memory
 */
template <class T, bool triangular = false> class PairMatrix {
  private:
    T default_value; // default value when resizing
    std::vector<std::vector<T>> matrix;

  public:
    void resize(size_t n) {
        matrix.resize(n);
        for (size_t i = 0; i < matrix.size(); i++) {
            if constexpr (triangular) {
                matrix[i].resize(i + 1, default_value);
            } else {
                matrix[i].resize(n, default_value);
            }
        }
    }

    PairMatrix(size_t n = 0, T val = T()) : default_value(val) { resize(n); }

    auto size() const { return matrix.size(); }

    inline const T &operator()(size_t i, size_t j) const {
        if constexpr (triangular) {
            if (j > i) {
                std::swap(i, j);
            }
        }
        assert(i < matrix.size());
        assert(j < matrix[i].size());
        return matrix[i][j];
    }

    void set(size_t i, size_t j, T val) {
        if (j > i) {
            std::swap(i, j);
        }
        if (i >= matrix.size()) {
            resize(i + 1);
        }
        if constexpr (!triangular) {
            matrix[j][i] = val;
        }
        matrix[i][j] = val;
    }

    /** Set a uniform value */
    void set(T val) {
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix.size(); j++) {
                set(i, j, val);
            }
        }
    }

    void setZero() { set(T()); }
};

TEST_CASE("[Faunus] PairMatrix") {
    int i = 2, j = 3; // particle type, for example

    SUBCASE("full matrix") {
        PairMatrix<double, false> m;
        m.set(i, j, 12.1);
        CHECK_EQ(m.size(), 4);
        CHECK_EQ(m(i, j), 12.1);
        CHECK_EQ(m(i, j), m(j, i));
        CHECK_EQ(m(0, 2), 0);
        CHECK_EQ(m(2, 0), 0);
    }

    SUBCASE("full matrix - default value") {
        PairMatrix<double, false> m(5, 3.1);
        for (size_t i = 0; i < 5; i++)
            for (size_t j = 0; j < 5; j++)
                CHECK_EQ(m(i, j), 3.1);
    }

    SUBCASE("triangular matrix - default value") {
        PairMatrix<double, true> m(5, 3.1);
        for (size_t i = 0; i < 5; i++)
            for (size_t j = 0; j < 5; j++)
                CHECK_EQ(m(i, j), 3.1);
    }

    SUBCASE("triangular matrix") {
        PairMatrix<double, true> m;
        m.set(i, j, 12.1);
        CHECK_EQ(m.size(), 4);
        CHECK_EQ(m(i, j), 12.1);
        CHECK_EQ(m(i, j), m(j, i));
        CHECK_EQ(m(0, 2), 0);
        CHECK_EQ(m(2, 0), 0);
    }
}

} // namespace Faunus