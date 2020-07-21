#pragma once
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
 */
template <class T, bool triangular = false> class PairMatrix {
  private:
    T default_value; // default value when resizing
    std::vector<std::vector<T>> m;

  public:
    void resize(size_t n) {
        m.resize(n);
        for (size_t i = 0; i < m.size(); i++) {
            if constexpr (triangular) {
                m[i].resize(i + 1, default_value);
            } else {
                m[i].resize(n, default_value);
            }
        }
    }

    PairMatrix(size_t n = 0, T val = T()) : default_value(val) { resize(n); }

    auto size() const { return m.size(); }

    inline const T &operator()(size_t i, size_t j) const {
        if constexpr (triangular) {
            if (j > i) {
                std::swap(i, j);
            }
        }
        assert(i < m.size());
        assert(j < m[i].size());
        return m[i][j];
    }

    void set(size_t i, size_t j, T val) {
        if (j > i) {
            std::swap(i, j);
        }
        if (i >= m.size()) {
            resize(i + 1);
        }
        if constexpr (!triangular) {
            m[j][i] = val;
        }
        m[i][j] = val;
    }
};
#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] PairMatrix") {
    int i = 2, j = 3; // particle type, for example

    SUBCASE("full matrix") {
        PairMatrix<double, false> m;
        m.set(i, j, 12.1);
        CHECK(m.size() == 4);
        CHECK(m(i, j) == 12.1);
        CHECK(m(i, j) == m(j, i));
        CHECK(m(0, 2) == 0);
        CHECK(m(2, 0) == 0);
    }

    SUBCASE("full matrix - default value") {
        PairMatrix<double, false> m(5, 3.1);
        for (size_t i = 0; i < 5; i++)
            for (size_t j = 0; j < 5; j++)
                CHECK(m(i, j) == 3.1);
    }

    SUBCASE("triangular matrix - default value") {
        PairMatrix<double, true> m(5, 3.1);
        for (size_t i = 0; i < 5; i++)
            for (size_t j = 0; j < 5; j++)
                CHECK(m(i, j) == 3.1);
    }

    SUBCASE("triangular matrix") {
        PairMatrix<double, true> m;
        m.set(i, j, 12.1);
        CHECK(m.size() == 4);
        CHECK(m(i, j) == 12.1);
        CHECK(m(i, j) == m(j, i));
        CHECK(m(0, 2) == 0);
        CHECK(m(2, 0) == 0);
    }
}
#endif

} // namespace Faunus