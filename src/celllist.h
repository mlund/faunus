#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cmath>
#include <array>
#include <Eigen/Core>

namespace Faunus {

// dynamic 3d matrix with contiguous memory
// https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays
template <typename T, typename Tindices = std::array<int, 3>> class DymamicArray3D {
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
    DymamicArray3D<double, Eigen::Vector3i> m;
    Eigen::Vector3i dim = {4, 10, 2};
    Eigen::Vector3i one = {1, 1, 1};
    m.resize(dim);
    m({2, 3, 1}) = 0.1;
    assert(m.data[m.index({2, 3, 1})] == 0.1);       // access element via index
    assert(&m({0, 0, 0}) == &m.data.front());        // access first element
    assert(&m(dim - one) == &m.data.back());         // access last element
    assert(m.index(dim - one) == m.data.size() - 1); // index of last element
    int cnt = 0;
    for (int k = 0; k < dim[0]; k++)
        for (int l = 0; l < dim[1]; l++)
            for (int m = 0; m < dim[2]; m++)
                cnt++;
    assert(cnt == m.data.size()); // count number of elements
}
#endif

/**
 * @brief Cuboidal cell list with periodic boundaries
 *
 * Maps cartesian points to a grid of arbitrary resolution that
 * stores particle index.
 *
 * - cartesian space is assumed to use all 8 octants (i.e. +/i round 0,0,0)
 * - grid space use only the first octant (all +)
 * - resolution and size is set by `resize`
 * - index of neighbors to a grid point
 *   is obtained with `neighbors()`
 * - index can be moved from one grid point to another with `move()`
 * - the list of particle index in each grid point is stored in a
 *   `std::set<int>` container.
 *
 * @todo
 * - Update cell list based on Change object
 * - Make a non-periodic version
 * - std::set --> set::vector?
 *
 * @date Malmo, March 2018
 */
template <typename CellPoint = Eigen::Vector3i> class CellList {
    typedef size_t Tindex;
    typedef Eigen::Vector3d Point;
    Point halfbox;
    double cellsize = 0;                                           // cell side length (angstrom)
    std::vector<std::vector<std::vector<std::set<Tindex>>>> cells; // dense storage

    std::vector<std::set<Tindex>> cells_dense;    // dense storage (to replace `cells`)
    std::map<int, std::set<Tindex>> cells_sparse; // sparse storage

    inline std::set<Tindex> &row_major(const CellPoint &c) {
        return cells_dense[c[0] * KLM[0] * KLM[2] + c[1] * KLM[2] + c[2]]; // to replace operator[] below
    }                                                                      // row-major ordering of dense storage

  public:
    bool sparse = false;
    CellPoint KLM = {0, 0, 0}; // max cell index K,L,M

    std::set<Tindex> &operator[](const CellPoint &c) {
        return cells[c[0]][c[1]][c[2]];
    } //!< returns set with all index in given cell (complexity: constant)

    CellPoint p2c(const Point &p) const {
        return ((p + halfbox) / cellsize).array().round().template cast<int>();
    } //!< cartesian point --> cell point

    Point c2p(const CellPoint &c) const {
        return (c.template cast<double>()) * cellsize - halfbox;
    } //!< cell point --> cartesian point

    void move(Tindex i, const CellPoint &src, const CellPoint &dst) {
        assert((*this)[src].count(i) == 1 && "i not present in old cell");
        assert((*this)[dst].count(i) == 0 && "i already in new cell");
        (*this)[src].erase(i);
        (*this)[dst].insert(i);
    } //!< move particle index i from one cell to another (complexity: log N)

    void resize(const Point &box, double cutoff) {
        clear();
        halfbox = 0.5 * box;
        cellsize = cutoff;
        KLM = (box / cellsize).array().round().template cast<int>();
        if (KLM.minCoeff() >= 3) {
            cells.resize(KLM[0] + 1);
            for (auto &k : cells) {
                k.resize(KLM[1] + 1);
                for (auto &l : k)
                    l.resize(KLM[2] + 1);
            }
        } else
            throw std::runtime_error("celllist error: too few grid point - cutoff or box too small");
    }

    void clear() {
        for (auto &k : cells)
            for (auto &l : k)
                for (auto &m : l)
                    m.clear();
    } //<! clear all index in cell list

    template <class Tpvec, class T = std::function<Point(const typename Tpvec::value_type &)>>
    void update(
        const Tpvec &p, T getpos = [](auto &i) { return i; }) {
        clear();
        for (Tindex i = 0; i < p.size(); i++)
            (*this)[p2c(getpos(p[i]))].insert(i);
    }

    void neighbors(const Eigen::Vector3i &c, std::vector<Tindex> &index, bool clear = true) const {
        if (clear)
            index.clear();
        int cnt = 0;
        std::array<int, 3> k = {{c[0] - 1, c[0], c[0] + 1}}, l = {{c[1] - 1, c[1], c[1] + 1}},
                           m = {{c[2] - 1, c[2], c[2] + 1}};
        if (k[0] < 0)
            k[0] = KLM[0];
        else if (k[2] > KLM[0])
            k[2] = 0;
        if (l[0] < 0)
            l[0] = KLM[1];
        else if (l[2] > KLM[1])
            l[2] = 0;
        if (m[0] < 0)
            m[0] = KLM[2];
        else if (m[2] > KLM[2])
            m[2] = 0;
        for (int _k : k)
            for (int _l : l)
                for (int _m : m) {
                    cnt++;
                    auto &s = cells[_k][_l][_m];
                    std::copy(s.begin(), s.end(), std::back_inserter(index));
                }
        assert(cnt == 27);
    } //!< Index from all 26+1 neighboring+own cells (complexity: N neighbors)
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] CellList") {
    typedef Eigen::Vector3d Point;
    Point box = {10, 20, 6};
    CellList<Eigen::Vector3i> l;
    l.resize(box, 2);
    CHECK(l.KLM == Eigen::Vector3i(5, 10, 3));
    CHECK(l.p2c({5, 10, 3}) == l.KLM);
    CHECK(l.p2c({-5, -10, -3}) == Eigen::Vector3i(0, 0, 0));
    CHECK(l.p2c({0, 0, 0}) == Eigen::Vector3i(3, 5, 2));

    std::vector<size_t> index; // index of neighbors (and self) in...
    std::vector<Point> vec;    // ...array of points

    vec = {{0, 0, 0}, {0, 3, 0}};
    l.update(vec);
    l.neighbors(l.p2c(vec[0]), index);
    CHECK(index.size() == 1);  // alone by myself...
    CHECK(index.front() == 0); // ...am I really me?

    vec = {{0, 0, 0}, {0, -3, 0}};
    l.update(vec);
    l.neighbors(l.p2c(vec[0]), index);
    CHECK(index.size() == 2); // now we're two
    l.neighbors(l.p2c(vec[1]), index);
    CHECK(index.size() == 2); // now we're two
}
#endif
} // namespace Faunus
