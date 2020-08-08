#pragma once
#include <doctest/doctest.h>
#include <vector>
#include <cmath>
#include <string>
#include <limits>

namespace Faunus {
/*
 * @brief Table for storing XY data with random access
 *
 * Functionality:
 *
 * - equidistant x-spacing
 * - supports centered and lower bound (default) binning via `centerbin`
 * - internal storage is `std::vector<Ty>`
 * - lookup complexity: constant
 * - range-based for loops
 * - minimum x value and resolution must be given upon construction
 * - dynamic memory allocation when adding x>=xmin
 * - stream (de)serialization
 * - unit tests
 * - *much* faster, leaner than `std::map`-based Table2D
 */
template <typename Tx = double, typename Ty = Tx, bool centerbin = false> class Equidistant2DTable {
  private:
    Tx _dxinv, _xmin;
    int offset;
    std::vector<Ty> vec;

    inline int to_bin(Tx x) const {
        if constexpr (centerbin) {
            return (x < 0) ? int(x * _dxinv - 0.5) : int(x * _dxinv + 0.5);
        } else {
            return std::floor((x - _xmin) * _dxinv);
        }
    } //!< x-value to vector index

    Tx from_bin(int i) const {
        assert(i >= 0);
        return i / _dxinv + _xmin;
    } //!< vector index to x-value

  public:
    class iterator {
      public:
        typedef Equidistant2DTable<Tx, Ty, centerbin> T;
        iterator(T &tbl, size_t i) : tbl(tbl), i(i) {}
        iterator operator++() {
            ++i;
            return *this;
        }
        bool operator!=(const iterator &other) const { return i != other.i; }
        auto operator*() { return tbl[i]; }

      protected:
        T &tbl;
        size_t i;
    }; // enable range-based for loops

    class const_iterator {
      public:
        typedef Equidistant2DTable<Tx, Ty, centerbin> T;
        const_iterator(const T &tbl, size_t i) : tbl(tbl), i(i) {}
        const_iterator operator++() {
            ++i;
            return *this;
        }
        bool operator!=(const const_iterator &other) const { return i != other.i; }
        auto operator*() const { return tbl[i]; }

      protected:
        const T &tbl;
        size_t i;
    }; // enable range-based for loops

    iterator begin() { return iterator(*this, 0); }
    iterator end() { return iterator(*this, size()); }
    const_iterator begin() const { return const_iterator(*this, 0); }
    const_iterator end() const { return const_iterator(*this, size()); }

    Equidistant2DTable() = default;

    /**
     * @brief Constructor
     * @param dx x spacing
     * @param xmin minimum x value
     * @param xmax maximum x value (for more efficient memory handling, only)
     */
    Equidistant2DTable(Tx dx, Tx xmin, Tx xmax = std::numeric_limits<Tx>::infinity()) { setResolution(dx, xmin, xmax); }

    void setResolution(Tx dx, Tx xmin, Tx xmax = std::numeric_limits<Tx>::infinity()) {
        _xmin = xmin;
        _dxinv = 1 / dx;
        offset = -to_bin(_xmin);
        if (xmax < std::numeric_limits<Tx>::infinity()) {
            vec.reserve(to_bin(xmax) + offset);
        }
    }

    double sumy() const {
        double sum = 0;
        for (auto &i : vec) {
            sum += double(i);
        }
        return sum;
    } //!< sum all y-values

    const std::vector<Ty> &yvec() const { return vec; } //!< vector with y-values

    std::vector<Tx> xvec() const {
        std::vector<Tx> v;
        v.reserve(vec.size());
        for (size_t i = 0; i < vec.size(); i++) {
            v.push_back(from_bin(i));
        }
        return v;
    }

    void clear() { vec.clear(); }
    size_t size() const { return vec.size(); }
    bool empty() const { return vec.empty(); }
    Tx dx() const { return 1 / _dxinv; }
    Tx xmin() const { return _xmin; } //!< minimum x-value

    Tx xmax() const { return (vec.empty()) ? _xmin : from_bin(vec.size() - 1); } //!< maximum stored x-value

    Ty &operator()(Tx x) {
        assert(x >= _xmin);
        int i = to_bin(x) + offset;
        if (i >= vec.size()) {
            vec.resize(i + 1, Ty());
        }
        return vec[i];
    } // return y value for given x

    std::pair<Tx, Ty &> operator[](size_t index) {
        assert(index >= 0);
        assert(index < size());
        return {from_bin(index), vec[index]};
    } // pair with x,y& value

    std::pair<Tx, const Ty &> operator[](size_t index) const {
        assert(size() > index);
        return {from_bin(index), vec[index]};
    } // const pair with x,y& value

    const Ty &operator()(Tx x) const {
        assert(x >= _xmin);
        int i = to_bin(x) + offset;
        assert(i < vec.size());
        return vec.at(i);
    } // return y value for given x

    // can be optinally used to customize streaming out, normalise etc.
    std::function<void(std::ostream &, Tx, Ty)> stream_decorator = nullptr;

    friend std::ostream &operator<<(std::ostream &o, const Equidistant2DTable<Tx, Ty, centerbin> &tbl) {
        o.precision(16);
        for (auto d : tbl) {
            if (tbl.stream_decorator == nullptr) {
                o << d.first << " " << d.second << "\n";
            } else
                tbl.stream_decorator(o, d.first, d.second);
        }
        return o;
    } // write to stream

    auto &operator<<(std::istream &in) {
        assert(stream_decorator == nullptr && "you probably don't want to load a decorated file");
        clear();
        Tx x;
        Ty y;
        std::string line;
        while (std::getline(in, line)) {
            std::stringstream o(line);
            while (o >> x) {
                if ((x >= _xmin) and (y << o)) {
                    operator()(x) = y;
                } else {
                    throw std::runtime_error("table load error: x smaller than xmin");
                }
            }
        }
        return *this;
    } // load from stream
};

TEST_CASE("[Faunus] Equidistant2DTable") {
    using doctest::Approx;

    SUBCASE("centered") {
        Equidistant2DTable<double, double, true> y(0.5, -3.0);
        CHECK(y.xmin() == Approx(-3.0));
        CHECK(y.dx() == Approx(0.5));
        y(-2.51) = 0.15;
        CHECK(y(-2.5) == Approx(0.15));
        y(-2.76) = 0.11;
        CHECK(y(-3) == Approx(0.11));
        y(0.4) = 0.3;
        CHECK(y(0.5) == Approx(0.3));
        y(1.3) = 0.5;
        CHECK(y(1.5) == Approx(0.5));
        CHECK(y.xmax() == Approx(1.5));
    }

    SUBCASE("lower bound") {
        Equidistant2DTable<double> y(0.5, -3.0);
        CHECK(y.xmin() == Approx(-3.0));
        CHECK(y.dx() == Approx(0.5));
        y(-2.51) = 0.15;
        CHECK(y(-3.0) == Approx(0.15));
        y(-2.76) = 0.11;
        CHECK(y(-3) == Approx(0.11));
        y(0.4) = 0.3;
        CHECK(y(0.0) == Approx(0.3));
        y(1.3) = 0.5;
        CHECK(y(1.0) == Approx(0.5));
        CHECK(y.xmax() == Approx(1.0));
    }
}

} // namespace Faunus