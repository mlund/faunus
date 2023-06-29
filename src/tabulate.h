#pragma once

#include <iostream>
#include <cstdlib>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>
#include <memory>
#include <concepts>

namespace Faunus {

/**
 * @brief Class templates for tabulation of pair potentials
 *
 * This contains general 1D tabulation templates that can
 * be used to tabulate any function f(x) in a specified
 * interval. This is used to tabulate pair potentials
 *
 * @note Code in this namespace has no dependencies and can be
 *       used for general tabulation outside Faunus.
 */

namespace Tabulate {

/* base class for all tabulators - no dependencies */
template <std::floating_point T = double> class TabulatorBase {
  protected:
    T utol = 1e-5, ftol = -1, umaxtol = -1, fmaxtol = -1;
    T numdr = 0.0001; // dr for derivative evaluation

    // First derivative with respect to x
    T f1(std::function<T(T)> f, T x) const { return (f(x + numdr * 0.5) - f(x - numdr * 0.5)) / (numdr); }

    // Second derivative with respect to x
    T f2(std::function<T(T)> f, T x) const { return (f1(f, x + numdr * 0.5) - f1(f, x - numdr * 0.5)) / (numdr); }

    void check() const {
        if (ftol != -1 && ftol <= 0.0) {
            std::cerr << "ftol=" << ftol << " too small\n" << std::endl;
            abort();
        }
        if (umaxtol != -1 && umaxtol <= 0.0) {
            std::cerr << "umaxtol=" << umaxtol << " too small\n" << std::endl;
            abort();
        }
        if (fmaxtol != -1 && fmaxtol <= 0.0) {
            std::cerr << "fmaxtol=" << fmaxtol << " too small\n" << std::endl;
            abort();
        }
    }

  public:
    struct data {
        std::vector<T> r2;      // r2 for intervals
        std::vector<T> c;       // c for coefficents
        T rmin2 = 0, rmax2 = 0; // useful to save these with table
        bool empty() const { return r2.empty() && c.empty(); }
        inline size_t numKnots() const { return r2.size(); }
    };

    void setTolerance(T _utol, T _ftol = -1, T _umaxtol = -1, T _fmaxtol = -1) {
        utol = _utol;
        ftol = _ftol;
        umaxtol = _umaxtol;
        fmaxtol = _fmaxtol;
    }

    void setNumdr(T _numdr) { numdr = _numdr; }
};

/**
 * @brief Andrea table with logarithmic search
 *
 * Tabulator with logarithmic search.
 * Code mainly from MolSim (Per Linse) with some upgrades
 * Reference: doi:10/frzp4d
 *
 * @note Slow on Intel compiler
 * @todo Hide data and functions; clean up r vs r2 mess.
 */
template <std::floating_point T = double> class Andrea : public TabulatorBase<T> {
  private:
    typedef TabulatorBase<T> base; // for convenience
    int mngrid = 1200;             // Max number of controlpoints
    int ndr = 100;                 // Max number of trials to decr dr
    T drfrac = 0.9;                // Multiplicative factor to decr dr

    std::vector<T> SetUBuffer(T, T zlow, T, T zupp, T u0low, T u1low, T u2low, T u0upp, T u1upp, T u2upp) {

        // Zero potential and force return no coefficients
        if (std::fabs(u0low) < 1e-9)
            if (std::fabs(u1low) < 1e-9)
                return {0, 0, 0, 0, 0, 0, 0};

        T dz1 = zupp - zlow;
        T dz2 = dz1 * dz1;
        T dz3 = dz2 * dz1;
        T w0low = u0low;
        T w1low = u1low;
        T w2low = u2low;
        T w0upp = u0upp;
        T w1upp = u1upp;
        T w2upp = u2upp;
        T c0 = w0low;
        T c1 = w1low;
        T c2 = w2low * 0.5;
        T a = 6 * (w0upp - c0 - c1 * dz1 - c2 * dz2) / dz3;
        T b = 2 * (w1upp - c1 - 2 * c2 * dz1) / dz2;
        T c = (w2upp - 2 * c2) / dz1;
        T c3 = (10 * a - 12 * b + 3 * c) / 6;
        T c4 = (-15 * a + 21 * b - 6 * c) / (6 * dz1);
        T c5 = (2 * a - 3 * b + c) / (2 * dz2);

        return {zlow, c0, c1, c2, c3, c4, c5};
    }

    /**
     * @returns boolean vector.
     * - `[0]==true`: tolerance is approved,
     * - `[1]==true` Repulsive part is found.
     */
    std::vector<bool> CheckUBuffer(std::vector<T>& ubuft, T rlow, T rupp, std::function<T(T)> f) const {

        // Number of points to control
        int ncheck = 11;
        T dr = (rupp - rlow) / (ncheck - 1);
        std::vector<bool> vb(2, false);

        for (int i = 0; i < ncheck; i++) {
            T r1 = rlow + dr * ((T)i);
            T r2 = r1 * r1;
            T u0 = f(r2);
            T u1 = base::f1(f, r2);
            T dz = r2 - rlow * rlow;
            T usum =
                ubuft.at(1) +
                dz * (ubuft.at(2) + dz * (ubuft.at(3) + dz * (ubuft.at(4) + dz * (ubuft.at(5) + dz * ubuft.at(6)))));

            T fsum = ubuft.at(2) +
                     dz * (2 * ubuft.at(3) + dz * (3 * ubuft.at(4) + dz * (4 * ubuft.at(5) + dz * (5 * ubuft.at(6)))));

            if (std::fabs(usum - u0) > base::utol)
                return vb;
            if (base::ftol != -1 && std::fabs(fsum - u1) > base::ftol)
                return vb;
            if (base::umaxtol != -1 && std::fabs(usum) > base::umaxtol)
                vb[1] = true;
            if (base::fmaxtol != -1 && std::fabs(usum) > base::fmaxtol)
                vb[1] = true;
        }
        vb[0] = true;
        return vb;
    }

  public:
    /**
     * @brief Get tabulated value at f(x)
     * @param d Table data
     * @param r2 value
     * @note Auto-vectorization in Clang: https://llvm.org/docs/Vectorizers.html
     */
    inline T eval(const typename base::data& d, T r2) const {
        size_t pos = std::lower_bound(d.r2.begin(), d.r2.end(), r2) - d.r2.begin() - 1;
        size_t pos6 = 6 * pos;
        assert((pos6 + 5) < d.c.size());
        T dz = r2 - d.r2[pos];
        if constexpr (true) { // loop version
            T sum = 0;
#pragma clang loop vectorize(enable) interleave(enable)
            for (size_t i = 5; i > 0; i--)
                sum = dz * (sum + d.c[pos6 + i]);
            return sum + d.c[pos6];
        } else // manually unrolled version
            return d.c[pos6] +
                   dz * (d.c[pos6 + 1] +
                         dz * (d.c[pos6 + 2] + dz * (d.c[pos6 + 3] + dz * (d.c[pos6 + 4] + dz * (d.c[pos6 + 5])))));
    }

    /**
     * @brief Get tabulated value at df(x)/dx
     * @param d Table data
     * @param r2 value
     */
    T evalDer(const typename base::data& d, T r2) const {
        size_t pos = std::lower_bound(d.r2.begin(), d.r2.end(), r2) - d.r2.begin() - 1;
        size_t pos6 = 6 * pos;
        T dz = r2 - d.r2[pos];
        return (d.c[pos6 + 1] +
                dz * (2.0 * d.c[pos6 + 2] +
                      dz * (3.0 * d.c[pos6 + 3] + dz * (4.0 * d.c[pos6 + 4] + dz * (5.0 * d.c[pos6 + 5])))));
    }

    /**
     * @brief Tabulate f(x) in interval ]min,max]
     */
    typename base::data generate(std::function<T(T)> f, double rmin, double rmax) {
        rmin = std::sqrt(rmin);
        rmax = std::sqrt(rmax);
        base::check();
        typename base::data td;
        td.rmin2 = rmin * rmin;
        td.rmax2 = rmax * rmax;

        T rumin = rmin;
        T rmax2 = rmax * rmax;
        T dr = rmax - rmin;
        T rupp = rmax;
        T zupp = rmax2;
        bool repul = false; // Stop tabulation if repul is true

        td.r2.push_back(zupp);

        int i;
        for (i = 0; i < mngrid; i++) {
            T rlow = rupp;
            T zlow;
            std::vector<T> ubuft;
            int j;

            dr = (rupp - rmin);

            for (j = 0; j < ndr; j++) {
                zupp = rupp * rupp;
                rlow = rupp - dr;
                if (rumin > rlow)
                    rlow = rumin;

                zlow = rlow * rlow;

                T u0low = f(zlow);
                T u1low = base::f1(f, zlow);
                T u2low = base::f2(f, zlow);
                T u0upp = f(zupp);
                T u1upp = base::f1(f, zupp);
                T u2upp = base::f2(f, zupp);

                ubuft = SetUBuffer(rlow, zlow, rupp, zupp, u0low, u1low, u2low, u0upp, u1upp, u2upp);
                std::vector<bool> vb = CheckUBuffer(ubuft, rlow, rupp, f);
                repul = vb[1];
                if (vb[0]) {
                    rupp = rlow;
                    break;
                }
                dr *= drfrac;
            }

            if (j >= ndr)
                throw std::runtime_error("Andrea spline: try to increase utol/ftol");
            if (ubuft.size() != 7)
                throw std::runtime_error("Andrea spline: wrong size of ubuft, min value + 6 coefficients");

            td.r2.push_back(zlow);
            for (size_t k = 1; k < ubuft.size(); k++)
                td.c.push_back(ubuft.at(k));

            // Entered a highly repulsive part, stop tabulation
            if (repul) {
                rumin = rlow;
                td.rmin2 = rlow * rlow;
            }
            if (rlow <= rumin || repul)
                break;
        }

        if (i >= mngrid)
            throw std::runtime_error("Andrea spline: try to increase utol/ftol");

        // create final reversed c and r2
        assert(td.c.size() % 6 == 0);
        assert(td.c.size() / (td.r2.size() - 1) == 6);
        assert(std::is_sorted(td.r2.rbegin(), td.r2.rend()));
        std::reverse(td.r2.begin(), td.r2.end());       // reverse all elements
        for (size_t i = 0; i < td.c.size() / 2; i += 6) // reverse knot order in packets of six
            std::swap_ranges(td.c.begin() + i, td.c.begin() + i + 6, td.c.end() - i - 6); // c++17 only
        return td;
    }
};
} // namespace Tabulate
} // namespace Faunus

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] Andrea") {
    using doctest::Approx;
    using namespace Faunus::Tabulate;

    auto f = [](double x) { return 0.5 * x * std::sin(x) + 2; };
    Andrea<double> spline;
    spline.setTolerance(2e-6, 1e-4); // ftol carries no meaning
    auto d = spline.generate(f, 0, 10);

    CHECK(d.r2.size() == 19);
    CHECK(d.c.size() == 108);
    CHECK(d.numKnots() == 19);
    CHECK(d.rmin2 == Approx(0.0));
    CHECK(d.rmax2 == Approx(10.0));
    CHECK(d.r2.at(0) == Approx(0.0));
    CHECK(d.r2.at(1) == Approx(0.212991));
    CHECK(d.r2.at(2) == Approx(0.782554));
    CHECK(d.r2.back() == Approx(10.0));

    CHECK(d.c.at(0) == Approx(2.0));
    CHECK(d.c.at(1) == Approx(0.0));
    CHECK(d.c.at(2) == Approx(0.5));
    CHECK(d.c.back() == Approx(-0.0441931));

    CHECK(spline.eval(d, 1e-9) == Approx(f(1e-9)));
    CHECK(spline.eval(d, 5) == Approx(f(5)));
    CHECK(spline.eval(d, 10) == Approx(f(10)));

    // Check if numerical derivation of *splined* function
    // matches the analytical solution in `evalDer()`.
    auto f_prime = [&](double x, double dx = 1e-10) {
        return (spline.eval(d, x + dx) - spline.eval(d, x - dx)) / (2 * dx);
    };
    double x = 1e-9;
    CHECK(spline.evalDer(d, x) == Approx(f_prime(x)));
    x = 1;
    CHECK(spline.evalDer(d, x) == Approx(f_prime(x)));
    x = 5;
    CHECK(spline.evalDer(d, x) == Approx(f_prime(x)));

    // Check if analytical spline derivative matches
    // derivative of original function
    auto f_prime_exact = [&](double x, double dx = 1e-10) { return (f(x + dx) - f(x - dx)) / (2 * dx); };
    x = 1e-9;
    CHECK(spline.evalDer(d, x) == Approx(f_prime_exact(x)));
    x = 1;
    CHECK(spline.evalDer(d, x) == Approx(f_prime_exact(x)));
    x = 5;
    CHECK(spline.evalDer(d, x) == Approx(f_prime_exact(x)));
}
#endif
