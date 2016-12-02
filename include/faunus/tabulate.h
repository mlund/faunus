#ifndef FAUNUS_TABULATE_H
#define FAUNUS_TABULATE_H

#ifndef SWIG

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <functional>
#include <map>
#include <algorithm>
#include <cmath>
#include <memory>

#include <faunus/potentials.h>

#endif

namespace Faunus
{

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

  namespace Tabulate
  {

    /* base class for all tabulators - no dependencies */
    template<typename T=double>
    class TabulatorBase
    {
    protected:
        T utol, ftol, umaxtol, fmaxtol, rmin, rmax;
        T numdr; // dr for derivative evaluation

        // First derivative with respect to x
        T f1( std::function<T( T )> &f, T x ) const
        {
            return (f(x + numdr * 0.5) - f(x - numdr * 0.5)) / (numdr);
        }

        // Second derivative with respect to x
        T f2( std::function<T( T )> &f, T x ) const
        {
            return (f1(f, x + numdr * 0.5) - f1(f, x - numdr * 0.5)) / (numdr);
        }

        void check() const
        {
            assert(rmin >= 0);
            assert(rmax > 0);
            assert(rmin < rmax);
            assert(utol >= 0.0000000000001);
            if ( ftol != -1 && ftol <= 0.0 )
            {
                std::cerr << "ftol=" << ftol << " too small\n" << std::endl;
                abort();
            }
            if ( umaxtol != -1 && umaxtol <= 0.0 )
            {
                std::cerr << "umaxtol=" << umaxtol << " too small\n" << std::endl;
                abort();
            }
            if ( fmaxtol != -1 && fmaxtol <= 0.0 )
            {
                std::cerr << "fmaxtol=" << fmaxtol << " too small\n" << std::endl;
                abort();
            }
        }

    public:
        struct data
        {
            std::vector<T> r2;  // r2 for intervals
            std::vector<T> c;   // c for coefficents
            T rmin2, rmax2;     // useful to save these with table

            bool empty() const
            {
                if ( r2.empty())
                    if ( c.empty())
                        return true;
                return false;
            }

        };

        void setTolerance( T _utol, T _ftol = -1, T _umaxtol = -1, T _fmaxtol = -1 )
        {
            utol = _utol;
            ftol = _ftol;
            umaxtol = _umaxtol;
            fmaxtol = _fmaxtol;
        }

        void setRange( T _rmin, T _rmax )
        {
            rmin = _rmin;
            rmax = _rmax;
        }

        void setNumdr( T _numdr )
        {
            numdr = _numdr;
        }

        TabulatorBase()
        {
            utol = 0.01;
            ftol = -1;
            umaxtol = -1;
            fmaxtol = -1;
            numdr = 0.0001; // dr for derivative evaluation                    
        }

        std::string info( char w = 20 )
        {
            using namespace Faunus::textio;
            std::ostringstream o("");
            o << pad(SUB, w, "Rmin") << rmin << std::endl
              << pad(SUB, w, "Rmax") << rmax << std::endl
              << pad(SUB, w, "Utol") << utol << std::endl
              << pad(SUB, w, "Ftol") << ftol << std::endl
              << pad(SUB, w, "Umaxtol") << umaxtol << std::endl
              << pad(SUB, w, "Fmaxtol") << fmaxtol << std::endl;
            return o.str();
        }
    };

    /**
     * @brief Andrea table with logarithmic search
     *
     * Tabulator with logarithmic search.
     * Code mainly from MolSim (Per Linse) with some upgrades
     * Reference: doi:10/frzp4d
     *
     * @note Slow on Intel compiler
     * @todo Hide data and functions
     */
    template<typename T=double>
    class Andrea : public TabulatorBase<T>
    {
    private:
        typedef TabulatorBase<T> base;// for convenience
        int mngrid; // Max number of controlpoints
        int ndr;    // Max number of trials to decr dr
        T drfrac;   // Multiplicative factor to decr dr

        std::vector<T> SetUBuffer( T rlow,
                                   T zlow,
                                   T rupp,
                                   T zupp,
                                   T u0low,
                                   T u1low,
                                   T u2low,
                                   T u0upp,
                                   T u1upp,
                                   T u2upp )
        {

            std::vector<T> ubuft;
            ubuft.push_back(zlow);

            // Zero potential and force return no coefficients
            if ( fabs(u0low) < 1e-9 && fabs(u1low) < 1e-9 )
            {
                ubuft.resize(7, 0);
                return ubuft;
            }

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

            ubuft.push_back(c0);
            ubuft.push_back(c1);
            ubuft.push_back(c2);
            ubuft.push_back(c3);
            ubuft.push_back(c4);
            ubuft.push_back(c5);

            return ubuft;
        }

        /**
         * @brief Description...
         * @param ubuft description...
         * @returns boolean vector.
         * - `[0]==true`: tolerance is approved,
         * - `[1]==true` Repulsive part is found.
         */
        std::vector<bool> CheckUBuffer( std::vector<T> &ubuft,
                                        T rlow,
                                        T rupp,
                                        std::function<T( T )> f ) const
        {

            // Number of points to control
            int ncheck = 11;
            T dr = (rupp - rlow) / (ncheck - 1);
            std::vector<bool> vb(2, false);

            for ( int i = 0; i < ncheck; i++ )
            {
                T r1 = rlow + dr * ((T) i);
                T r2 = r1 * r1;
                T u0 = f(r2);
                T u1 = base::f1(f, r2);
                T dz = r2 - rlow * rlow;
                T usum = ubuft.at(1) +
                    dz * (ubuft.at(2) +
                        dz * (ubuft.at(3) +
                            dz * (ubuft.at(4) +
                                dz * (ubuft.at(5) +
                                    dz * ubuft.at(6)))));

                T fsum = ubuft.at(2) +
                    dz * (2 * ubuft.at(3) +
                        dz * (3 * ubuft.at(4) +
                            dz * (4 * ubuft.at(5) +
                                dz * (5 * ubuft.at(6)))));

                if ( std::abs(usum - u0) > base::utol )
                    return vb;
                if ( base::ftol != -1 && std::abs(fsum - u1) > base::ftol )
                    return vb;
                if ( base::umaxtol != -1 && std::abs(usum) > base::umaxtol )
                    vb[1] = true;
                if ( base::fmaxtol != -1 && std::abs(usum) > base::fmaxtol )
                    vb[1] = true;
            }
            vb[0] = true;
            return vb;
        }

    public:
        Andrea() : base()
        {
            mngrid = 1200;
            ndr = 100;
            drfrac = 0.9;
        }

        /**
         * @brief Get tabulated value at f(x)
         * @param d Table data
         * @param r2 x value
         */
        T eval( const typename base::data &d, T r2 ) const
        {
            auto low = std::lower_bound(d.r2.begin(), d.r2.end(), r2);
            size_t pos = (low - d.r2.begin() - 1);
            T min = d.r2[pos];
            T dz = r2 - min;
            int pos6 = 6 * pos;
            T usum = d.c[pos6 + 0] +
                dz * (d.c[pos6 + 1] +
                    dz * (d.c[pos6 + 2] +
                        dz * (d.c[pos6 + 3] +
                            dz * (d.c[pos6 + 4] +
                                dz * (d.c[pos6 + 5])))));
            return usum;
        }

        /**
         * @brief Tabulate f(x)
         */
        typename base::data generate( std::function<T( T )> f )
        {
            base::check();
            typename base::data td;
            td.rmax2 = base::rmax * base::rmax;
            td.rmin2 = base::rmin * base::rmin;

            T minv = base::rmin;
            T maxv = base::rmax;
            T rumin = minv;
            T maxv2 = maxv * maxv;
            T dr = maxv - minv;
            T rupp = maxv;
            T zupp = maxv2;
            bool repul = false; // Stop tabulation if repul is true

            td.r2.push_back(zupp);

            int i;
            for ( i = 0; i < mngrid; i++ )
            {
                T rlow = rupp;
                T zlow;
                std::vector<T> ubuft;
                int j;

                dr = (rupp - minv);

                for ( j = 0; j < ndr; j++ )
                {
                    zupp = rupp * rupp;
                    rlow = rupp - dr;
                    if ( rumin > rlow )
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
                    if ( vb[0] == true )
                    {
                        rupp = rlow;
                        break;
                    }
                    dr *= drfrac;
                }
                assert(j < ndr && "Try to increase utol/ftol");
                assert(ubuft.size() == 7 && "Wrong size of ubuft, minvalue + 6 coefficients");
                td.r2.push_back(zlow);
                for ( size_t k = 1; k < ubuft.size(); k++ )
                    td.c.push_back(ubuft.at(k));

                // Entered a highly repulsive part, stop tabulation
                if ( repul == true )
                {
                    rumin = rlow;
                    td.rmin2 = rlow * rlow;
                }
                if ( rlow <= rumin || repul == true )
                    break;
            }
            assert(i < mngrid && "Try to increase utol/ftol");

            // Sort td
            typename base::data tdsort;
            tdsort.rmax2 = td.rmax2;
            tdsort.rmin2 = td.rmin2;
            tdsort.r2.push_back(td.r2.at(td.r2.size() - 1));
            for ( int i = td.r2.size() - 2; i >= 0; i-- )
            {
                tdsort.r2.push_back(td.r2.at(i));
                for ( int j = 0; j < 6; j++ )
                    tdsort.c.push_back(td.c.at(6 * i + j));
            }
            return tdsort;
        }

        /**
         * @brief ...
         */
        typename base::data generate_full( std::function<T( T )> f )
        {
            typename base::data tg = generate(f);

            // Assume zero at from max to "infinity"
            tg.rmax2 = 1e9;
            tg.r2.push_back(pc::infty);
            tg.c.resize(tg.c.size() + 6, 0);

            // Assume infinity from min to zero
            tg.rmin2 = 0;
            tg.r2.insert(tg.r2.begin(), 0);
            vector<T> c(6, 0);
            c[0] = 100000;
            tg.c.insert(tg.c.begin(), c.begin(), c.end());

            return tg;
        }

        /**
         * @brief ...
         */
        typename base::data generate_empty()
        {
            typename base::data tg;
            // Assume zero at from zero to "infinity"
            tg.rmin2 = 0.0;
            tg.rmax2 = 1e10;
            tg.r2.push_back(0);
            tg.r2.push_back(1e10);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            return tg;
        }

        std::string print( typename base::data &d )
        {
            std::ostringstream o;
            o << "Size of r2: " << d.r2.size() << endl
              << "rmax2 r2=" << d.rmax2 << " r=" << sqrt(d.rmax2) << endl
              << "rmin2 r2=" << d.rmin2 << " r=" << sqrt(d.rmin2) << endl;
            for ( size_t i = 0; i < d.r2.size(); i++ )
            {
                o << i << ": r2=" << d.r2.at(i)
                  << " r=" << std::sqrt(d.r2.at(i)) << endl;
                if ( i != (d.r2.size() - 1))
                {
                    o << "coeffs:";
                    for ( int j = 0; j < 6; j++ )
                        o << " " << d.c.at(i * 6 + j) << ",";
                }
                o << endl;
            }
            return o.str();
        }
    };

    /**
     * @brief Andrea table optimized for intel compiler
     *
     * Tabulator with linear search starting from large distances
     * Code mainly from MolSim (Per Linse) with some upgrades
     * Reference: doi:10/frzp4d
     *
     * @todo Inherit from "Andrea"
     */
    template<typename T=double>
    class AndreaIntel : public TabulatorBase<T>
    {
    private:
        typedef TabulatorBase<T> base; // for convenience, only
        int mngrid;// Max number of controlpoints
        int ndr;   // Max number of trials to decr dr
        T drfrac;  // Multiplicative factor to decr dr

        std::vector<T> SetUBuffer( T rlow,
                                   T zlow,
                                   T rupp,
                                   T zupp,
                                   T u0low,
                                   T u1low,
                                   T u2low,
                                   T u0upp,
                                   T u1upp,
                                   T u2upp )
        {

            std::vector<T> ubuft;

            ubuft.push_back(zlow);
            // Zero potential and force return no coefficients
            if ( u0low == 0.0 && u1low == 0.0 )
            {
                for ( int i = 0; i < 6; i++ )
                    ubuft.push_back(0.0);
                return ubuft;
            }

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

            ubuft.push_back(c0);
            ubuft.push_back(c1);
            ubuft.push_back(c2);
            ubuft.push_back(c3);
            ubuft.push_back(c4);
            ubuft.push_back(c5);

            return ubuft;
        }

        std::vector<bool> CheckUBuffer( std::vector<T> &ubuft,
                                        T rlow,
                                        T rupp,
                                        std::function<T( T )> f )
        {
            // vb[0]: Tolerance is approved
            // vb[1]: A repulsive part is found
            std::vector<bool> vb;
            vb.push_back(false);
            vb.push_back(false);

            // Number of points to control
            int ncheck = 11;
            T dr = (rupp - rlow) / (ncheck - 1);

            for ( int i = 0; i < ncheck; i++ )
            {
                T r1 = rlow + dr * ((T) i);
                T r2 = r1 * r1;
                T u0 = f(r2);
                T u1 = base::f1(f, r2);
                T dz = r2 - rlow * rlow;
                T usum = ubuft.at(1) +
                    dz * (ubuft.at(2) +
                        dz * (ubuft.at(3) +
                            dz * (ubuft.at(4) +
                                dz * (ubuft.at(5) +
                                    dz * ubuft.at(6)))));

                T fsum = ubuft.at(2) +
                    dz * (2 * ubuft.at(3) +
                        dz * (3 * ubuft.at(4) +
                            dz * (4 * ubuft.at(5) +
                                dz * (5 * ubuft.at(6)))));

                if ( std::abs(usum - u0) > base::utol )
                    return vb;
                if ( base::ftol != -1 && std::abs(fsum - u1) > base::ftol )
                    return vb;
                if ( base::umaxtol != -1 && std::abs(usum) > base::umaxtol )
                    vb.at(1) = true;
                if ( base::fmaxtol != -1 && std::abs(usum) > base::fmaxtol )
                    vb.at(1) = true;
            }
            vb.at(0) = true;
            return vb;
        }

    public:
        AndreaIntel() : base()
        {
            mngrid = 1200;
            ndr = 100;
            drfrac = 0.9;
        }

        // gets tabulated value at r2
        T eval( const typename base::data &d, T r2 ) const
        {

            // Linear search starting from large distances
            unsigned int i;
            for ( i = 0; i < d.r2.size(); i++ )
                if ( r2 > d.r2[i] )
                    break;
            int pos = i;

            T min = d.r2[pos];
            T dz = r2 - min;
            int pos6 = 6 * (pos - 1);

            T usum =
                d.c[pos6 + 0] +
                    dz * d.c[pos6 + 1] +
                    dz * dz * d.c[pos6 + 2] +
                    dz * dz * dz * d.c[pos6 + 3] +
                    dz * dz * dz * dz * d.c[pos6 + 4] +
                    dz * dz * dz * dz * dz * d.c[pos6 + 5];
            return usum;
        }

        typename base::data generate( std::function<T( T )> f )
        {
            base::check();

            typename base::data td;
            td.rmax2 = base::rmax * base::rmax;
            td.rmin2 = base::rmin * base::rmin;

            T minv = base::rmin;
            T maxv = base::rmax;
            T rumin = minv;
            T maxv2 = maxv * maxv;
            T dr = maxv - minv;
            T rupp = maxv;
            T zupp = maxv2;
            bool repul = false; // Stop tabulation if repul is true
            td.r2.push_back(zupp);

            int i;
            for ( i = 0; i < mngrid; i++ )
            {
                T rlow = rupp;
                T zlow;
                std::vector<T> ubuft;
                dr = (rupp - minv);
                int j;
                for ( j = 0; j < ndr; j++ )
                {
                    zupp = rupp * rupp;
                    rlow = rupp - dr;
                    if ( rumin > rlow )
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
                    repul = vb.at(1);
                    if ( vb.at(0) == true )
                    {
                        rupp = rlow;
                        break;
                    }
                    else
                    {
                    }
                    dr *= drfrac;
                }
                assert(j < ndr && "Try to increase utol/ftol");
                td.r2.push_back(zlow);
                assert(ubuft.size() == 7 && "Wrong size of ubuft, minvalue + 6 coefficients");
                for ( unsigned int k = 1; k < ubuft.size(); k++ )
                    td.c.push_back(ubuft.at(k));

                // Entered a highly repulsive part, stop tabulation
                if ( repul == true )
                {
                    rumin = rlow;
                    td.rmin2 = rlow * rlow;
                }
                if ( rlow <= rumin || repul == true )
                    break;
            }
            assert(i < mngrid && "Increase utol or ftol");
            return td;
        }

        typename base::data generate_full( std::function<T( T )> f )
        {
            typename base::data tg = generate(f);

            // Assume infinity from min to zero
            tg.rmin2 = 0;
            tg.r2.push_back(0);
            tg.c.push_back(100000.0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);

            // Assume zero at from max to "infinity"
            auto it = tg.r2.begin();
            tg.r2.insert(it, pc::infty);
            tg.rmax2 = 10000000000000.0;
            it = tg.c.begin();
            tg.c.insert(it, 0);
            it = tg.c.begin();
            tg.c.insert(it, 0);
            it = tg.c.begin();
            tg.c.insert(it, 0);
            it = tg.c.begin();
            tg.c.insert(it, 0);
            it = tg.c.begin();
            tg.c.insert(it, 0);
            it = tg.c.begin();
            tg.c.insert(it, 0);
            return tg;
        }

        typename base::data generate_empty()
        {
            typename base::data tg;
            // Assume zero at from zero to "infinity"
            tg.rmin2 = 0.0;
            tg.rmax2 = 10000000000000.0;
            tg.r2.push_back(10000000000000.0);
            tg.r2.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            return tg;
        }

        std::string print( typename base::data &d )
        {
            std::ostringstream o("Size of r2: " + d.r2.size());
            o << "rmax2 r2=" << d.rmax2 << " r=" << std::sqrt(d.rmax2) << std::endl;
            o << "rmin2 r2=" << d.rmin2 << " r=" << std::sqrt(d.rmin2) << std::endl;
            for ( unsigned int i = 0; i < d.r2.size(); i++ )
            {
                o << i << ": r2=" << d.r2.at(i) << " r=" << std::sqrt(d.r2.at(i)) << std::endl;
                if ( i != (d.r2.size() - 1))
                {
                    o << "coeffs:";
                    for ( unsigned int j = 0; j < 6; j++ )
                        o << " " << d.c.at(i * 6 + j) << ",";
                }
                o << std::endl;
            }
            return o.str();
        }
    };

    /**
     * @brief Hermite cubic spline
     * @todo Hide data and functions
     */
    template<typename T=double>
    struct Hermite : public TabulatorBase<T>
    {
        typedef TabulatorBase<T> base; // for convenience, only

        int mngrid;    // Max number of controlpoints
        int ndr;        // Max number of trials to decr dr
        T drfrac;  // Multiplicative factor to decr dr

        Hermite() : base()
        {
            mngrid = 1200;
            ndr = 100;
            drfrac = 0.9;
        }

        // gets tabulated value at r2
        T eval( const typename base::data &d, T r2 ) const
        {
            auto low = std::lower_bound(d.r2.begin(), d.r2.end(), r2);
            int pos = (low - d.r2.begin() - 1);
            T min = d.r2[pos];
            T dz = r2 - min;
            int pos4 = 4 * pos;
            T usum = d.c[pos4 + 0] +
                dz * (d.c[pos4 + 1] +
                    dz * (d.c[pos4 + 2] +
                        dz * (d.c[pos4 + 3])
                    ));
            return usum;
        }

        std::vector<T> SetUBuffer( T rlow,
                                   T zlow,
                                   T rupp,
                                   T zupp,
                                   T u0low,
                                   T u1low,
                                   T u0upp,
                                   T u1upp )
        {

            std::vector<T> ubuft(1, zlow);

            // Zero potential and force return no coefficients
            if ( u0low == 0.0 && u1low == 0.0 )
            {
                ubuft.resize(5, 0);
                return ubuft;
            }

            T dz1 = zupp - zlow;
            T dz2 = dz1 * dz1;
            T dz3 = dz2 * dz1;
            T w0low = u0low;
            T w1low = u1low;
            T w0upp = u0upp;
            T w1upp = u1upp;
            T c0 = w0low;
            T c1 = w1low;
            T c2 = 3 * (w0upp - w0low) / dz2 - (w1upp + 2 * w1low) / dz1;
            T c3 = (w1upp + w1low) / dz2 - 2 * (w0upp - w0low) / dz3;

            ubuft.push_back(c0);
            ubuft.push_back(c1);
            ubuft.push_back(c2);
            ubuft.push_back(c3);
            return ubuft;
        }

        std::vector<bool> CheckUBuffer( std::vector<T> &ubuft,
                                        T rlow,
                                        T rupp,
                                        std::function<T( T )> f )
        {
            std::vector<bool> vb;
            vb.push_back(false); // vb[0]: Tolerance is approved
            vb.push_back(false); // vb[1]: A repulsive part is found

            // Number of points to control
            int ncheck = 101;
            T dr = (rupp - rlow) / (ncheck - 1);

            for ( int i = 0; i < ncheck; i++ )
            {
                T r1 = rlow + dr * ((T) i);
                T r2 = r1 * r1;
                T u0 = f(r2);
                T u1 = base::f1(f, r2);
                T dz = r2 - rlow * rlow;
                T usum = ubuft.at(1) +
                    dz * (ubuft.at(2) +
                        dz * (ubuft.at(3) +
                            dz * (ubuft.at(4))
                        ));
                T fsum = ubuft.at(2) +
                    dz * (2.0 * ubuft.at(3) +
                        dz * (3.0 * ubuft.at(4))
                    );

                if ( std::abs(usum - u0) > base::utol )
                    return vb;
                if ( base::ftol != -1 && std::abs(fsum - u1) > base::ftol )
                    return vb;
                if ( base::umaxtol != -1 && std::abs(usum) > base::umaxtol )
                    vb.at(1) = true;
                if ( base::fmaxtol != -1 && std::abs(usum) > base::fmaxtol )
                    vb.at(1) = true;
            }
            vb.at(0) = true;
            return vb;
        }

        typename base::data generate( std::function<T( T )> f )
        {
            base::check();

            typename base::data td;
            td.rmax2 = base::rmax * base::rmax;
            td.rmin2 = base::rmin * base::rmin;

            T minv = base::rmin;
            T maxv = base::rmax;
            T rumin = minv;
            T maxv2 = maxv * maxv;
            T dr = maxv - minv;
            T rupp = maxv;
            T zupp = maxv2;
            bool repul = false; // Stop tabulation if repul is true
            td.r2.push_back(zupp);

            int i;
            for ( i = 0; i < mngrid; i++ )
            {
                T rlow = rupp;
                T zlow;
                std::vector<T> ubuft;
                int j;

                dr = (rupp - minv);

                for ( j = 0; j < ndr; j++ )
                {
                    zupp = rupp * rupp;
                    rlow = rupp - dr;
                    if ( rumin > rlow )
                        rlow = rumin;
                    zlow = rlow * rlow;
                    T u0low = f(zlow);
                    T u1low = base::f1(f, zlow);
                    T u0upp = f(zupp);
                    T u1upp = base::f1(f, zupp);

                    ubuft = SetUBuffer(rlow, zlow, rupp, zupp, u0low, u1low, u0upp, u1upp);
                    std::vector<bool> vb = CheckUBuffer(ubuft, rlow, rupp, f);
                    repul = vb.at(1);
                    if ( vb.at(0) == true )
                    {
                        rupp = rlow;
                        break;
                    }
                    else
                    {
                    }
                    dr *= drfrac;
                }
                assert(j < ndr && "Try to increase utol");
                td.r2.push_back(zlow);
                assert(ubuft.size() == 5 && "Wrong size of ubuft, minvalue + 4 coefficients");
                for ( unsigned int k = 1; k < ubuft.size(); k++ )
                    td.c.push_back(ubuft.at(k));

                // Entered a highly repulsive part, stop tabulation
                if ( repul == true )
                {
                    rumin = rlow;
                    td.rmin2 = rlow * rlow;
                }
                if ( rlow <= rumin || repul == true )
                    break;
            }
            assert(i < mngrid && "Try to increase utol/ftol");

            // Sort td
            typename base::data tdsort;
            tdsort.rmax2 = td.rmax2;
            tdsort.rmin2 = td.rmin2;
            tdsort.r2.push_back(td.r2.at(td.r2.size() - 1));
            for ( int i = td.r2.size() - 2; i >= 0; i-- )
            {
                tdsort.r2.push_back(td.r2.at(i));
                for ( int j = 0; j < 4; j++ )
                    tdsort.c.push_back(td.c.at(4 * i + j));
            }
            return tdsort;
        }

        typename base::data generate_full( std::function<T( T )> f )
        {

            typename base::data tg = generate(f);

            // Assume zero at from max to "infinity"
            tg.rmax2 = 1e9;
            tg.r2.push_back(pc::infty);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);

            // Assume infinity from min to zero
            auto it = tg.r2.begin();
            tg.rmin2 = 0;
            tg.r2.insert(it, 0);

            it = tg.c.begin();
            tg.c.insert(it, 0.0);
            it = tg.c.begin();
            tg.c.insert(it, 0.0);
            it = tg.c.begin();
            tg.c.insert(it, 0.0);
            it = tg.c.begin();
            tg.c.insert(it, 100000.0);
            return tg;
        }

        typename base::data generate_empty()
        {
            typename base::data tg;

            // Assume zero at from zero to "infinity"
            tg.rmin2 = 0.0;
            tg.rmax2 = 1e9;
            tg.r2.push_back(0);
            tg.r2.push_back(1e9);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            tg.c.push_back(0);
            return tg;
        }

        std::string print( typename base::data &d )
        {
            std::ostringstream o("Size of r2: " + d.r2.size());
            o << "rmax2 r2=" << d.rmax2 << " r=" << std::sqrt(d.rmax2) << std::endl;
            o << "rmin2 r2=" << d.rmin2 << " r=" << std::sqrt(d.rmin2) << std::endl;
            for ( unsigned int i = 0; i < d.r2.size(); i++ )
            {
                o << i << ": r2=" << d.r2.at(i) << " r="
                  << std::sqrt(d.r2.at(i)) << std::endl;
                if ( i != (d.r2.size() - 1))
                {
                    o << "coeffs:";
                    for ( unsigned int j = 0; j < 6; j++ )
                        o << " " << d.c.at(i * 6 + j) << ",";
                }
                o << std::endl;
            }
            return o.str();
        }
    };

    /**
     * @brief Linear interpolation spline
     * @todo Hide data and functions
     */
    template<typename T=double>
    struct Linear : public TabulatorBase<T>
    {
        typedef TabulatorBase<T> base; // for convenience, only

        int mngrid;    // Max number of controlpoints
        int ndr;        // Max number of trials to decr dr
        T drfrac;  // Multiplicative factor to decr dr

        Linear() : base()
        {
            mngrid = 1200;
            ndr = 100;
            drfrac = 0.9;
        }

        // gets tabulated value at r2
        T eval( const typename base::data &d, T r2 ) const
        {
            auto low = std::lower_bound(d.r2.begin(), d.r2.end(), r2);
            int pos = (low - d.r2.begin() - 1);
            T min = d.r2[pos];
            T dz = r2 - min;
            int pos2 = 2 * pos;
            T usum = d.c[pos2 + 0] +
                dz * d.c[pos2 + 1];
            return usum;
        }

        std::vector<T> SetUBuffer( T zlow,
                                   T zupp,
                                   T u0low,
                                   T u0upp )
        {

            std::vector<T> ubuft;
            ubuft.push_back(zlow);

            T dz = zupp - zlow;
            T w0low = u0low;
            T w0upp = u0upp;
            T c0 = w0low;
            T c1 = (w0upp - w0low) / dz;

            ubuft.push_back(c0);
            ubuft.push_back(c1);

            return ubuft;
        }

        std::vector<bool> CheckUBuffer( std::vector<T> &ubuft,
                                        T rlow,
                                        T rupp,
                                        std::function<T( T )> f )
        {

            std::vector<bool> vb(2);
            vb[0] = false; // vb[0]: Tolerance is approved
            vb[1] = false; // vb[1]: A repulsive part is found

            // Number of points to control
            int ncheck = 11;
            T dr = (rupp - rlow) / (ncheck - 1);

            for ( int i = 0; i < ncheck; i++ )
            {
                T r1 = rlow + dr * ((T) i);
                T r2 = r1 * r1;
                T u0 = f(r2);
                T dz = r2 - rlow * rlow;
                T usum = ubuft.at(1) +
                    dz * ubuft.at(2);
                if ( std::abs(usum - u0) > base::utol )
                    return vb;
                if ( base::umaxtol != -1 && std::abs(usum) > base::umaxtol )
                    vb.at(1) = true;
            }
            vb.at(0) = true;
            return vb;
        }

        typename base::data generate( std::function<T( T )> f )
        {
            base::check();
            typename base::data td;
            td.rmax2 = base::rmax * base::rmax;
            td.rmin2 = base::rmin * base::rmin;
            T minv = base::rmin;
            T maxv = base::rmax;
            T rumin = minv;
            T maxv2 = maxv * maxv;
            T dr = maxv - minv;
            T rupp = maxv;
            T zupp = maxv2;
            bool repul = false; // Stop tabulation if repul is true
            td.r2.push_back(zupp);

            int i;
            for ( i = 0; i < mngrid; i++ )
            {
                T rlow = rupp;
                T zlow;
                std::vector<T> ubuft;
                int j;

                dr = (rupp - minv);

                for ( j = 0; j < ndr; j++ )
                {
                    zupp = rupp * rupp;
                    rlow = rupp - dr;
                    if ( rumin > rlow )
                        rlow = rumin;
                    zlow = rlow * rlow;
                    T u0low = f(zlow);
                    T u0upp = f(zupp);
                    ubuft = SetUBuffer(zlow, zupp, u0low, u0upp);
                    std::vector<bool> vb = CheckUBuffer(ubuft, rlow, rupp, f);
                    repul = vb.at(1);
                    if ( vb.at(0) == true )
                    {
                        rupp = rlow;
                        break;
                    }
                    else
                    {
                    }
                    dr *= drfrac;
                }
                assert(j < ndr && "Try to increase utol/ftol");
                td.r2.push_back(zlow);
                assert(ubuft.size() == 3 && "Wrong size of ubuft, minvalue + 2 coefficients");
                for ( unsigned int k = 1; k < ubuft.size(); k++ )
                    td.c.push_back(ubuft.at(k));

                // Entered a highly repulsive part, stop tabulation
                if ( repul == true )
                {
                    rumin = rlow;
                    td.rmin2 = rlow * rlow;
                }
                if ( rlow <= rumin || repul == true )
                    break;
            }
            // mngrid not enough (increase utol or ftol)
            assert(i < mngrid && "Try to increase utol or ftol");

            // Sort td
            typename base::data tdsort;
            tdsort.rmax2 = td.rmax2;
            tdsort.rmin2 = td.rmin2;
            tdsort.r2.push_back(td.r2.at(td.r2.size() - 1));
            for ( int i = td.r2.size() - 2; i >= 0; i-- )
            {
                tdsort.r2.push_back(td.r2.at(i));
                for ( int j = 0; j < 2; j++ )
                {
                    tdsort.c.push_back(td.c.at(2 * i + j));
                }
            }
            return tdsort;
        }

        typename base::data generate_full( std::function<T( T )> f )
        {
            typename base::data tg = generate(f);
            // Assume zero at from max to "infinity"
            tg.rmax2 = 10000000000000.0;
            tg.r2.push_back(pc::infty);
            tg.c.push_back(0.0);
            tg.c.push_back(0.0);
            // Assume infinity from min to zero
            auto it = tg.r2.begin();
            tg.rmin2 = 0.0;
            tg.r2.insert(it, 0.0);

            it = tg.c.begin();
            tg.c.insert(it, 0.0);
            it = tg.c.begin();
            tg.c.insert(it, 100000.0);

            return tg;
        }

        typename base::data generate_empty()
        {
            typename base::data tg;
            // Assume zero at from zero to "infinity"
            tg.rmin2 = 0.0;
            tg.rmax2 = 10000000000000.0;
            tg.r2.push_back(0.0);
            tg.r2.push_back(10000000000000.0);
            tg.c.push_back(0.0);
            tg.c.push_back(0.0);
            return tg;
        }

        std::string print( typename base::data &d )
        {
            std::ostringstream o("Size of r2: " + d.r2.size());
            o << "rmax2 r2=" << d.rmax2 << " r=" << std::sqrt(d.rmax2) << std::endl
              << "rmin2 r2=" << d.rmin2 << " r=" << std::sqrt(d.rmin2) << std::endl;
            for ( unsigned int i = 0; i < d.r2.size(); i++ )
            {
                o << i << ": r2=" << d.r2.at(i) << " r="
                  << std::sqrt(d.r2.at(i)) << std::endl;
                if ( i != (d.r2.size() - 1))
                {
                    o << "coeffs:";
                    for ( unsigned int j = 0; j < 2; j++ )
                        o << " " << d.c.at(i * 2 + j) << ",";
                }
                o << std::endl;
            }
            return o.str();
        }
    };

  } //Tabulate namespace

#ifdef FAUNUS_POTENTIAL_H

  namespace Potential
  {

    /**
     * @brief Tabulated potential between all particle types
     */
    template<typename Tpairpot, typename Ttabulator=Tabulate::Andrea<double> >
    class PotentialTabulate : public Tpairpot
    {
    private:
        Ttabulator tab;
        typedef opair<int> Tpair;
        std::map<Tpair, typename Ttabulator::data> m;

    public:
        PotentialTabulate( Tmjson &j ) : Tpairpot(j)
        {
            string sec = Tpairpot::jsonsec;
            tab.setRange(
                j[sec]["tab_rmin"] | 1.0,
                j[sec]["tab_rmax"] | 100.0);
            tab.setTolerance(
                j[sec]["tab_utol"] | 0.01,
                j[sec]["tab_ftol"] | -1.0,
                j[sec]["tab_umaxtol"] | -1.0,
                j[sec]["tab_fmaxtol"] | -1.0);
        }

        template<class Tparticle>
        double operator()( const Tparticle &a, const Tparticle &b, double r2 )
        {
            Tpair ab(a.id, b.id);
            auto it = m.find(ab);
            if ( it != m.end())
                return tab.eval(it->second, r2);
            std::function<double( double )> f = [=]( double r2 ) { return Tpairpot(*this)(a, b, r2); };
            m[ab] = tab.generate_full(f);
            return (*this)(a, b, r2);
        }
    };

    /**
     * @brief Construct a distance dependent potential from json entry
     *
     * This will scan the given json entry for pair potentials and,
     * for the given particle ids, return a function of the potential (in kT)
     * that takes only the squared distance, r2, as input.
     * This is useful for splining.
     * If more than one potential are found, these will be added together.
     */
    template<typename Tparticle=PointParticle, typename Tfunc=std::function<double( double )> >
    std::pair<Tfunc, string> mixPairPotential( Tmjson &j, int id1, int id2 )
    {

        // auxiliary lambda function to sum two potential functions
        auto sum = []( const Tfunc &a, const Tfunc &b )
        {
            return [=]( double r2 ) { return a(r2) + b(r2); };
        };

        string brief; // sum of brief() information strings

        // starting potential function, returning zero
        Tfunc f = [=]( double r2 ) { return 0; };

        // the two interacting particles
        Tparticle a, b;
        a = atom[id1];
        b = atom[id2];

        // search for pair potentials in json entry and
        // return them as a function of separation, only
        if ( j.find("coulomb") != j.end())
        {
            Coulomb pot(j);
            brief += " + " + pot.brief();
            f = sum(f, [=]( double r2 ) { return pot(a, b, r2); });
        }
        if ( j.find("debyehuckel") != j.end())
        {
            DebyeHuckel pot(j);
            brief += " + " + pot.brief();
            f = sum(f, [=]( double r2 ) { return pot(a, b, r2); });
        }
        if ( j.find("harmonic") != j.end())
        {
            Harmonic pot(j);
            brief += " + " + pot.brief();
            f = sum(f, [=]( double r2 ) { return pot(a, b, r2); });
        }
        if ( j.find("lennardjones") != j.end())
        {
            LennardJonesLB pot(j);
            brief += " + " + pot.brief();
            f = sum(f, [=]( double r2 ) { return pot(a, b, r2); });
        }
        if ( j.find("fromdisk") != j.end())
        {
            FromDisk<> pot(j);
            brief += " + " + pot.brief();
            f = sum(f, [=]( double r2 ) { return pot(a, b, r2); });
        }


        // delete initial " + "
        if ( !brief.empty())
            brief = brief.substr(3);

        return {f, brief};
    }

    /**
     * @brief Arbitrary pair potential between atom types using splines
     *
     * This will take the given pair-potential(s) and generate
     * splined versions for all possible atom type combinations.
     * If more than one potential is given, these will be added
     * together, then tabulated.
     *
     * ~~~~
     * "pairpotentialmap" : {
     *     "spline"  : { "rmin":1e-6, "rmax":100, "utol":0.01 },
     *     "default" : {
     *         "coulomb" : { "epsr":2.0 },
     *         "lennardjones" : {}
     *     },
     *     "Na Cl" : { "coulomb" : { "epsr":2.0 } }
     * }
     * ~~~~
     *
     * All pair-potentials are tabulated in constructor
     */
    template<typename Ttabulator=Tabulate::Andrea<double>, typename Tparticle=PointParticle>
    class PotentialMapSpline : public PairPotentialBase
    {
    private:
        typedef opair<int> Tpair;
        std::map<string, std::set<Tpair> > nfo;
        typedef typename Ttabulator::data Tdata;
        PairMatrix <Tdata> m;
        Ttabulator tab;
        bool verbose;
        double rmin, rmax;

        string _brief() override { return string("splined pair-potentials"); }

    public:
        PotentialMapSpline( Tmjson &in, const string sec = "pairpotentialmap" ) : PairPotentialBase(sec)
        {

            PairPotentialBase::name = "pairpotentialmap";

            auto j = in[sec]["spline"];

            rmin = j["rmin"] | 1.0;
            rmax = j["rmax"] | 100.0;
            tab.setRange(rmin, rmax);
            tab.setTolerance(
                j["utol"] | 0.01,
                j["ftol"] | -1,
                j["umaxtol"] | -1,
                j["fmaxtol"] | -1);

            verbose = in[sec]["verbose"] | false;

            m.resize(atom.size());

            // loop over all pairs in json entry
            for ( auto i = in[sec].begin(); i != in[sec].end(); ++i )
            {
                auto v = textio::words2vec<string>(i.key());
                if ( v.size() == 2 )
                {
                    Tpair p(
                        atom[v.at(0)].id,
                        atom[v.at(1)].id);
                    auto d = mixPairPotential(i.value(), p.first, p.second);
                    m.set(p.first, p.second, tab.generate(d.first));
                    nfo[d.second].insert(p);
                }
            }

            // for all non-assigned pairs, fallback to default potential
            for ( auto &i : atom )
                for ( auto &j : atom )
                    if ( i.id > 0 && j.id > 0 )
                        if ( m(i.id, j.id).empty())
                        {
                            auto d = mixPairPotential(in[sec]["default"], i.id, j.id);
                            m.set(i.id, j.id, tab.generate(d.first));
                            nfo[d.second].insert(Tpair(i.id, j.id));
                        }

            save("pairpotentials.dat", rmin, rmax);
        }

        double operator()( const Tparticle &a, const Tparticle &b, double r2 )
        {
            if ( r2 < rmin * rmin )
                return pc::infty;
            if ( r2 > rmax * rmax )
                return 0;
            assert(!m(a.id, b.id).empty());
            return tab.eval(m(a.id, b.id), r2);
        }

        void save( const string &file, double rmin = 1e-3, double rmax = 100 )
        {
            std::ofstream f(file.c_str());
            Tparticle a, b;
            for ( size_t i = 0; i < m.size(); i++ )
                for ( size_t j = 0; j < m.size(); j++ )
                    if ( j >= i )
                        if ( !m(i, j).empty())
                        {
                            f << "# " << atom[i].name << "-" << atom[j].name << endl;
                            a = atom[i];
                            b = atom[j];
                            for ( double r = rmin; r <= rmax; r += 0.1 )
                                f << r << " " << (*this)(a, b, r * r) << endl;
                            f << endl;
                        }
            f.close();
        }

        string info( char w = 0 ) override
        {
            std::ostringstream o;
            o << tab.info() << endl;
            for ( auto &i : nfo )
            {
                o << i.first << ":\n";
                for ( auto &p : i.second )
                    o << atom[p.first].name + "-" + atom[p.second].name << " ";
                o << endl << endl;
            }
            return o.str();
        }

    };

    /**
     * @brief Custom tabulated potentials between specific particle types
     *
     * If the pair is not recognized, i.e. not added with the
     * `add()` function, the `Tdefault` pair potential is used.
     * If the pair is found then a tabulation will be used.
     */
    template<typename Tdefault, typename Ttabulator=Tabulate::Andrea<double>, typename Tparticle=PointParticle>
    class PotentialMapTabulated : public PotentialMap<Tdefault>
    {
    private:
        typedef opair<int> Tpair;
        typedef PotentialMap <Tdefault> base;
        double rmin2, rmax2;
        int print;
        Ttabulator tab;
        std::map<Tpair, typename Ttabulator::data> mtab;

    public:
        PotentialMapTabulated( InputMap &in ) : base(in)
        {
            rmin2 = in.get<double>("tab_rmin", 1.0);
            rmax2 = in.get<double>("tab_rmax", 100.0);
            rmin2 = rmin2 * rmin2;
            rmax2 = rmax2 * rmax2;
            tab.setRange(
                in.get<double>("tab_rmin", 1.0),
                in.get<double>("tab_rmax", 100.0));
            tab.setTolerance(
                in.get<double>("tab_utol", 0.01),
                in.get<double>("tab_ftol", -1),
                in.get<double>("tab_umaxtol", -1),
                in.get<double>("tab_fmaxtol", -1));
            print = in.get<int>("tab_print", 0);
        }

        double operator()( const Tparticle &a, const Tparticle &b, double r2 )
        {
            auto ab = Tpair(a.id, b.id);
            auto it = mtab.find(ab);
            if ( it != mtab.end())
            {
                if ( r2 < it->second.rmax2 )
                    if ( r2 > it->second.rmin2 )
                        return tab.eval(it->second, r2);
                return base::m[ab](a, b, r2); // fall back to original
            }
            return Tdefault::operator()(a, b, r2); // fall back to default
        }

        template<class Tpairpot>
        void add( AtomData::Tid id1, AtomData::Tid id2, Tpairpot pot )
        {
            Tparticle a, b;
            a = atom[id1];
            b = atom[id2];
            base::add(a.id, b.id, pot);
            std::function<double( double )> f = [=]( double r2 ) { return Tpairpot(pot)(a, b, r2); };
            mtab[Tpair(id1, id2)] = tab.generate(f);
        }

        std::string info( char w = 20 )
        {
            using namespace Faunus::textio;
            std::ostringstream o(base::info(w));
            o << tab.info(w) << std::endl;
            for ( auto &i : mtab )
            {
                auto ab = Tpair(i.first.first, i.first.second);
                auto it = mtab.find(ab);
                o << pad(SUB,
                         w,
                         "Nbr of elements in table (" + atom[i.first.first].name + "<->" + atom[i.first.second].name
                             + "): ")
                  << it->second.r2.size() << endl;
            }
            o << endl;
            if ( print == 1 )
                print_tabulation();
            return o.str();
        }

        void print_tabulation( int n = 1000 )
        {
            for ( auto &i : mtab )
            {
                auto ab = Tpair(i.first.first, i.first.second);
                auto it = mtab.find(ab);

                Tparticle a, b;
                a = atom[i.first.first];
                b = atom[i.first.second];

                std::ofstream
                    ff1(std::string(atom[i.first.first].name + "." + atom[i.first.second].name + ".real.dat").c_str());
                ff1.precision(10);

                std::ofstream
                    ff2(std::string(atom[i.first.first].name + "." + atom[i.first.second].name + ".tab.dat").c_str());
                ff2.precision(10);

                double max = it->second.r2.at(it->second.r2.size() - 2);
                double min = it->second.rmin2;
                double dr = (max - min) / (double) n;
                for ( int j = 1; j < n; j++ )
                {
                    double r2 = min + dr * ((double) j);
                    ff1 << sqrt(r2) << " " << base::m[ab](a, b, r2) << endl;
                    ff2 << sqrt(r2) << " " << tab.eval(it->second, r2) << endl;
                }
            }
        }
    };

  }//Potential namespace
#endif
}//Faunus namespace
#endif
