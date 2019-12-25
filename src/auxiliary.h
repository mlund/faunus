#pragma once
#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <regex>
#include <chrono>

#include "average.h"

/**
 * @file auxiliary.h
 *
 * This file contains auxiliary functionality that
 * have no dependencies other than STL and can hence
 * be copied to other projects.
 */
namespace Faunus
{
    template<typename T=double>
        inline int anint( T x ) {
            return int(x > 0.0 ? x + 0.5 : x - 0.5);
        } //!< Round to int

    /**
     * @brief Iterate over pairs in container, call a function on the elements, and sum results
     * @tparam T Floating point type. Default: `double)`
     * @param begin Begin iterator
     * @param end End iterator
     * @param f Function to apply to the pair
     * @param aggregator Function to aggregate the result from each pair. Default: `std::plus<T>`
     */
    template <typename Titer, typename Tfunction, typename T = double, typename Taggregate_function>
    T for_each_unique_pair(Titer begin, Titer end, Tfunction f, Taggregate_function aggregator = std::plus<T>()) {
        T x = T();
        for (auto i = begin; i != end; ++i)
            for (auto j = i; ++j != end;)
                x = aggregator(x, f(*i, *j));
        return x;
    }
#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] for_each_pair")
    {
        int x;
        std::vector<int> a = {1,2,3};
        x = for_each_unique_pair(
            a.begin(), a.end(), [](int i, int j) { return i * j; }, std::plus<int>());
        CHECK(x==2+3+6);
        a.resize(1);
        x = for_each_unique_pair(
            a.begin(), a.end(), [](int i, int j) { return i * j; }, std::plus<int>());
        CHECK(x==0);
    }
#endif


    /** @brief Erase from container `a` all values found in container `b` */
    template<typename T>
        T erase_range( T a, const T &b )
        {
            a.erase(
                    std::remove_if(a.begin(), a.end(),
                        [b]( typename T::value_type i ) { return std::find(b.begin(), b.end(), i) != b.end(); }),
                    a.end());
            return a;
        }

    /**
     * @brief Ordered pair where `first<=second`
     *
     * Upon construction, the smallest element is placed in `first`
     * so that `opair<int>(i,j)==opair<int>(j,i)` is always true.
     *
     * @todo Add std::pair copy operator
     */
    template<class T>
        struct opair : public std::pair<T, T>
    {
        typedef std::pair<T, T> base;

        opair() {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-overflow"

        opair( const T &a, const T &b ) : base(a, b)
        {
            if ( a > b )
                std::swap(base::first, base::second);
        }

#pragma GCC diagnostic pop

        bool find( const T &i ) const
        {
            assert(base::first <= base::second);
            if ( i != base::first )
                if ( i != base::second )
                    return false;
            return true;
        }
    };

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] opair")
    {
        opair<int> a = {1,2}, b = {2,1};
        CHECK( (a.first == 1 && a.second == 2) );
        CHECK( (b.first == 1 && b.second == 2) );
        CHECK( a == b );
        CHECK( a.find(1) );
        CHECK( a.find(2) );
        CHECK( a.find(3)==false );
    }
#endif

    /**
     * @brief Approximation of erfc-function
     * @param x Value for which erfc should be calculated
     * @details Reference for this approximation is found in Abramowitz and Stegun,
     *          Handbook of mathematical functions, eq. 7.1.26
     *
     * @f[
     *     \erf(x) = 1 - (a_1t + a_2t^2 + a_3t^3 + a_4t^4 + a_5t^5)e^{-x^2} + \epsilon(x)
     * @f]
     * @f[
     *     t = \frac{1}{1 + px}
     * @f]
     * @f[
     *     |\epsilon(x)| \le 1.5\times 10^{-7}
     * @f]
     *
     * @warning Needs modification if x < 0
     */
    template<typename T>
        T erfc_x( T x )
        {
            static_assert(std::is_floating_point<T>::value, "type must be floating point");
            assert(x >= 0 && "x cannot be negative");
            T t = 1.0 / (1.0 + 0.3275911 * x);
            const T a1 = 0.254829592;
            const T a2 = -0.284496736;
            const T a3 = 1.421413741;
            const T a4 = -1.453152027;
            const T a5 = 1.061405429;
            return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))) * std::exp(-x * x);
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] erfc_x")
    {
        double infty = std::numeric_limits<double>::infinity();
        using doctest::Approx;
        CHECK( erfc_x(infty) == Approx(0));
        CHECK( std::erfc(-infty) == Approx(2-erfc_x(infty)));
        CHECK( erfc_x(0.0) == Approx(1.0) );
        CHECK( 2-erfc_x(0.2) == Approx(std::erfc(-0.2)));
        CHECK( erfc_x(0.2) == Approx( std::erfc(0.2) ));
    }
#endif

    /**
     * @brief Approximate 1 - erfc_x
     * @param x Value for which erf should be calculated
     */
    template<typename T>
        T erf_x( T x ) { return (1 - erfc_x(x)); }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"

    /**
     * @brief Fast inverse square-root approximation
     *
     * Modified to work with both double and float and with one (less precise) or
     * two (more precise) iterations. Template conditionals should be optimized out
     * at compile time.
     *
     * @note Code comments supposedly from the original Quake III Arena code
     */
    template <typename T, char iterations = 2> inline T inv_sqrt(T x) {
        static_assert(std::is_floating_point<T>::value, "T must be floating point");
        static_assert(iterations == 1 or iterations == 2, "itarations must equal 1 or 2");
        typedef typename std::conditional<sizeof(T) == 8, std::int64_t, std::int32_t>::type Tint;
        T y = x;
        T x2 = y * 0.5;
        Tint i = *(Tint *)&y;                                              // evil floating point bit level hacking
        i = (sizeof(T) == 8 ? 0x5fe6eb50c7b537a9 : 0x5f3759df) - (i >> 1); // what the fuck?
        y = *(T *)&i;
        y = y * (1.5 - (x2 * y * y)); // 1st iteration
        if (iterations > 1)
            y = y * (1.5 - (x2 * y * y)); // 2nd iteration, this can be removed
        return y;
    }

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE_TEMPLATE("inv_sqrt", T, double, float) {
        std::vector<T> vals = {0.23, 3.3, 10.2, 100.45, 512.06};
        for (auto x : vals)
            CHECK(inv_sqrt<T>(x) == doctest::Approx(1.0 / std::sqrt(x)));
    }
#endif

    /**
     * @brief n'th integer power of float
     *
     * On GCC/Clang this will use the fast `__builtin_powi` function.
     *
     * See also:
     * - https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp
     * - https://martin.ankerl.com/2007/10/04/optimized-pow-approximation-for-java-and-c-c/
     */
        template <class T> inline constexpr T powi(T x, unsigned int n) {
#if defined(__GNUG__)
            return __builtin_powi(x, n);
#else
            return n > 0 ? x * powi(x, n - 1) : 1;
#endif
        }
#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] powi") {
            using doctest::Approx;
            double x = 3.1;
            CHECK(powi(x, 0) == Approx(1));
            CHECK(powi(x, 1) == Approx(x));
            CHECK(powi(x, 2) == Approx(x * x));
            CHECK(powi(x, 4) == Approx(x * x * x * x));
        }
#endif

        /**
         * @brief Approximate exp() function
         * @note see [Cawley 2000](http://dx.doi.org/10.1162/089976600300015033)
         * @warning Does not work in big endian systems, nor on gcc
         *
         * Update 2019: http://www.federicoperini.info/wp-content/uploads/FastExp-Final.pdf
         */
        template <class Tint = std::int32_t> double exp_cawley(double y) {
            static_assert(2 * sizeof(Tint) == sizeof(double), "Approximate exp() requires 4-byte integer");
            union
            {
                double d;
                struct { Tint j, i; } n;  // little endian
                //struct { int i, j; } n;  // bin endian
            } eco;
            eco.n.i = 1072632447 + (Tint) (y * 1512775.39519519);
            eco.n.j = 0;
            return eco.d;
        }

    inline double exp_untested( double y )
    {
        typedef std::int32_t Tint;
        static_assert(2 * sizeof(Tint) == sizeof(double),
                "Approximate exp() requires 4-byte integer");
        double d(0);
        *((Tint *) (&d) + 0) = 0;
        *((Tint *) (&d) + 1) = (Tint) (1512775 * y + 1072632447);
        return d;
    }

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] exp_cawley")
    {
        double infty = std::numeric_limits<double>::infinity();
        using doctest::Approx;
        WARN( exp_cawley(-infty) == Approx(0)); // clang=OK; GCC=not OK
        CHECK( exp_cawley(0) == Approx(0.9710078239));
        CHECK( exp_cawley(2) == Approx(7.3096199036));
        CHECK( exp_cawley(-2) == Approx(0.13207829));
    }

    TEST_CASE("[Faunus] exp_untested")
    {
        double infty = std::numeric_limits<double>::infinity();
        using doctest::Approx;
        CHECK( exp_untested(-infty) == Approx(0)); // clang=OK; GCC=not OK
        CHECK( exp_untested(0) == Approx(0.9710078239));
        CHECK( exp_untested(2) == Approx(7.3096199036));
        CHECK( exp_untested(-2) == Approx(0.13207829));
    }
#endif

#pragma GCC diagnostic pop

        /**
         * @brief Evaluate n'th degree Legendre polynomial
         *
         * Example:
         * @code
         * Legendre<float> l(10);
         * auto P = l.eval(1.3)
         * std::cout << P[3]; --> third order value
         * @endcode
         *
         * @author Mikael Lund
         * @date Canberra 2005-2006
         * @note Since C++17 there's `std::legendre` but this seems more efficient
         *       if a range of degrees are needed
         */
        template <typename T = double> class Legendre {
          private:
            size_t n;         //!< Maximum Legendre order
            std::vector<T> P; //!< Legendre terms stored here
            std::vector<T> y; //!< Lookup table for 1+1/i (overkill?)
          public:
            /** @brief Construct w. polynomial order>=0 */
            Legendre(size_t max_order) : n(max_order) {
                P.resize(n + 1);
                P[0] = 1.0;
                y.resize(n + 1);
                for (size_t i = 1; i < n; i++)
                    y[i] = 1.0 + 1.0 / T(i);
            }

            /** @brief Evaluate polynomials at x */
            const std::vector<T> &eval(T x) {
                if (n > 0) {
                    P[1] = x;
                    for (size_t i = 1; i < n; ++i) {
                        P[i + 1] = ((y[i] + 1.0) * x * P[i] - P[i - 1]) / y[i];
                    }
                }
                return P;
            }
        };
#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] Legendre") {
            using doctest::Approx;
            Legendre<double> l(3);
            double x = 2.2;
            auto P = l.eval(x);
            CHECK(P[0] == Approx(1));
            CHECK(P[1] == Approx(x));
            CHECK(P[2] == Approx(0.5 * (3 * x * x - 1)));
            CHECK(P[3] == Approx(0.5 * (5 * x * x * x - 3 * x)));
        }
#endif

    /**
     * @brief Round a floating point to integer for use with histogram binning.
     *
     * This will round to nearest integer, assuming a certain resolution
     * (default 1). This can be useful for binning floating point data
     * into a histogram or table of resolution `dx` (second argument).
     *
     * Example:
     *
     * ~~~~
     *     double x=21.3;
     *     to_bin(x);     // -> 21
     *     to_bin(x,2);   // -> 11
     *     to_bin(x,0.5); // -> 43
     * ~~~~
     *
     */
    template<class T, class Tint=int>
        Tint to_bin( T x, T dx = 1 ) {
            return (x < 0) ? Tint(x / dx - 0.5) : Tint(x / dx + 0.5);
        }

    /**
     * @brief Round and bin numbers for use with tables and histograms
     *
     * This will round a float to nearest value divisible by `dx` or
     * convert to an integer corresponding to a binning value. In the
     * latter case, a minimum value must be specified upon construction.
     *
     * Example:
     *
     * ~~~{.cpp}
     * Quantize<double> Q(0.5, -0.5);
     * Q = 0.61;
     * double x = Q;      // --> 0.5
     * int bin = Q;       // --> 2
     * bin = Q(-0.5);     // --> 0
     * x = Q.frombin(2);  // --> 0.5
     *
     * std::vector<bool> v(10);
     * v[ Q(0.61) ] = false; // 2nd element set to false
     * ~~~
     *
     */
    template<class Tfloat=double>
        class Quantize {
            private:
                Tfloat xmin, dx, x;
            public:
                /**
                 * @brief Constructor
                 * @param dx resolution
                 * @param xmin minimum value if converting to integral type (binning)
                 */
                Quantize(Tfloat dx, Tfloat xmin=0) : xmin(xmin), dx(dx) {}

                /** @brief Assigment operator */
                Quantize& operator=(Tfloat val) {
                    assert( val>=xmin );
                    if (val >= 0)
                        x = int(val / dx + 0.5) * dx;
                    else
                        x = int(val / dx - 0.5) * dx;
                    assert(x>=xmin);
                    return *this;
                }

                Quantize& frombin(unsigned int i) {
                    x = i*dx + xmin;
                    return *this;
                }

                /** @brief Assignment with function operator */
                Quantize& operator()(Tfloat val) {
                    *this = val;
                    return *this;
                }

                /** @brief Implicit convertion to integral (bin) or float (rounded) */
                template<typename T> operator T()
                {
                    if (std::is_integral<T>::value)
                        return T( (x-xmin) / dx + 0.5);
                    return x;
                }
        };

    /**
     * @brief Use xy data from disk as function object
     *
     * This is a very basic structure to store xy data
     * in a map to be used as a function object. Linear
     * interpolation is used between data points. The lookup
     * complexity is `log(N)` and no particular spacing
     * between data points is expected. Upon loading, data
     * is sorted and possible duplicate points trimmed.
     *
     * Example:
     *
     *     InterpolTable<double> f("xy.dat");
     *     double y;
     *     y = f(15);  // --> 2.5
     *     y = f(-1);  // --> NaN
     *     y = f(21);  // --> NaN
     *
     * where the `xy.dat` file may look like this (comments not allowed!),
     *
     *      0.0   1.0
     *     10.0   2.0
     *     20.0   3.0
     *
     * @date Malmo 2014
     */
    template<typename T=double>
        class InterpolTable
        {
            private:
                typedef std::pair<T, T> Tpair; // xy data stored as pairs
                std::vector<Tpair> t;         // in a vector
            public:
                InterpolTable() {}

                InterpolTable( const std::string &filename ) { load(filename); }

                bool load( const std::string &filename )
                {
                    Tpair a;
                    t.clear();
                    std::ifstream in(filename);
                    if ( in )
                    {
                        while ( in >> a.first >> a.second )
                            t.push_back(a);
                        std::sort(t.begin(), t.end()); // sort
                        t.erase(std::unique(t.begin(), t.end()), t.end()); // remove duplicates
                        if ( !t.empty())
                            return true;
                    }
                    return false;
                }

                T xmin() const { return t.front().first; }

                T xmax() const { return t.back().first; }

                T operator()( T x ) const
                {
                    assert(!t.empty() && "Table is empty");
                    if ( x > t.back().first )
                        return std::numeric_limits<T>::quiet_NaN();
                    if ( x < t[0].first )
                        return std::numeric_limits<T>::quiet_NaN();
                    auto it = std::lower_bound(t.begin(), t.end(), Tpair(x, 0));
                    if ( it == t.begin())
                        return it->second;
                    auto it2 = it;
                    --it2;
                    return it2->second + (it->second - it2->second) * (x - it2->first) / (it->first - it2->first);
                }
        };

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
    template<class T, bool triangular=false>
        class PairMatrix {
            private:
                T __val; // default value when resizing
                std::vector<std::vector<T>> m;
            public:
                void resize(size_t n) {
                    m.resize(n);
                    for (size_t i=0; i<m.size(); i++)
                        if (triangular)
                            m[i].resize(i+1, __val);
                        else
                            m[i].resize(n, __val);
                }

                PairMatrix(size_t n=0, T val=T()) : __val(val) {
                    resize(n);
                }

                auto size() const { return m.size(); }

                inline const T &operator()(size_t i, size_t j) const {
                    if (triangular)
                        if (j>i)
                            std::swap(i,j);
                    assert(i < m.size());
                    assert(j < m[i].size());
                    return m[i][j];
                }

                void set(size_t i, size_t j, T val) {
                    if (j>i)
                        std::swap(i, j);
                    if (i>=m.size())
                        resize(i+1);
                    if (triangular==false)
                        m[j][i] = val;
                    m[i][j] = val;
                }
        };
#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] PairMatrix")
    {
        int i=2,j=3; // particle type, for example

        SUBCASE("full matrix") {
            PairMatrix<double,false> m;
            m.set(i,j,12.1);
            CHECK(m.size()==4);
            CHECK(m(i,j)==12.1);
            CHECK(m(i,j)==m(j,i));
            CHECK(m(0,2)==0);
            CHECK(m(2,0)==0);
        }

        SUBCASE("full matrix - default value") {
            PairMatrix<double,false> m(5,3.1);
            for (size_t i=0; i<5; i++)
                for (size_t j=0; j<5; j++)
                    CHECK(m(i,j)==3.1);
        }

        SUBCASE("triangular matrix - default value") {
            PairMatrix<double,true> m(5,3.1);
            for (size_t i=0; i<5; i++)
                for (size_t j=0; j<5; j++)
                    CHECK(m(i,j)==3.1);
        }

        SUBCASE("triangular matrix") {
            PairMatrix<double,true> m;
            m.set(i,j,12.1);
            CHECK(m.size()==4);
            CHECK(m(i,j)==12.1);
            CHECK(m(i,j)==m(j,i));
            CHECK(m(0,2)==0);
            CHECK(m(2,0)==0);
        }
    }
#endif

    template<typename Tcoeff=double, typename base=Eigen::Matrix<Tcoeff, Eigen::Dynamic, Eigen::Dynamic>>
        class Table : public base
    {
        private:
            typedef std::vector<double> Tvec;
            Tvec _bw, _lo, _hi;
            int _rows, _cols;
        public:
            Table( const Tvec &bw = {1, 1}, const Tvec &lo = {0, 0}, const Tvec &hi = {2, 2} ) {
                reInitializer(bw, lo, hi);
            }

            // required for assignment from Eigen::Matrix and Eigen::Array objects
            template<typename OtherDerived>
                Table(const Eigen::MatrixBase<OtherDerived>& other) : base(other) {}
            template<typename OtherDerived>
                Table& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
                    this->base::operator=(other);
                    return *this;
                }
            template<typename OtherDerived>
                Table(const Eigen::ArrayBase<OtherDerived>& other) : base(other) {}
            template<typename OtherDerived>
                Table& operator=(const Eigen::ArrayBase<OtherDerived>& other) {
                    this->base::operator=(other);
                    return *this;
                }

            void reInitializer( const Tvec &bw, const Tvec &lo, const Tvec &hi ) {
                assert( bw.size()==1 || bw.size()==2 );
                assert( bw.size()==lo.size() && lo.size()==hi.size() );
                _bw = bw;
                _lo = lo;
                _hi = hi;
                _rows = (_hi[0] - _lo[0]) / _bw[0] + 1.;
                if (bw.size()==2)
                    _cols = (_hi[1] - _lo[1]) / _bw[1] + 1;
                else
                    _cols = 1;
                base::resize(_rows, _cols);
                base::setZero();
            }

            void round( Tvec &v ) const {
                for ( Tvec::size_type i = 0; i != v.size(); ++i )
                    v[i] = (v[i] >= 0) ? int(v[i] / _bw[i] + 0.5) * _bw[i] : int(v[i] / _bw[i] - 0.5) * _bw[i];
            }

            void to_index(Tvec &v) const {
                for (Tvec::size_type i = 0; i != v.size(); ++i) {
                    v[i] = (v[i] >= 0) ? int(v[i] / _bw[i] + 0.5) : int(v[i] / _bw[i] - 0.5);
                    v[i] = v[i] - _lo[i] / _bw[i];
                }
                v.resize(2, 0);
            }

            Tcoeff &operator[](const Tvec &v) {
                return base::operator()(v[0], v[1]);
            }

            bool isInRange(const Tvec &v) const {
                bool b = true;
                for ( Tvec::size_type i = 0; i != v.size(); ++i )
                    b = b && v[i] >= _lo[i] && v[i] <= _hi[i];
                return b;
            }

            Tvec hist2buf( int ) const {
                Tvec sendBuf;
                for ( int i = 0; i < _cols; ++i )
                    for ( int j = 0; j < _rows; ++j )
                        sendBuf.push_back(base::operator()(j, i));
                return sendBuf;
            }

            void buf2hist(const Tvec &v) {
                assert(!v.empty());
                base::setZero();
                int p = v.size() / this->size(), n = 0;
                double nproc = p;
                while ( p-- > 0 )
                {
                    for ( int i = 0; i < _cols; ++i )
                        for ( int j = 0; j < _rows; ++j )
                        {
                            base::operator()(j, i) += v.at(n) / nproc;
                            ++n;
                        }
                }
            }

            base getBlock( const Tvec &v ) { // {xmin,xmax} or {xmin,xmax,ymin,ymax}
                Tvec w(4,0);
                switch (v.size()) {
                    case (1):w[0] = w[1] = (v[0] - _lo[0]) / _bw[0];
                             w[3] = _cols - 1;
                             break;
                    case (2):w[0] = (v[0] - _lo[0]) / _bw[0];
                             w[1] = (v[1] - _lo[0]) / _bw[0];
                             break;
                    case (3):w[0] = w[1] = _rows - 1;
                             w[2] = w[3] = _cols - 1;
                             break;
                    case (4):w[0] = (v[0] - _lo[0]) / _bw[0];
                             w[1] = (v[1] - _lo[0]) / _bw[0];
                             w[2] = (v[2] - _lo[1]) / _bw[1];
                             w[3] = (v[3] - _lo[1]) / _bw[1];
                }
                return this->block(w[0], w[2], w[1] - w[0] + 1, w[3] - w[2] + 1); // xmin,ymin,rows,cols
            }

            Tcoeff avg( const Tvec &v ) const {
                return this->getBlock(v).mean();
            }

            void save( const std::string &filename, Tcoeff scale = 1, Tcoeff translate = 0 ) const {
                Eigen::VectorXd v1(_cols + 1), v2(_rows + 1);
                v1(0) = v2(0) = base::size();
                for ( int i = 1; i != _cols + 1; ++i )
                {
                    v1(i) = (i - 1) * _bw[1] + _lo[1];
                }
                for ( int i = 1; i != _rows + 1; ++i )
                    v2(i) = (i - 1) * _bw[0] + _lo[0];
                base m(_rows + 1, _cols + 1);
                m.leftCols(1) = v2;
                m.topRows(1) = v1.transpose();
                m.bottomRightCorner(_rows, _cols) = *this;
                if ( scale != 1 )
                    m.bottomRightCorner(_rows, _cols) *= scale;
                if ( translate != 0 )
                    m.bottomRightCorner(_rows, _cols) += base::Constant(_rows, _cols, translate);
                std::ofstream f(filename.c_str());
                if ( _cols == 1 )
                    f << "#";
                f.precision(10);
                if ( f )
                    f << m;
            }

            void saveRow( const std::string &filename, const Tvec &v, Tcoeff scale = 1, Tcoeff translate = 0 )
            {
                if ( this->isInRange(v))
                {
                    auto b = this->getBlock(v);
                    int size = b.size();
                    Eigen::VectorXd w(size);
                    for ( int i = 0; i != size; ++i )
                        w(i) = i * _bw[1] + _lo[1];
                    base m(size, 2);
                    m.leftCols(1) = w;
                    m.bottomRightCorner(size, 1) = b.transpose();
                    if ( scale != 1 )
                        m.bottomRightCorner(size, 1) *= scale;
                    if ( translate != 0 )
                        m.bottomRightCorner(size, 1) += base::Constant(size, 1, translate);
                    std::ofstream f(filename.c_str());
                    f.precision(10);
                    if ( f )
                        f << m;
                }
            }

            void load( const std::string &filename )
            {
                std::ifstream f(filename.c_str());
                if ( f )
                {
                    int i = 0, j = -1;
                    std::string line;
                    getline(f, line);
                    while ( getline(f, line))
                    {
                        if ( i > _rows - 1 || j > _cols - 1 )
                        {
                            std::cerr << "Error: input file '" + filename + "' is larger than expected\n";
                            exit(1);
                        }
                        j = -1;
                        std::istringstream iss(line);
                        Tcoeff a, b;
                        iss >> a;
                        while ( iss >> b )
                            base::operator()(i, ++j) = b;
                        ++i;
                    }
                    if ( i != _rows || j != _cols - 1 )
                    {
                        std::cerr << "Error: input file '" + filename + "' is smaller than expected\n";
                        exit(1);
                    }
                }
            }
    };

    /**
     * @brief General class for handling 2D tables - xy data, for example.
     * @date Lund 2011
     * @note `Tx` is used as the `std::map` key and which may be
     * problematic due to direct floating point comparison (== operator).
     * We have not experienced any issues with this, though. This uses
     * `std::map` and table lookup is of complexity logarithmic with N.
     */
    template<typename Tx, typename Ty>
        class Table2D
        {
            protected:
                typedef std::map<Tx, Ty> Tmap;

                Ty count()
                {
                    Ty cnt = 0;
                    for ( auto &m : map )
                        cnt += m.second;
                    return cnt;
                }

                Tx dx;
                Tmap map;
                std::string name;
            private:
                Tx round( Tx x ) { return (x >= 0) ? int(x / dx + 0.5) * dx : int(x / dx - 0.5) * dx; }

                virtual double get( Tx x ) { return operator()(x); }

            public:
                enum type { HISTOGRAM, XYDATA };
                type tabletype;

                /** @brief Sum of all y values (same as `count()`) */
                Ty sumy() const
                {
                    Ty sum = 0;
                    for ( auto &m : map )
                        sum += m.second;
                    return sum;
                }

                /**
                 * @brief Constructor
                 * @param resolution Resolution of the x axis
                 * @param key Table type: HISTOGRAM or XYDATA
                 */
                Table2D( Tx resolution = 0.2, type key = XYDATA )
                {
                    tabletype = key;
                    setResolution(resolution);
                }

                /** @brief Convert to map */
                std::map<std::string, std::vector < double>> to_map() {
                    std::map<std::string, std::vector < double>>
                        m;
                    m["x"].reserve(map.size());
                    m["y"].reserve(map.size());
                    for ( auto &i : map )
                    {
                        m["x"].push_back(i.first);
                        m["y"].push_back(get(i.first));
                    }
                    return m;
                }

                void clear() { map.clear(); }

                void setResolution( Tx resolution )
                {
                    assert(resolution > 0);
                    dx = resolution;
                    map.clear();
                }

                void setResolution( std::vector<Tx> &resolution )
                {
                    assert(resolution[0] > 0);
                    dx = resolution[0];
                    map.clear();
                }

                virtual ~Table2D() {}

                /** @brief Access operator - returns reference to y(x) */
                Ty &operator()( Tx x )
                {
                    return map[round(x)];
                }

                /** @brief Access operator - returns reference to y(x) */
                Ty &operator()( std::vector<Tx> &x )
                {
                    return map[round(x[0])];
                }

                /** @brief Find key and return corresponding value otherwise zero*/
                Ty find( std::vector<Tx> &x )
                {
                    Ty value = 0;
                    auto it = map.find(round(x[0]));
                    if ( it != map.end())
                        value = it->second;
                    return value;
                }

                /** @brief Save table to disk */
                template<class T=double>
                    void save( const std::string &filename, T scale = 1, T translate = 0 )
                    {
                        if ( tabletype == HISTOGRAM )
                        {
                            if ( !map.empty())
                                map.begin()->second *= 2;   // compensate for half bin width
                            if ( map.size() > 1 )
                                (--map.end())->second *= 2; // -//-
                        }

                        if ( !map.empty())
                        {
                            std::ofstream f(filename.c_str());
                            f.precision(10);
                            if ( f )
                            {
                                for ( auto &m : map )
                                    f << m.first << " " << (m.second + translate) * scale << "\n";
                            }
                        }

                        if ( tabletype == HISTOGRAM )
                        {
                            if ( !map.empty())
                                map.begin()->second /= 2;   // restore half bin width
                            if ( map.size() > 1 )
                                (--map.end())->second /= 2; // -//-
                        }
                    }

                /** @brief Save normalized table to disk */
                template<class T=double>
                    void normSave( const std::string &filename )
                    {
                        if ( tabletype == HISTOGRAM )
                        {
                            if ( !map.empty())
                                map.begin()->second *= 2;   // compensate for half bin width
                            if ( map.size() > 1 )
                                (--map.end())->second *= 2; // -//-
                        }

                        if ( !map.empty())
                        {
                            std::ofstream f(filename.c_str());
                            f.precision(10);
                            Ty cnt = count() * dx;
                            if ( f )
                            {
                                for ( auto &m : map )
                                    f << m.first << " " << m.second / cnt << "\n";
                            }
                        }

                        if ( tabletype == HISTOGRAM )
                        {
                            if ( !map.empty())
                                map.begin()->second /= 2;   // restore half bin width
                            if ( map.size() > 1 )
                                (--map.end())->second /= 2; // -//-
                        }
                    }

                /** @brief Sums up all previous elements and saves table to disk */
                template<class T=double>
                    void sumSave( std::string filename, T scale = 1 )
                    {
                        if ( tabletype == HISTOGRAM )
                        {
                            if ( !map.empty())
                                map.begin()->second *= 2;   // compensate for half bin width
                            if ( map.size() > 1 )
                                (--map.end())->second *= 2; // -//-
                        }

                        if ( !map.empty())
                        {
                            std::ofstream f(filename.c_str());
                            f.precision(10);
                            if ( f )
                            {
                                Ty sum_t = 0.0;
                                for ( auto &m : map )
                                {
                                    sum_t += m.second;
                                    f << m.first << " " << sum_t * scale << "\n";
                                }
                            }
                        }

                        if ( tabletype == HISTOGRAM )
                        {
                            if ( !map.empty())
                                map.begin()->second /= 2;   // restore half bin width
                            if ( map.size() > 1 )
                                (--map.end())->second /= 2; // -//-
                        }
                    }

                const Tmap& getMap() const
                {
                    return map;
                }

                Tmap& getMap()
                {
                    return map;
                }

                Tx getResolution()
                {
                    return dx;
                }

                /*! Returns average */
                Tx mean()
                {
                    assert(!map.empty());
                    Tx avg = 0;
                    for ( auto &m : map )
                        avg += m.first * m.second;
                    return avg / count();
                }

                /*! Returns standard deviation */
                Tx std()
                {
                    assert(!map.empty());
                    Tx std2 = 0;
                    Tx avg = mean();
                    for ( auto &m : map )
                        std2 += m.second * (m.first - avg) * (m.first - avg);
                    return sqrt(std2 / count());
                }

                /*! Returns iterator of minumum y */
                typename Tmap::const_iterator min()
                {
                    assert(!map.empty());
                    Ty min = std::numeric_limits<Ty>::max();
                    typename Tmap::const_iterator it;
                    for ( auto m = map.begin(); m != map.end(); ++m )
                        if ( m->second < min )
                        {
                            min = m->second;
                            it = m;
                        }
                    return it;
                }

                /*! Returns iterator of maximum y */
                typename Tmap::const_iterator max()
                {
                    assert(!map.empty());
                    Ty max = std::numeric_limits<Ty>::min();
                    typename Tmap::const_iterator it;
                    for ( auto m = map.begin(); m != map.end(); ++m )
                        if ( m->second > max )
                        {
                            max = m->second;
                            it = m;
                        }
                    return it;
                }

                /*! Returns x at minumum x */
                Tx minx()
                {
                    assert(!map.empty());
                    Tx x = 0;
                    for ( auto &m : map )
                    {
                        x = m.first;
                        break;
                    }
                    return x;
                }

                /*! Returns average in interval */
                Ty avg( const std::vector<Tx> &limits )
                {
                    Ty avg = 0;
                    int cnt = 0;
                    assert(!map.empty());
                    for ( auto &m : map )
                    {
                        if ( m.first >= limits[0] && m.first <= limits[1] )
                        {
                            avg += m.second;
                            ++cnt;
                        }
                    }
                    if ( cnt > 0 )
                        avg /= cnt;
                    return avg;
                }

                /**
                 * @brief Convert table2D to vector of floats
                 */
                std::vector<double> hist2buf( int &size )
                {
                    std::vector<double> sendBuf;
                    assert(!map.empty());
                    for ( auto &m : map )
                    {
                        sendBuf.push_back(m.first);
                        sendBuf.push_back(m.second);
                    }
                    sendBuf.resize(size, -1);
                    return sendBuf;
                }

                /**
                 * @brief Convert vector of floats to table2D
                 */
                void buf2hist( std::vector<double> &v )
                {
                    this->clear();
                    assert(!v.empty());
                    std::map<double, Average<double>> all;
                    for ( int i = 0; i < int(v.size()) - 1; i += 2 )
                        if ( v.at(i + 1) != -1 )
                            all[v.at(i)] += v.at(i + 1);
                    for ( auto &m : all )
                        this->operator()(m.first) = m.second.avg();
                }

                /**
                 * @brief Load table from disk
                 * @note The first line - used for comments - is ignored.
                 * @todo Implement end bin compensation as in the save()
                 * function when loading HISTOGRAMs
                 */
                bool load( const std::string &filename )
                {
                    std::ifstream f(filename.c_str());
                    if ( f )
                    {
                        map.clear();
                        while ( !f.eof())
                        {
                            Tx x;
                            double y;
                            f >> x >> y;
                            operator()(x) = y;
                        }
                        if ( tabletype == HISTOGRAM )
                        {
                            if ( !map.empty())
                                map.begin()->second /= 2;   // restore half bin width
                            if ( map.size() > 1 )
                                (--map.end())->second /= 2; // -//-
                        }
                        return true;
                    }
                    return false;
                }

                /**
                 * @brief Convert table to matrix
                 */
                Eigen::MatrixXd tableToMatrix()
                {
                    assert(!this->map.empty() && "Map is empty!");
                    Eigen::MatrixXd table(2, map.size());
                    table.setZero();
                    int I = 0;
                    for ( auto &m : this->map )
                    {
                        table(0, I) = m.first;
                        table(1, I) = m.second;
                        I++;
                    }
                    return table;
                }
        };

    /**
     * @brief Subtract two tables
     */
    template<class Tx, class Ty>
        Table2D<Tx, Ty> operator-( Table2D<Tx, Ty> &a, Table2D<Tx, Ty> &b )
        {
            assert(a.tabletype == b.tabletype && "Table a and b needs to be of same type");
            Table2D<Tx, Ty> c(std::min(a.getResolution(), b.getResolution()), a.tabletype);
            auto a_map = a.getMap();
            auto b_map = b.getMap();

            if ( a.tabletype == Table2D<Tx, Ty>::HISTOGRAM )
            {
                if ( !a_map.empty())
                    a_map.begin()->second *= 2;   // compensate for half bin width
                if ( a_map.size() > 1 )
                    (--a_map.end())->second *= 2; // -//-
                if ( !b_map.empty())
                    b_map.begin()->second *= 2;   // compensate for half bin width
                if ( b_map.size() > 1 )
                    (--b_map.end())->second *= 2; // -//-
            }

            for ( auto &m1 : a_map )
            {
                for ( auto &m2 : b_map )
                {
                    c(m1.first) = m1.second - m2.second;
                    break;
                }
            }

            if ( a.tabletype == Table2D<Tx, Ty>::HISTOGRAM )
            {
                if ( !a_map.empty())
                    a_map.begin()->second /= 2;   // compensate for half bin width
                if ( a_map.size() > 1 )
                    (--a_map.end())->second /= 2; // -//-
                if ( !b_map.empty())
                    b_map.begin()->second /= 2;   // compensate for half bin width
                if ( b_map.size() > 1 )
                    (--b_map.end())->second /= 2; // -//-
            }
            return c;
        }

    /**
     * @brief Addition two tables
     */
    template<class Tx, class Ty>
        Table2D<Tx, Ty> operator+( Table2D<Tx, Ty> &a, Table2D<Tx, Ty> &b )
        {
            assert(a.tabletype == b.tabletype && "Table a and b needs to be of same type");
            Table2D<Tx, Ty> c(std::min(a.getResolution(), b.getResolution()), a.tabletype);
            auto a_map = a.getMap();
            auto b_map = b.getMap();

            if ( a.tabletype == Table2D<Tx, Ty>::HISTOGRAM )
            {
                if ( !a_map.empty())
                    a_map.begin()->second *= 2;   // compensate for half bin width
                if ( a_map.size() > 1 )
                    (--a_map.end())->second *= 2; // -//-
                if ( !b_map.empty())
                    b_map.begin()->second *= 2;   // compensate for half bin width
                if ( b_map.size() > 1 )
                    (--b_map.end())->second *= 2; // -//-
            }

            for ( auto &m : a_map )
            {
                c(m.first) += m.second;
            }
            for ( auto &m : b_map )
            {
                c(m.first) += m.second;
            }

            if ( a.tabletype == Table2D<Tx, Ty>::HISTOGRAM )
            {
                if ( !a_map.empty())
                    a_map.begin()->second /= 2;   // compensate for half bin width
                if ( a_map.size() > 1 )
                    (--a_map.end())->second /= 2; // -//-
                if ( !b_map.empty())
                    b_map.begin()->second /= 2;   // compensate for half bin width
                if ( b_map.size() > 1 )
                    (--b_map.end())->second /= 2; // -//-
            }

            return c;
        }

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
    template<typename Tx=double, typename Ty=Tx, bool centerbin=false>
        class Equidistant2DTable {
            private:
                Tx _dxinv, _xmin;
                int offset;
                std::vector<Ty> vec;

                inline int to_bin(Tx x) const {
                    if (centerbin)
                        return (x<0) ? int(x*_dxinv-0.5) : int(x*_dxinv+0.5);
                    else
                        return std::floor( (x-_xmin)*_dxinv );
                } //!< x-value to vector index

                Tx from_bin(int i) const {
                    assert(i>=0);
                    return i/_dxinv + _xmin;
                } //!< vector index to x-value

            public:

                class iterator {
                    public:
                        typedef Equidistant2DTable<Tx,Ty,centerbin> T;
                        iterator(T &tbl, size_t i): tbl(tbl), i(i) {}
                        iterator operator++() { ++i; return *this; }
                        bool operator!=(const iterator& other) const { return i != other.i; }
                        auto operator*() { return tbl[i]; }
                    protected:
                        T &tbl;
                        size_t i;
                }; // enable range-based for loops

                class const_iterator {
                    public:
                        typedef Equidistant2DTable<Tx,Ty,centerbin> T;
                        const_iterator(const T &tbl, size_t i): tbl(tbl), i(i) {}
                        const_iterator operator++() { ++i; return *this; }
                        bool operator!=(const const_iterator& other) const { return i != other.i; }
                        auto operator*() const { return tbl[i]; }
                    protected:
                        const T &tbl;
                        size_t i;
                }; // enable range-based for loops

                iterator begin() { return iterator(*this, 0); }
                iterator end() { return iterator(*this, size()); }
                const_iterator begin() const { return const_iterator(*this, 0); }
                const_iterator end() const { return const_iterator(*this, size()); }

                Equidistant2DTable() {}

                /**
                 * @brief Constructor
                 * @param dx x spacing
                 * @param xmin minimum x value
                 * @param xmax maximum x value (for more efficient memory handling, only)
                 */
                Equidistant2DTable(Tx dx, Tx xmin, Tx xmax=std::numeric_limits<Tx>::infinity()) {
                    setResolution(dx, xmin, xmax);
                }

                void setResolution(Tx dx, Tx xmin, Tx xmax=std::numeric_limits<Tx>::infinity()) {
                    _xmin = xmin;
                    _dxinv = 1/dx;
                    offset = -to_bin(_xmin);
                    if (xmax<std::numeric_limits<Tx>::infinity())
                        vec.reserve( to_bin(xmax) + offset );
                }

                double sumy() const {
                    double sum=0;
                    for (auto &i : vec)
                        sum += double(i);
                    return sum;
                } //!< sum all y-values

                const std::vector<Ty> &yvec() const { return vec; } //!< vector with y-values

                std::vector<Tx> xvec() const {
                    std::vector<Tx> v;
                    v.reserve( vec.size() );
                    for (size_t i=0; i<vec.size(); i++)
                        v.push_back( from_bin(i) );
                    return v;
                }

                void clear() { vec.clear(); }
                size_t size() const { return vec.size(); }
                bool empty() const { return vec.empty(); }
                Tx dx() const { return 1/_dxinv; }
                Tx xmin() const { return _xmin; } //!< minimum x-value

                Tx xmax() const {
                    return (vec.empty()) ? _xmin : from_bin(vec.size()-1);
                } //!< maximum stored x-value

                Ty& operator()(Tx x) {
                    assert(x>=_xmin);
                    int i = to_bin(x) + offset;
                    if (i>=vec.size())
                        vec.resize( i+1, Ty() );
                    return vec[i];
                } // return y value for given x

                std::pair<Tx,Ty&> operator[](size_t index) {
                    assert(index>=0);
                    assert(index<size());
                    return { from_bin(index), vec[index] };
                } // pair with x,y& value

                std::pair<Tx,const Ty&> operator[](size_t index) const {
                    assert(index>=0);
                    assert(index<size());
                    return { from_bin(index), vec[index] };
                } // const pair with x,y& value

                const Ty& operator()(Tx x) const {
                    assert(x>=_xmin);
                    int i = to_bin(x) + offset;
                    assert(i<vec.size());
                    return vec.at(i);
                } // return y value for given x

                // can be optinally used to customize streaming out, normalise etc.
                std::function<void(std::ostream&,Tx,Ty)> stream_decorator=nullptr;

                friend std::ostream &operator<<( std::ostream &o, const Equidistant2DTable<Tx,Ty,centerbin> &tbl ) {
                    o.precision(16);
                    for (auto d : tbl)
                        if (tbl.stream_decorator==nullptr)
                            o << d.first << " " << d.second << "\n";
                        else
                            tbl.stream_decorator(o, d.first, d.second);
                    return o;
                } // write to stream

                auto &operator<<( std::istream &in ) {
                    assert(stream_decorator==nullptr && "you probably don't want to load a decorated file");
                    clear();
                    Tx x;
                    Ty y;
                    std::string line;
                    while ( std::getline(in, line)) {
                        std::stringstream o(line);
                        while (o >> x) {
                            if ((x>=_xmin) and (y << o))
                                operator()(x) = y;
                            else
                                throw std::runtime_error("table load error: x smaller than xmin");
                        }
                    }
                    return *this;
                } // load from stream

        };

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Equidistant2DTable")
    {
        using doctest::Approx;

        SUBCASE("centered") {
            Equidistant2DTable<double, double, true> y(0.5, -3.0);
            CHECK( y.xmin() == Approx(-3.0) );
            CHECK( y.dx() == Approx(0.5) );
            y(-2.51) = 0.15;
            CHECK( y(-2.5) == Approx(0.15));
            y(-2.76) = 0.11;
            CHECK( y(-3) == Approx(0.11));
            y(0.4) = 0.3;
            CHECK( y(0.5) == Approx(0.3));
            y(1.3) = 0.5;
            CHECK( y(1.5) == Approx(0.5));
            CHECK( y.xmax() == Approx(1.5) );
        }

        SUBCASE("lower bound") {
            Equidistant2DTable<double> y(0.5, -3.0);
            CHECK( y.xmin() == Approx(-3.0) );
            CHECK( y.dx() == Approx(0.5) );
            y(-2.51) = 0.15;
            CHECK( y(-3.0) == Approx(0.15));
            y(-2.76) = 0.11;
            CHECK( y(-3) == Approx(0.11));
            y(0.4) = 0.3;
            CHECK( y(0.0) == Approx(0.3));
            y(1.3) = 0.5;
            CHECK( y(1.0) == Approx(0.5));
            CHECK( y.xmax() == Approx(1.0) );
        }

    }
#endif

    /**
     * @brief Timer for measuring relative time consumption
     *
     * Time t=0 is set upon construction whereafter combined `start()`/
     * `stop()` calls can be made multiple times. The result is
     * the fraction of total time, consumed in between start/stop calls.
     */
    template<typename Tunit = std::chrono::microseconds>
        class TimeRelativeOfTotal
        {
            private:
                Tunit delta;
                std::chrono::steady_clock::time_point t0, tx;
            public:
              TimeRelativeOfTotal() : delta(0) { t0 = std::chrono::steady_clock::now(); }

              operator bool() const { return delta.count() != 0 ? true : false; }

              void start() { tx = std::chrono::steady_clock::now(); }

              void stop() { delta += std::chrono::duration_cast<Tunit>(std::chrono::steady_clock::now() - tx); }

              double result() const {
                  auto now = std::chrono::steady_clock::now();
                  auto total = std::chrono::duration_cast<Tunit>(now - t0);
                  return delta.count() / double(total.count());
                }
        };

    /**
     * @brief Simple struct to measure duration of code
     *
     * Example:
     *
     * ~~~ cpp
     * Stopwatch w;
     * w.start();
     * // do something expensive...
     * w.stop(); // --> stdio
     * ~~~
     *
     * `stop()` takes an optional template parameter to specify
     * the time resolution. Default is milliseconds.
     */
    struct Stopwatch {
        using clock = std::chrono::high_resolution_clock;
        std::chrono::time_point<clock> beg, end;
        inline Stopwatch() { start(); }
        inline void start() { beg = clock::now(); }
        template<typename T=std::chrono::milliseconds>
            void stop(bool print=true) {
                end = clock::now();
                if (print)
                    std::cout << "duration: " << std::chrono::duration_cast<T>(end-beg).count() << std::endl;
            }
    };

    /** @brief Count number of white-space separated words in a string */
    inline size_t numWords( const std::string &s )
    {
        return std::distance(std::istream_iterator<std::string>(
                    std::istringstream(s) >> std::ws), std::istream_iterator<std::string>());
    }

    /**
     * @brief Convert whitespace separated words into vector of given type
     *
     * Example:
     *
     * ~~~~
     * auto v = textio::words2vec<double>( "0.2 1 100" );
     * for (auto i : v)
     *   cout << 2*i << " "; // -> 0.4 2 200
     * ~~~~
     */
    template<class T>
        std::vector<T> words2vec( const std::string &w )
        {
            std::vector<T> v(numWords(w));
            std::stringstream s(w);
            size_t i = 0;
            while ( i < v.size())
            {
                s >> v[i++];
            }
            return v;
        } // space separated string to vector

    template<class T>
        std::string vec2words( const std::vector<T> &v ) {
            std::ostringstream o;
            if (!v.empty()) {
                o << v.front();
                for (size_t i=1; i<v.size(); i++)
                    o << " " << v[i];
            }
            return o.str();
        } // vector to space separated string w values

    /** @brief Convert string to lower case */
    inline std::string lowercase( std::string s )
    {
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);
        return s;
    }

    /** @brief Uppercase first letter in string */
    inline std::string toupper_first( std::string s ) {
        if (!s.empty())
            s[0] = std::toupper( s[0] );
        return s;
    }

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Text manipulation")
    {
        CHECK( numWords("a b c") == 3 );
        CHECK( vec2words<double>({1.0, -1.2, 0}) == "1 -1.2 0" );
        CHECK( words2vec<double>("1 -1.2 0") == std::vector<double>({1.0, -1.2, 0}) );
        CHECK( lowercase("aBc") == "abc" );
        CHECK( toupper_first("abc") == "Abc" );
    }
#endif

    template <typename T> struct BasePointerVector {
        std::vector<std::shared_ptr<T>> vec; //!< Vector of shared pointers to base class

        auto begin() noexcept { return vec.begin(); }
        auto begin() const noexcept { return vec.begin(); }
        auto end() noexcept { return vec.end(); }
        auto end() const noexcept { return vec.end(); }
        auto empty() const noexcept { return vec.empty(); }
        auto size() const noexcept { return vec.size(); }
        auto& back() noexcept { return vec.back(); }
        auto& back() const noexcept { return vec.back(); }
        auto& front() noexcept { return vec.front(); }
        auto& front() const noexcept { return vec.front(); }
        auto& at(size_t n) { return vec.at(n); }
        auto& at(size_t n) const { return vec.at(n); }

        template <typename Tderived, class... Args, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
        void emplace_back(Args &... args) {
            vec.push_back(std::make_shared<Tderived>(args...));
        } //!< Create an (derived) instance and append a pointer to it to the vector

        template <typename Tderived, class... Args, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
        void emplace_back(const Args &... args) {
            vec.push_back(std::make_shared<Tderived>(args...));
        } //!< Create an (derived) instance and append a pointer to it to the vector

        template <typename Tderived, class Arg, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
        void push_back(std::shared_ptr<Arg> arg) {
            vec.push_back(arg);
        } //!< Append a pointer to a (derived) instance to the vector

        template <typename Tderived, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
        auto& operator=(const BasePointerVector<Tderived> &d) {
            vec.assign(d.vec.begin(), d.vec.end());
            return *this;
        } //!< Allow assignment to a vector of ancestors

        template <typename Tderived, class = std::enable_if_t<std::is_base_of<T, Tderived>::value>>
        auto find() const {
            BasePointerVector<Tderived> _v;
            for (auto base : vec) {
                auto derived = std::dynamic_pointer_cast<Tderived>(base);
                if (derived)
                    _v.template push_back<Tderived>(derived);
            }
            return _v;
        } //!< Pointer list to all matching type

        template <typename U>
        friend void to_json(json &, const BasePointerVector<U> &); //!< Allow serialization to JSON
    }; //!< Helper class for storing vectors of base pointers

    template <typename T> void to_json(json &j, const BasePointerVector<T> &b) {
        try {
            for (auto shared_ptr : b.vec) {
                j.push_back(*shared_ptr);
            }
        } catch (const std::exception &e) {
            throw std::runtime_error("error converting to json: "s + e.what());
        }
    }

    template <typename T> void from_json(const json &j, BasePointerVector<T> &b) {
        try {
            for (auto it : j) {
                std::shared_ptr<T> ptr = it;
                b.template push_back<T>(ptr);
            }
        } catch (const std::exception &e) {
            throw std::runtime_error("error converting from json: "s + e.what());
        }
    }
} //namespace Faunus
