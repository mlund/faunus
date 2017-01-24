#ifndef FAUNUS_AUXILIARY_H
#define FAUNUS_AUXILIARY_H

#ifndef SWIG
#include <utility>
#include <algorithm>
#include <utility>
#include <regex>
#include <cstdint>
#include <chrono>

#ifdef FAU_HASHTABLE
#include <unordered_map>
#include <functional>
#endif
#endif

/**
 * @file auxiliary.h
 *
 * This file contains auxiliary functionality that
 * have no dependencies other than STL and can hence
 * be copied to other projects.
 */
namespace Faunus
{

  /**
   * @brief Iterate over pairs in container, call a function on the elements, and sum results
   */
  template<typename Tit, typename Tfunc, typename T=double, typename Top>
  T for_each_pair( const Tit &begin, const Tit &end, Tfunc f, Top operation = std::plus<T>())
  {
      T x;
      for ( auto i = begin; i != end; ++i )
          for ( auto j = i; ++j != end; )
              x = operation(x, f(*i, *j));
      return x;
  }

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
}//namespace

/*
 * The following is an extension to `std::hash` in order to
 * construct hash tables (`std::unordered_map`), using `opair<T>` as keys.
 * See more at:
 * <http://en.wikipedia.org/wiki/Unordered_associative_containers_(C%2B%2B)#Usage_example>
 */
#ifdef FAU_HASHTABLE
namespace std {
  template<typename T>
    struct hash< Faunus::opair<T> > {
      size_t operator()(const Faunus::opair<T> &x) const {
        return hash<T>()(x.first) ^ hash<T>()(x.second);
      }
    };
}//namespace
#endif

namespace Faunus
{
  /**
   * @brief Store data for pairs
   */
  template<typename Tdata, typename T=int>
  class pair_list
  {
  protected:
      typedef opair<T> Tpair; // ordered pair
      std::map<Tpair, Tdata> list;   // main pair list
      std::multimap<T, T> mlist; // additional map for faster access
  public:
      /** @brief Associate data with a pair */
      void add( T i, T j, Tdata d )
      {
          list[Tpair(i, j)] = d;
          mlist.insert(std::pair<T, T>(i, j));
          mlist.insert(std::pair<T, T>(j, i));
      }

      /** @brief Access data of a pair */
      Tdata &operator()( T i, T j )
      {
          Tpair ij(i, j);
          assert(list[ij] != nullptr); //debug
          return list[ij];
      }

      /** @brief Clears all data */
      void clear()
      {
          list.clear();
          mlist.clear();
      }

      decltype(list) &getBondList() { return list; }
  };

  template<class Tdata, class T=int, class Tbase=std::map<opair<T>, Tdata> >
  struct map_ij : private Tbase
  {
      using Tbase::begin;
      using Tbase::end;

      Tdata &operator()( T i, T j )
      {
          return Tbase::operator[](opair<T>(i, j));
      }

      Tdata &operator()( T i, T j ) const
      {
          return Tbase::operator[](opair<T>(i, j));
      }

      typename Tbase::const_iterator find( T i, T j ) const
      {
          return Tbase::find(opair<T>(i, j));
      }

      typename Tbase::iterator find( T i, T j )
      {
          return Tbase::find(opair<T>(i, j));
      }
  };

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
      return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))) * exp(-x * x);
  }

  /**
   * @brief Approximate 1 - erfc_x
   * @param x Value for which erf should be calculated 
   */
  template<typename T>
  T erf_x( T x ) { return (1 - erfc_x(x)); }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"

  /**
   * @brief Quake inverse square root approximation
   */
  template<class Tint=std::int32_t>
  float invsqrtQuake( float number )
  {
      static_assert(sizeof(Tint) == 4, "Integer size must be 4 bytes for quake invsqrt.");
      float y = number;
      float x2 = y * 0.5F;
      Tint i = *(Tint *) &y;
      i = 0x5f3759df - (i >> 1);
      y = *(float *) &i;
      y = y * (1.5F - (x2 * y * y));   // 1st iteration
      //      y  = y * ( 1.5F - ( x2 * y * y ) );   // 2nd iteration, this can be removed
      return y;
  }

  /**
   * @brief n'th integer power of float
   *
   * On GCC/Clang this will use the fast `__builtin_powi` function.
   * If not, and `n<7`, a simple loop (that can be unrolled at compile
   * time) is performed. If none of the above, `std::pow` is used.
   */
  template<int n, typename T=double>
  T _powi( T &x )
  {
#if defined(__GNUG__)
      return __builtin_powi(x, n);
#else
      if (n>6)
        return std::pow(x,n);
      while (--n>0)
        x*=x;
      return x;
#endif
  }

  /**
   * @brief Approximate exp() function
   * @note see [Cawley 2000](http://dx.doi.org/10.1162/089976600300015033)
   * @warning Does not work in big endian systems!
   */
  template<class Tint=std::int32_t>
  double exp_cawley( double y )
  {
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

  template<class Tint=std::int32_t>
  double exp_untested( double y )
  {
      static_assert(2 * sizeof(Tint) == sizeof(double),
                    "Approximate exp() requires 4-byte integer");
      double d;
      *((Tint *) (&d) + 0) = 0;
      *((Tint *) (&d) + 1) = (Tint) (1512775 * y + 1072632447);
      return d;
  }

#pragma GCC diagnostic pop

  /**
   * @brief Convert string to float, int, bool
   *
   * Examples:
   *
   *     double x = str2val("10.0");       // -> 10.0
   *     double y = str2val("",0.5);       // -> 0.5 (fallback)
   *     bool z   = str2val<bool>("true"); // -> true
   *
   * Boolean text can be "yes/true/on" - if not, the default
   * fallback, `false` is returned. Matching is case insensitive.
   */
  template<class T>
  T str2val( const std::string &s, T fallback = T())
  {
      return (!s.empty()) ? T(std::stod(s)) : fallback;
  }

  template<>
  inline bool str2val<bool>( const std::string &s, bool fallback )
  {
      if ( std::regex_match(s, std::regex("yes|true|on|aye", std::regex_constants::icase)))
          return true;
      return fallback;
  }

  /**
   * @brief Fast sine calculation in range (-pi:pi)
   * @warning Invalid beyond boundaries!
   * @note <http://devmaster.net/forums/topic/4648-fast-and-accurate-sinecosine>
   */
  template<class T>
  T sinApprox( T x )
  {
      const T B = 4 / std::acos(-1);
      const T C = -B / std::acos(-1);

      T y = B * x + C * x * abs(x);

#ifdef EXTRA_PRECISION
      //  const float Q = 0.775;
      const T P = 0.225;

      y = P * (y * abs(y) - y) + y;   // Q * y + P * y * abs(y)
#endif
      return y;
  }

  template<class T>
  T cosApprox( T x )
  {
      const T shift = 0.5 * std::acos(-1);
      return sinApprox(x + shift);
  }

  /**
   * @brief Evaluate n'th degree Legendre polynomium
   *
   * Example:
   * @code
   * legendre<float> l(10);
   * l.eval(1.3);
   * cout << l.p[3]
   * @endcode
   *
   * @author Mikael Lund
   * @date Canberra 2005-2006
   */
  template<typename T=double>
  class legendre
  {
  private:
      int n;           //!< Legendre order
      void resize( int order )
      {
          n = order;
          assert(n >= 0);
          P.resize(n + 1);
          P[0] = 1.;
      }

  public:
      /** @brief Construct w. polynomium order>=0 */
      legendre( int order ) { resize(order); }

      /** @brief Legendre terms stored here */
      std::vector<T> P;

      /** @brief Evaluate polynomium at x */
      void eval( T x )
      {
          if ( n > 0 )
          {
              P[1] = x;
              for ( int i = 1; i < n; ++i )
                  P[i + 1] = (2 * i + 1) / (i + 1) * x * P[i] - i / (i + 1) * P[i - 1];
          }
      }
  };

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
  Tint to_bin( T x, T dx = 1 )
  {
      return (x < 0) ? Tint(x / dx - 0.5) : Tint(x / dx + 0.5);
  }

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
   * This will maintain a symmetric, dynamic NxN matrix for storing data
   * about pairs.
   * Use the `set()` function for setting values and the function
   * operator for access:
   *
   *     int i=2,j=3; // particle type, for example
   *     PairMatrix<double> m;
   *     m.set(i,j,12.0);
   *     cout << m(i,j);         // -> 12.0
   *     cout << m(i,j)==m(j,i); // -> true
   */
  template<class T=double>
  class PairMatrix
  {
  public:
      vector <vector<T>> m; // symmetric matrix (mem.wasteful - fast access)
      void resize( size_t n )
      {
          m.resize(n);
          for ( auto &i : m )
              i.resize(n, T());
      }

      PairMatrix( size_t n = 0 ) { resize(n); }

      size_t size() const { return m.size(); }

      const T &operator()( size_t i, size_t j ) const
      {
          assert(i < m.size());
          assert(j < m[i].size());
          //assert( m[i][j]==m[j][i] );
          return m[i][j];
      }

      void set( size_t i, size_t j, T val )
      {
          size_t n = std::max(i, j);
          if ( n >= m.size())
              resize(n + 1);
          m[i][j] = m[j][i] = val;
      }
  };

  template<typename Tcoeff=double, typename base=Eigen::Matrix<Tcoeff, Eigen::Dynamic, Eigen::Dynamic>>
  class Table : public base
  {
  private:
      typedef std::vector<double> Tvec;
      Tvec _bw, _lo, _hi;
      int _rows, _cols;
  public:
      Table( const Tvec &bw = {1, 1}, const Tvec &lo = {0, 0}, const Tvec &hi = {2, 2} )
      {
          _bw = bw;
          _lo = lo;
          _hi = hi;
          _rows = (_hi[0] - _lo[0]) / _bw[0] + 1.;
          _cols = (_hi[1] - _lo[1]) / _bw[1] + 1.;
          base::resize(_rows, _cols);
          base::setZero();
      }

      void reInitializer( Tvec &bw, Tvec &lo, Tvec &hi )
      {
          _bw = bw;
          _lo = lo;
          _hi = hi;
          _rows = (_hi[0] - _lo[0]) / _bw[0] + 1.;
          _cols = (_hi[1] - _lo[1]) / _bw[1] + 1.;
          base::resize(_rows, _cols);
          base::setZero();
      }

      void clear() { base::setZero(); }

      void round( Tvec &v )
      {
          for ( Tvec::size_type i = 0; i != v.size(); ++i )
              v[i] = (v[i] >= 0) ? int(v[i] / _bw[i] + 0.5) * _bw[i] : int(v[i] / _bw[i] - 0.5) * _bw[i];
      }

      void to_index( Tvec &v )
      {
          for ( Tvec::size_type i = 0; i != v.size(); ++i )
          {
              v[i] = (v[i] >= 0) ? int(v[i] / _bw[i] + 0.5) : int(v[i] / _bw[i] - 0.5);
              v[i] = v[i] - _lo[i] / _bw[i];
          }
          v.resize(2, 0);
      }

      Tcoeff &operator[]( const Tvec &v )
      {
          return base::operator()(v[0], v[1]);
      }

      bool isInRange( const Tvec &v )
      {
          bool b = true;
          for ( Tvec::size_type i = 0; i != v.size(); ++i )
              b = b && v[i] >= _lo[i] && v[i] <= _hi[i];
          return b;
      }

      Tvec hist2buf( int &size )
      {
          Tvec sendBuf;
          for ( int i = 0; i < _cols; ++i )
              for ( int j = 0; j < _rows; ++j )
                  sendBuf.push_back(base::operator()(j, i));
          return sendBuf;
      }

      void buf2hist( Tvec &v )
      {
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

      base getBlock( const Tvec &v )
      { // {xmin,xmax} or {xmin,xmax,ymin,ymax}
          Tvec w = {0, 0, 0, 0};
          switch (v.size())
          {
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

      Tcoeff avg( const Tvec &v )
      {
          return this->getBlock(v).mean();
      }

      void translate( Tcoeff s )
      {
          *this += base::Constant(_rows, _cols, s);
      }

      void save( const string &filename, Tcoeff scale = 1, Tcoeff translate = 0 )
      {
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

      void saveRow( const string &filename, const Tvec &v, Tcoeff scale = 1, Tcoeff translate = 0 )
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

      void load( const string &filename )
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
      string name;
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
      std::map<string, vector < double>> to_map() {
          std::map<string, vector < double>>
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
      void save( const string &filename, T scale = 1, T translate = 0 )
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
      void normSave( const string &filename )
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
      void sumSave( string filename, T scale = 1 )
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
      vector<double> hist2buf( int &size )
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
      void buf2hist( vector<double> &v )
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
      bool load( const string &filename )
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

  template<typename Tx, typename Ty>
  class Table3D
  {
  protected:
      typedef std::map<std::pair<Tx, Tx>, Ty> Tmap;
      Tx dx1, dx2;
      Tmap map;
      string name;

      Ty count()
      {
          Ty cnt = 0;
          for ( auto &m : map )
              cnt += m.second;
          return cnt;
      }

  private:
      Tx round1( Tx x ) { return (x >= 0) ? int(x / dx1 + 0.5) * dx1 : int(x / dx1 - 0.5) * dx1; }

      Tx round2( Tx x ) { return (x >= 0) ? int(x / dx2 + 0.5) * dx2 : int(x / dx2 - 0.5) * dx2; }

      virtual double get( Tx x1, Tx x2 ) { return operator()(x1, x2); }

  public:
      enum type { HISTOGRAM, XYDATA };
      type tabletype;

      /**
       * @brief Constructor
       * @param resolution Resolution of the x axis
       * @param key Table type: HISTOGRAM or XYDATA
       */
      Table3D( Tx resolution1 = 1, Tx resolution2 = 1, type key = XYDATA )
      {
          tabletype = key;
          setResolution(resolution1, resolution2);
      }

      void clear() { map.clear(); }

      void setResolution( Tx resolution1, Tx resolution2 )
      {
          assert(resolution1 > 0 && resolution2 > 0);
          dx1 = resolution1;
          dx2 = resolution2;
          map.clear();
      }

      void setResolution( std::vector<Tx> &resolution )
      {
          assert(resolution.at(0) > 0 && resolution.at(1) > 0);
          dx1 = resolution[0];
          dx2 = resolution[1];
          map.clear();
      }

      virtual ~Table3D() {}

      /** @brief Access operator - returns reference to y(x) */
      Ty &operator()( Tx x1, Tx x2 )
      {
          return map[std::make_pair(round1(x1), round2(x2))];
      }

      Ty &operator()( std::vector<Tx> &x )
      {
          return map[std::make_pair(round1(x[0]), round2(x[1]))];
      }

      /** @brief Find key and return corresponding value otherwise zero*/
      Ty find( std::vector<Tx> &x )
      {
          Ty value = 0;
          auto it = map.find(std::make_pair(round1(x[0]), round2(x[1])));
          if ( it != map.end())
              value = it->second;
          return value;
      }

      /** @brief Save table to disk */
      void save( string filename, Ty scale = 1, Ty translate = 0 )
      {
          if ( tabletype == HISTOGRAM )
          { // compensate for half bin width
              auto first = map.begin();
              auto last = map.rbegin();
              if ( !map.empty())
              {
                  for ( auto it = first; it != map.end(); ++it )
                  {
                      if ( it->first.first == first->first.first )
                          it->second *= 2;
                      else if ( it->first.second == first->first.second )
                          it->second *= 2;
                  }
              }
              if ( map.size() > 1 )
              {
                  for ( auto it = last; it != map.rend(); ++it )
                  {
                      if ( it->first.first == last->first.first )
                          it->second *= 2;
                      else if ( it->first.second == last->first.second )
                          it->second *= 2;
                  }
              }
          }

          if ( !map.empty())
          {
              std::ofstream f(filename.c_str());
              f.precision(10);
              if ( f )
              {
                  for ( auto &m : map )
                      f << m.first.first << " " << m.first.second
                        << " " << (m.second + translate) * scale << "\n";
              }
          }

          if ( tabletype == HISTOGRAM )
          { // restore half bin width
              auto first = map.begin();
              auto last = map.rbegin();
              if ( !map.empty())
              {
                  for ( auto it = first; it != map.end(); ++it )
                  {
                      if ( it->first.first == first->first.first )
                          it->second /= 2;
                      else if ( it->first.second == first->first.second )
                          it->second /= 2;
                  }
              }
              if ( map.size() > 1 )
              {
                  for ( auto it = last; it != map.rend(); ++it )
                  {
                      if ( it->first.first == last->first.first )
                          it->second /= 2;
                      else if ( it->first.second == last->first.second )
                          it->second /= 2;
                  }
              }
          }
      }

      /** @brief Save normalized table to disk */
      void normSave( string filename )
      {
          if ( tabletype == HISTOGRAM )
          { // compensate for half bin width
              auto first = map.begin();
              auto last = map.rbegin();
              if ( !map.empty())
              {
                  for ( auto it = first; it != map.end(); ++it )
                  {
                      if ( it->first.first == first->first.first )
                          it->second *= 2;
                      else if ( it->first.second == first->first.second )
                          it->second *= 2;
                  }
              }
              if ( map.size() > 1 )
              {
                  for ( auto it = last; it != map.rend(); ++it )
                  {
                      if ( it->first.first == last->first.first )
                          it->second *= 2;
                      else if ( it->first.second == last->first.second )
                          it->second *= 2;
                  }
              }
          }

          if ( !map.empty())
          {
              std::ofstream f(filename.c_str());
              f.precision(10);
              Ty cnt = count() * dx1 * dx2;
              if ( f )
              {
                  for ( auto &m : map )
                      f << m.first.first << " " << m.first.second
                        << " " << m.second / cnt << "\n";
              }
          }

          if ( tabletype == HISTOGRAM )
          { // restore half bin width
              auto first = map.begin();
              auto last = map.rbegin();
              if ( !map.empty())
              {
                  for ( auto it = first; it != map.end(); ++it )
                  {
                      if ( it->first.first == first->first.first )
                          it->second /= 2;
                      else if ( it->first.second == first->first.second )
                          it->second /= 2;
                  }
              }
              if ( map.size() > 1 )
              {
                  for ( auto it = last; it != map.rend(); ++it )
                  {
                      if ( it->first.first == last->first.first )
                          it->second /= 2;
                      else if ( it->first.second == last->first.second )
                          it->second /= 2;
                  }
              }
          }
      }

      Tmap getMap()
      {
          return map;
      }

      std::pair<Tx, Tx> getResolution()
      {
          return std::make_pair(dx1, dx2);
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

      /*! Returns average in interval */
      Ty avg( const std::vector<Tx> &limits )
      {
          Ty avg = 0;
          int cnt = 0;
          assert(!map.empty());
          for ( auto &m : map )
          {
              if ( m.first.first >= limits[0] && m.first.first <= limits[1]
                  && m.first.second >= limits[2] && m.first.second <= limits[3] )
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
       * @brief Convert table3D to vector of floats
       */
      vector<double> hist2buf( int &size )
      {
          std::vector<double> sendBuf;
          assert(!map.empty());
          for ( auto &m : map )
          {
              sendBuf.push_back(m.first.first);
              sendBuf.push_back(m.first.second);
              sendBuf.push_back(m.second);
          }
          sendBuf.resize(size, -1);
          return sendBuf;
      }

      /**
       * @brief Convert vector of floats to table3D
       */
      void buf2hist( vector<double> &v )
      {
          this->clear();
          assert(!v.empty());
          std::map<std::pair<double, double>, Average<double>> all;
          for ( int i = 0; i < int(v.size()) - 2; i += 3 )
              if ( v.at(i + 2) != -1 )
                  all[std::make_pair(v.at(i), v.at(i + 1))] += v.at(i + 2);
          for ( auto &m : all )
              this->operator()(m.first.first, m.first.second) = m.second.avg();
      }

      /**
       * @brief Load table from disk
       */
      bool load( const string &filename )
      {
          std::ifstream f(filename.c_str());
          if ( f )
          {
              map.clear();
              while ( !f.eof())
              {
                  Tx x1, x2;
                  Ty y;
                  f >> x1 >> x2 >> y;
                  operator()(x1, x2) = y;
              }
              if ( tabletype == HISTOGRAM )
              { // restore half bin width
                  if ( !map.empty())
                  {
                      auto first = map.begin();
                      auto last = map.rbegin();
                      for ( auto it = first; it != map.end(); ++it )
                      {
                          if ( it->first.first == first->first.first )
                              it->second /= 2;
                          else if ( it->first.second == first->first.second )
                              it->second /= 2;
                      }
                      if ( map.size() > 1 )
                      {
                          for ( auto it = last; it != map.rend(); ++it )
                          {
                              if ( it->first.first == last->first.first )
                                  it->second /= 2;
                              else if ( it->first.second == last->first.second )
                                  it->second /= 2;
                          }
                      }
                  }
              }
              return true;
          }
          return false;
      }
  };

  /**
   * @brief Subtract two tables
   */
  template<class Tx, class Ty>
  Table3D<Tx, Ty> operator-( Table3D<Tx, Ty> &a, Table3D<Tx, Ty> &b )
  {
      assert(a.tabletype == b.tabletype && "Table a and b needs to be of same type");
      Table3D<Tx, Ty> c(std::min(a.getResolution().first, b.getResolution().first),
                        std::min(a.getResolution().second, b.getResolution().second), a.tabletype);
      auto a_map = a.getMap();
      auto b_map = b.getMap();

      if ( a.tabletype == Table3D<Tx, Ty>::HISTOGRAM )
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
              c(m1.first.first, m1.first.second) = m1.second - m2.second;
              break;
          }
      }

      if ( a.tabletype == Table3D<Tx, Ty>::HISTOGRAM )
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
  Table3D<Tx, Ty> operator+( Table3D<Tx, Ty> &a, Table3D<Tx, Ty> &b )
  {
      assert(a.tabletype == b.tabletype && "Table a and b needs to be of same type");
      Table3D<Tx, Ty> c(std::min(a.getResolution().first, b.getResolution().first),
                        std::min(a.getResolution().second, b.getResolution().second), a.tabletype);
      auto a_map = a.getMap();
      auto b_map = b.getMap();

      if ( a.tabletype == Table3D<Tx, Ty>::HISTOGRAM )
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
          c(m.first.first, m.first.second) += m.second;
      }
      for ( auto &m : b_map )
      {
          c(m.first.first, m.first.second) += m.second;
      }

      if ( a.tabletype == Table3D<Tx, Ty>::HISTOGRAM )
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

  template<typename Tx, typename Ty=unsigned long int>
  class Histogram : public Table2D<Tx, Ty>
  {
  public:
      Histogram( Tx resolution = 0.2 ) : Table2D<Tx, Ty>(resolution, Table2D<Tx, Ty>::HISTOGRAM)
      {
          static_assert(std::is_integral<Ty>::value, "Histogram must be of integral type");
          static_assert(std::is_unsigned<Ty>::value, "Histogram must be unsigned");
      }
  };

  /**
   * @brief Finds pointer to element in tuple with specified type. `nullptr` if not found.
   *
   * Example:
   *
   *     std::tuple<int,float,bool> t( 1, 2.1, false );
   *     std::cout << *TupleFindType::get<float>( t ); // -> 2.1
   */
  class TupleFindType
  {
  private:
      template<class T>
      struct findtype
      {
          T *ptr;

          findtype() : ptr(nullptr) {}

          template<class E>
          void operator()( E &t, typename std::enable_if<std::is_same<T, E>::value>::type * = 0 ) { ptr = &t; }

          template<class E>
          void operator()( E &t, typename std::enable_if<!std::is_same<T, E>::value>::type * = 0 ) {}
      };

      template<class Tuple, std::size_t N>
      struct TuplePrinter
      {
          template<class Tfunc>
          static void print( Tuple &t, Tfunc &f )
          {
              TuplePrinter<Tuple, N - 1>::print(t, f);
              f(std::get<N - 1>(t));
          }
      };

      template<class Tuple>
      struct TuplePrinter<Tuple, 1>
      {
          template<class Tfunc>
          static void print( Tuple &t, Tfunc &f ) { f(std::get<0>(t)); }
      };

      template<class... Args, class Tfunc>
      static void for_each( std::tuple<Args...> &t, Tfunc &f )
      {
          TuplePrinter<decltype(t), sizeof...(Args)>::print(t, f);
      }

  public:
      template<class T, class... Args>
      static T *get( std::tuple<Args...> &t )
      {
          findtype<T> func;   // apply function object `func` ...
          for_each(t, func);  // ... on all elements
          return func.ptr;
      }
  };

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
      TimeRelativeOfTotal() : delta(0)
      {
          t0 = std::chrono::steady_clock::now();
      }

      void start() { tx = std::chrono::steady_clock::now(); }

      void stop()
      {
          delta += std::chrono::duration_cast<Tunit>
              (std::chrono::steady_clock::now() - tx);
      }

      double result() const
      {
          auto now = std::chrono::steady_clock::now();
          auto total = std::chrono::duration_cast<Tunit>(now - t0);
          return delta.count() / double(total.count());
      }
  };

}//namespace
#endif
