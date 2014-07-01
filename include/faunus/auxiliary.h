#ifndef FAUNUS_AUXILIARY_H
#define FAUNUS_AUXILIARY_H

#ifndef SWIG
#include <utility>
#include <algorithm>
#include <utility>
#include <regex>

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
namespace Faunus {

  /**
   * @brief Iterate over pairs in container, call a function on the elements, and sum results
   */
  template<typename Tit, typename Tfunc, typename T=double, typename Top>
    T for_each_pair(const Tit &begin, const Tit &end, Tfunc f, Top operation=std::plus<T>())
    {
      T x;
      for (auto i=begin; i!=end; ++i)
        for (auto j=i; ++j!=end;)
          x = operation(x, f(*i,*j));
      return x;
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
    struct opair : public std::pair<T,T> {
      typedef std::pair<T,T> base;
      opair() {}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-overflow"
      opair(const T &a, const T &b) : base(a,b) {
        if (a>b)
          std::swap(base::first,base::second);
      }
#pragma GCC diagnostic pop
      bool find(const T &i) const {
        assert(base::first<=base::second);
        if (i!=base::first)
          if (i!=base::second)
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

namespace Faunus {
  /**
   * @brief Store data for pairs
   */
  template<typename Tdata, typename T=int>
    class pair_list {
      protected:
        typedef opair<T> Tpair; // ordered pair
        std::map<Tpair,Tdata> list;   // main pair list
        std::multimap<T,T> mlist; // additional map for faster access
      public:
        /** @brief Associate data with a pair */
        void add(T i, T j, Tdata d) {
          list[ Tpair(i,j) ] = d; 
          mlist.insert( std::pair<T,T>(i,j) );
          mlist.insert( std::pair<T,T>(j,i) );
        }

        /** @brief Access data of a pair */
        Tdata& operator() (T i, T j) {
          Tpair ij(i,j);
          assert( list[ij] != nullptr ); //debug
          return list[ij];
        }

        /** @brief Clears all data */
        void clear() {
          list.clear();
          mlist.clear();
        }

        decltype(list)& getBondList() { return list; }
    };

  template<class Tdata, class T=int, class Tbase=std::map<opair<T>,Tdata> >
    struct map_ij : private Tbase {
      using Tbase::begin;
      using Tbase::end;
      Tdata& operator() (T i, T j)  {
        return Tbase::operator[]( opair<T>(i,j) );
      }
      Tdata& operator() (T i, T j) const {
        return Tbase::operator[]( opair<T>(i,j) );
      }
      typename Tbase::const_iterator find(T i, T j) const {
        return Tbase::find(opair<T>(i,j));
      }
      typename Tbase::iterator find(T i, T j) {
        return Tbase::find(opair<T>(i,j));
      }
    };

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
  /**
   * @brief Quake inverse square root approximation
   */
  template<class Tint=int>
    float invsqrtQuake(float number) {
      static_assert(sizeof(Tint)==4,
          "Integer size must be 4 bytes for quake invsqrt.");
      float y  = number;
      float x2 = y * 0.5F;
      Tint i  = * ( Tint * ) &y;
      i  = 0x5f3759df - ( i >> 1 );
      y  = * ( float * ) &i;
      y  = y * ( 1.5F - ( x2 * y * y ) );   // 1st iteration
      //      y  = y * ( 1.5F - ( x2 * y * y ) );   // 2nd iteration, this can be removed
      return y;
    }

  /**
   * @brief Approximate exp() function
   * @note see [Cawley 2000](http://dx.doi.org/10.1162/089976600300015033)
   * @warning Does not work in big endian systems!
   */
  template<class Tint=int>
    double exp_cawley(double y) {
      static_assert(2*sizeof(Tint)==sizeof(double),
          "Approximate exp() requires 4-byte integer" );
      //static double EXPA=1048576/std::log(2);
      union {
        double d;
        struct { Tint j, i; } n;  // little endian
        //struct { int i, j; } n;  // bin endian
      } eco;
      eco.n.i = 1072632447 + (Tint)(y*1512775.39519519);
      eco.n.j = 0;
      return eco.d;
    }

  template<class Tint>
    double exp_untested(double y) {
      static_assert(2*sizeof(Tint)==sizeof(double),
          "Approximate exp() requires 4-byte integer");
      double d;
      *((Tint*)(&d) + 0) = 0;
      *((Tint*)(&d) + 1) = (Tint)(1512775 * y + 1072632447);
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
    T str2val(const std::string &s, T fallback=T()) {
      return (!s.empty()) ? T(std::stod(s)) : fallback;
    }
  template<>
    inline bool str2val<bool>(const std::string &s, bool fallback) {
      if (std::regex_match(s, std::regex("yes|true|on", std::regex_constants::icase) ) )
        return true;
      return fallback;
    }

  /**
   * @brief Fast sine calculation in range (-pi:pi)
   * @warning Invalid beyond boundaries!
   * @note <http://devmaster.net/forums/topic/4648-fast-and-accurate-sinecosine>
   */
  template<class T>
    T sinApprox(T x) {
      const T B = 4/std::acos(-1);
      const T C = -B/std::acos(-1);

      T y = B * x + C * x * abs(x);

#ifdef EXTRA_PRECISION
      //  const float Q = 0.775;
      const T P = 0.225;

      y = P * (y * abs(y) - y) + y;   // Q * y + P * y * abs(y)
#endif
      return y;
    }

  template<class T>
    T cosApprox(T x) {
      const T shift = 0.5*std::acos(-1);
      return sinApprox(x+shift);
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
    class legendre {
      private:
        int n;           //!< Legendre order
        void resize(int order) {
          n=order;
          assert(n>=0);
          P.resize(n+1);
          P[0]=1.;
        }
      public:
        /** @brief Construct w. polynomium order>=0 */
        legendre(int order) { resize(order); }

        /** @brief Legendre terms stored here */
        std::vector<T> P;

        /** @brief Evaluate polynomium at x */
        void eval(T x) {
          if (n>0) {
            P[1] = x;
            for (int i=1; i<n; ++i)
              P[i+1] = (2*i+1)/(i+1)*x*P[i] - i/(i+1)*P[i-1];
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
    Tint to_bin(T x, T dx=1) {
      return (x<0) ? Tint( x/dx-0.5 ) : Tint( x/dx+0.5 );
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
    class InterpolTable {
      private:
        typedef std::pair<T,T> Tpair; // xy data stored as pairs
        std::vector<Tpair> t;         // in a vector
      public:
        InterpolTable(const std::string &filename) {
          Tpair a;
          std::ifstream in(filename);
          while (in >> a.first >> a.second)
            t.push_back(a);
          std::sort(t.begin(), t.end()); // sort
          t.erase( std::unique( t.begin(), t.end() ), t.end() ); // remove duplicates 
        }

        T operator()(T x) const {
          assert(!t.empty() && "Table is empty");
          if (x>t.back().first) return std::numeric_limits<T>::quiet_NaN();
          if (x<t[0].first) return std::numeric_limits<T>::quiet_NaN();
          auto it = std::lower_bound(t.begin(), t.end(), Tpair(x,0));
          if (it==t.begin())
            return it->second;
          auto it2=it;
          --it2;
          return it2->second + (it->second-it2->second)*(x-it2->first)/(it->first-it2->first);
        }
    }; // end of InterpolTable

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
    class PairMatrix {
      public:
        vector< vector<T> > m; // symmetric matrix (mem.wasteful - fast access)
        void resize(size_t n) {
          m.resize(n);
          for (auto &i : m)
            i.resize(n,0);
        }
        PairMatrix(size_t n=0) {
          resize(n);
        }
        const T& operator()(size_t i, size_t j) const {
          assert( i<m.size() );
          assert( j<m[i].size() );
          assert( m[i][j]==m[j][i] );
          return m[i][j]; 
        }
        void set(size_t i, size_t j, T val) {
          size_t n=std::max(i,j);
          if (n>=m.size())
            resize(n+1);
          m[i][j]=m[j][i]=val;
        }
    };

}//namespace
#endif

