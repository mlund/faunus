#ifndef FAUNUS_AUXILIARY_H
#define FAUNUS_AUXILIARY_H

#ifndef SWIG
#include <utility>
#ifdef FAU_HASHTABLE
#include <unordered_map>
#include <functional>
#endif
#endif

namespace Faunus {

  /**
   * @brief Iterate over pairs in container and return sum of function
   */
  template<typename Tit, typename Tfunc, typename T=double>
    T for_each_pair_sum(const Tit &begin, const Tit &end, Tfunc f, T sum=T())
    {
      for (auto i=begin; i!=end; ++i)
        for (auto j=i; ++j!=end;)
          sum+=f(*i,*j);
      return sum;
    }

  /**
   * @brief Ordered pair where `first<=second`
   *
   * Upon construction, the smallest element is placed in `first`
   * so that `opair<int>(i,j)==opair<int>(j,i)` is always true.
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
   * @note see Cawley 2000; doi:10.1162/089976600300015033
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

  /** @brief Convert string to float, int, bool */
  template<class T>
    T str2val(const std::string &s, T fallback=T()) {
      return (!s.empty()) ? T(std::stod(s)) : fallback;
    }
  template<>
    inline bool str2val<bool>(const std::string &s, bool fallback) {
      if (s=="yes" || s=="true") return true;
      if (s=="no" || s=="false") return false;
      return fallback;
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
}//namespace
#endif

