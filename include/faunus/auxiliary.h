#ifndef FAU_auxiliary
#define FAU_auxiliary

#ifndef SWIG
#include <utility>
#ifdef FAU_HASHTABLE
#include <unordered_map>
#include <functional>
#endif
#endif

namespace Faunus {
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
   * This is a list of pairs with associated data where the latter should de derived
   * from the base Tbase. When adding data with the add() function, a copy of the data
   * is created and stored internally.
   */
  template<class Tbase, typename Tij=int, typename Tpair=opair< Tij > >
    class pair_list {
      protected:
#ifdef FAU_HASHTABLE
        typedef std::unordered_map< Tpair, std::shared_ptr<Tbase> > Tlist;
#else
        typedef std::map< Tpair, std::shared_ptr<Tbase> > Tlist;
#endif
        Tlist list;
      public:
        /**
         * @brief Associate data with a pair using an internal copy.
         *
         * Data is added by making an internal COPY of the given `Tderived` object.
         * For large lists, consider adding a pointer instead using the alternative
         * `add()` function.
         */
        template<typename Tderived>
          void add(Tij i, Tij j, Tderived data) {
            list[ Tpair(i,j) ] = std::shared_ptr<Tderived>( new Tderived(data) ); 
          }

        /**
         * @brief Associate data with a pair using pointers.
         */
        template<typename Tderived>
          void add(Tij i, Tij j, std::shared_ptr<Tderived>& sptr) {
            list[ Tpair(i,j) ] = sptr; 
          }

        /**
         * @brief Access data of a pair
         */
        Tbase& operator() (Tij i, Tij j) {
          Tpair pair(i,j);
          assert( list[pair] != nullptr ); //debug
          return *list[pair];
        }

        /** @brief Clears all data */
        void clear() { list.clear(); }
    };

  template<typename Tij, typename Tval>
    struct map_ij {
      typedef opair<Tij> Tkey;
      std::map<Tkey,Tval> list;
      Tval& operator() (Tij i, Tij j) {
        return list[ Tpair(i,j) ];
      }
      //std::map<Tkey,Tval>::iterator
    };

  template<typename Tparticle, typename Tij=int, typename Tpair=opair< Tij > >
    class pair_list_functor {
      protected:
        typedef std::function<double(const Tparticle&, const Tparticle&, double)> Tfunctor;
#ifdef FAU_HASHTABLE
        typedef std::unordered_map< Tpair, Tfunctor > Tlist;
#else
        typedef std::map< Tpair, Tfunctor > Tlist;
#endif
        Tlist list;
      public:
        /**
         * @brief Associate data with a pair using an internal copy.
         *
         * Data is added by making an internal COPY of the given `Tderived` object.
         * For large lists, consider adding a pointer instead using the alternative
         * `add()` function.
         */
        void add(Tij i, Tij j, Tfunctor f) {
          list[ Tpair(i,j) ] = f; 
        }

        /**
         * @brief Access data of a pair
         */
        Tfunctor& operator() (Tij i, Tij j) {
          Tpair pair(i,j);
          assert( list[pair] != nullptr ); //debug
          return list[pair];
        }

        /** @brief Clears all data */
        void clear() { list.clear(); }
    };

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
  /**
   * @brief Quake inverse square root approximation
   */
  inline float invsqrtQuake(float number) {
    assert(sizeof(int)==4 && "Integer size must be 4 bytes for quake invsqrt. Are you using a 32bit system?");
    float y  = number;
    float x2 = y * 0.5F;
    int i  = * ( int * ) &y;
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
  inline double exp_cawley(double y) {
    assert(2*sizeof(int)==sizeof(double) && "Approximate exp() requires 4-byte integer");
    //static double EXPA=1048576/std::log(2);
    union {
      double d;
      struct { int j, i; } n;  // little endian
      //struct { int i, j; } n;  // bin endian
    } eco;
    eco.n.i = 1072632447 + (int)(y*1512775.39519519);
    eco.n.j = 0;
    return eco.d;
  }

  inline double exp_untested(double y) {
    assert(2*sizeof(int)==sizeof(double) && "Approximate exp() requires 4-byte integer");
    double d;
    *((int*)(&d) + 0) = 0;
    *((int*)(&d) + 1) = (int)(1512775 * y + 1072632447);
    return d;
  }
#pragma GCC diagnostic pop

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

