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
   * This is a template for storing permutable pairs of data `T`.
   * That is `(a,b)==(b,a)`. `operator<()` is implemented so the pairtype
   * can be used in STL maps etc.
   * @code
   * pair_permutable<int> a(2, 10);
   * pair_permutable<int> b(10, 2);
   * a==b;      // true
   * b.first;   // = 2
   * b.second;  // = 10
   * a.find(10);// true
   * a.find(3); // false
   * @endcode
   */
  template<class T> class pair_permutable {
    public:
      T first, second;
      pair_permutable() {}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-overflow"
      pair_permutable(T a, T b) : first(a), second(b) {
        if (first>second)
          std::swap(first,second);
      }
#pragma GCC diagnostic pop

      bool operator<(const pair_permutable<T> &a) const {
        assert(first<=second && a.first<=a.second);
        if (first<a.first)
          return true;
        if (first==a.first)
          if (second<a.second)
            return true;
        return false;
      }
    
      bool operator==(const pair_permutable<T> &a) const {
        assert(first<=second && a.first<=a.second);
        if (first==a.first)
          if (second==a.second)
            return true;
        return false;
      }

      bool find(const T &a) const {
        if (a!=first)
          if (a!=second)
            return false;
        return true;
      }
  };
}//namespace

/*
 * The following is an extension to std::hash in order to construct hash tables,
 * for example std::unordered_map, using Faunus::pair_permutable<> as keys.
 * See more at:
 * <http://en.wikipedia.org/wiki/Unordered_associative_containers_(C%2B%2B)#Usage_example>
 */
#ifdef FAU_HASHTABLE
namespace std {
  template<typename T>
    struct hash< Faunus::pair_permutable<T> > {
      size_t operator()(const Faunus::pair_permutable<T> &x) const {
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
  template<class Tbase, typename Tij=int, typename Tpair=pair_permutable< Tij > >
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
      public:
        std::vector<T> p;      //!< Legendre terms stored here
        //! Evaluate polynomium at x
        void eval(T x) {
          if (n > 0) {
            p[1] = x;
            T di;
            for (int i=1; i<n; ++i) {
              di=static_cast<T>(i);
              p[i+1] = (2*di + 1) / (di+1) * x * p[i] - di/(di+1)*p[i-1];
            }
          }
        }

        /**
         * @brief Constructor
         * @param order Order of the polynomium
         */
        legendre(int order=0) { resize(order); }

        void resize(int order) {
          n=order;
          assert(n>=0);
          p.resize(n+1);
          p[0]=1.;
        }
    };
}//namespace
#endif

