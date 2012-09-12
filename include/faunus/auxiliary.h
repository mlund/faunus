#ifndef FAU_auxiliary
#define FAU_auxiliary

#ifndef SWIG
#include <utility>
#endif

namespace Faunus {

  /*!
   * This is a template for storing permutable pairs of data T.
   * That is (a,b)==(b,a). The < operator is implemented so the pairtype
   * can be used in STL maps etc.
   * \code
   * pair_permutable<int> a(2, 10);
   * pair_permutable<int> b(10, 2);
   * a==b;      // true
   * b.first;   // = 2
   * b.second;  // = 10
   * a.find(10);// true
   * a.find(3); // false
   * \endcode
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

      // NOT USED BY stl::map.find()
      // see http://www.velocityreviews.com/forums/t290085-std-map-mystring-mystring-comparison-operator.html
      bool operator==(const pair_permutable<T> &a) const {
        if (a.first==first)
          if (a.second==second)
            return true;
        if (a.first==second)
          if (a.second==first)
            return true;
        return false;
      }

      bool operator<(const pair_permutable<T> &a) const {
        if (first<a.first)
          return true;
        if (first==a.first)
          if (second<a.second)
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

  /*!
   * This is a list of pairs with associated data where the latter should de derived
   * from the base Tbase. When adding data with the add() function, a copy of the data
   * is created and stored internally.
   */
  template<class Tbase, typename Tij=int, typename Tpair=pair_permutable< Tij > >
    class pair_list {
      protected:
        std::map< Tpair, std::shared_ptr<Tbase> > list;
      public:
        template<typename Tderived>
          void add(Tij i, Tij j, Tderived data) {
            list[ Tpair(i,j) ] = std::shared_ptr<Tderived>( new Tderived(data) ); 
          }
        Tbase& operator() (Tij i, Tij j) {
          Tpair pair(i,j);
          assert( list[pair] != nullptr ); //debug
          return *list[pair];
        }
        void clear() {
          list.clear();
        }
    };


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
  /*!
   * \brief Quake inverse square root approximation
   */
  inline float invsqrtQuake(float number) {
    assert(sizeof(int)==4 && "Integer size must be 4 bytes for quake invsqrt. Are you using a 32bit system?");
    const float threehalfs = 1.5F;
    float x2 = number * 0.5F;
    float y  = number;
    int i  = * ( int * ) &y;
    i  = 0x5f3759df - ( i >> 1 );
    y  = * ( float * ) &i;
    y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
    //      y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
    return y;
  }
#pragma GCC diagnostic pop

//#define EXPA (1048576/0.693147180559945309417)
//#define EXPC 60801
  /*!
   * \brief Approximate exp() function
   * \note see Cawley 2000; doi:10.1162/089976600300015033
   * \warning Does not work in big endian systems!
   */
  inline double exp_cawley(double y) {
    static double EXPA=1048576/std::log(2);
    static double EXPC=60801;
    union {
      double d;
      struct { int j, i; } n;  // little endian
      //struct { int i, j; } n;  // bin endian
    } eco;
    eco.n.i = (int) (EXPA*(y)) + (1072693248 - EXPC);
    eco.n.j = 0;
    return eco.d;
  }

} // end of namespace
#endif
