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

} // end of namespace
#endif
