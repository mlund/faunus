#ifndef FAU_AVERAGE_H
#define FAU_AVERAGE_H

namespace Faunus {
  /*!
   * \brief Class to collect average values
   * \author Mikael Lund
   * \date 2002-2007
   */
  template<class T>
    class average {
      public:
        average();
        T sum;                          ///< Total sum
        unsigned long long int cnt;     ///< Number of values
        T avg();                        ///< Return average of all sets in average::v
        void add(T);                    ///< Add value to current set.
        void reset();                   ///< Clear all data
        void operator+=(T);             ///< Add value to current set. 
        average operator+(average &);   ///< Add two averages
        //! Example of average class
        //! \example average-test.C
    };
  template<class T> average<T>::average() { reset(); }
  template<class T> T average<T>::avg() { return sum/cnt; }
  template<class T>
    void average<T>::reset() {
      sum=0;
      cnt=0;
    }
  template<class T>
    average<T> average<T>::operator+(average &a) {
      average<T> r;
      r.cnt = cnt + a.cnt;
      r.sum = sum + a.sum;
      return r;
    }
  template<class T> void average<T>::operator+=(T x) { add(x); }
  template<class T>
    void average<T>::add(T x) {
      sum += x;
      cnt++;
    }
}
#endif
