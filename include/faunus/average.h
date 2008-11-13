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
      private:
        T sqsum;                        ///< Square sum
      public:
        average();
        T sum;                          ///< Total sum
        unsigned long long int cnt;     ///< Number of values
        T avg();                        ///< Return average
        T rms();                        ///< Root-mean-square
        T stdev();                      ///< Standard deviation
        virtual void add(T);            ///< Add value to current set.
        void reset();                   ///< Clear all data
        average & operator+=(T);        ///< Add value to current set. 
        average operator+(const average &);///< Add two averages
        bool operator==(const average &);  ///< Comparison operator
        //! Example of average class
        //! \example average-test.C
    };
  template<class T> average<T>::average() { reset(); }
  template<class T> T average<T>::avg() { return sum/cnt; }
  template<class T>
    void average<T>::reset() {
      sum=sqsum=0;
      cnt=0;
    }
  template<class T>
    average<T> average<T>::operator+(const average &a) {
      average<T> r = *this;
      r.cnt += a.cnt;
      r.sum += a.sum;
      r.sqsum+=a.sqsum;
      return r;
    }
  template<class T> average<T> & average<T>::operator+=(T x) {
    add(x);
    return *this;
  }
  template<class T>
    void average<T>::add(T x) {
      sum+=x;
      sqsum+=x*x;
      cnt++;
    }
  template<class T>
    bool average<T>::operator==(const average &a) { return (*this==a); }
  template<class T>
    T average<T>::rms() { return sqrt(sqsum/cnt); }
  template<class T>
    T average<T>::stdev() { return sqrt( sqsum/cnt - pow(sum/cnt,2) ); }
 
}//namespace
#endif
