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
        T avg();                        ///< Return average
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
      sum=0;
      cnt=0;
    }
  template<class T>
    average<T> average<T>::operator+(const average &a) {
      average<T> r = *this;
      r.cnt += a.cnt;
      r.sum += a.sum;
      return r;
    }
  template<class T> average<T> & average<T>::operator+=(T x) {
    add(x);
    return *this;
  }
  template<class T>
    void average<T>::add(T x) {
      sum += x;
      cnt++;
    }
  template<class T>
    bool average<T>::operator==(const average &a) { return (*this==a); };

  /*!
   * \brief Class to collect average values
   *        - extended to include root-mean-square.
   * \todo Untested
   * \author Mikael Lund
   * \date 2002-2007
   */
  template<class T>
    class average_ext : public average<T> {
      using average<T>::cnt;
      using average<T>::sum;
      private:
        T sqsum;
      public:
        average_ext() { sqsum=0; }
        void add(T x) {
          average<T>::add(x);
          sqsum+=x*x;
        }
        T rms() { return sqrt(sqsum/cnt); } //!< Root-mean-square
        T stddev() { return sqrt( sqsum/cnt - pow(sum/cnt,2) ); } //!< Standard deviation (sigma)
    };
}
#endif
