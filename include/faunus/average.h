#ifndef FAU_AVERAGE_H
#define FAU_AVERAGE_H

#include "faunus/common.h"

namespace Faunus {
  /*!
   * \brief Class to collect average values
   * \author Mikael Lund
   * \date 2002-2007
   */
  template<class T>
    class average {
      protected:
        //   vector v;
      public:
        T sqsum;                        ///< Square sum
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
        bool operator<(const average &);
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
    bool average<T>::operator < (const average &a) { return avg() < a.sum/a.cnt; }
  template<class T>
    T average<T>::rms() { return sqrt(sqsum/cnt); }
  template<class T>
    T average<T>::stdev() { return sqrt( sqsum/cnt - pow(sum/cnt,2) ); }


  /*!
   * \brief Class to keep track of block correlations
   * \author Mikael Lund
   * \date October 2009
   *
   * The sampling is performed in blocks of length n specified in
   * the constructor.
   *
   * \f$ c_i = \frac{ \langle x_0x_i\rangle_{i<n} - \langle x\rangle^2 }
   * { \langle x^2\rangle - \langle x\rangle^2  } \f$
   *
   * Example\n
   * \code
   * correlation<double> ci(50); // energy correlation
   * ci += uinit + du;           // place in MC loop
   * ...
   * for (int i=0; i<ci.size(); i++)
   *   cout << i << " " << ci[i] << end;
   * \endcode
   * ci will eventually fall off from one (full correlation) to
   * zero (uncorrelated).
   */

  template<class T>
    class correlation {
      private:
        unsigned int n,            //!< Length of each correlation measurement
                     cnt;          //!< Internal counter for each correlation set
        vector< average<T> > xixj; //!< Average correlation product, <xixj>
        average<T> xmean;          //!< Average values, <x>
        T xi;                      //!< Reference value (i=0) for each correlation set
      public:
        correlation(unsigned int=50);
        correlation & operator+=(T);    //!< Sample value
        T operator[] (unsigned int);    //!< Get correlation at i
        unsigned int size();            //!< Get block length
    };

  /*!
   * \param len Sample length
   */
  template<class T>
    correlation<T>::correlation(unsigned int len) {
      cnt=0;
      n = len; 
      xixj.resize(n);
    }

  template<class T>
    correlation<T> & correlation<T>::operator+=(T xj) {
      xmean+=xj;            // update total average
      if (cnt==n)           // end of block? reset counter.
        cnt=0;
      if (cnt==0)           // block start? save reference value.
        xi=xj;
      xixj.at(cnt) += xi*xj;// update average
      cnt++;
      return *this;
    }

  template<class T>
    T correlation<T>::operator[] (unsigned int i) {
      T xm=xmean.avg();
      T x2m=xmean.sqsum/xmean.cnt;
      return ( xixj.at(i).avg() - xm*xm ) / ( x2m - xm*xm ); 
    }

  template<class T>
    unsigned int correlation<T>::size() {
      return xixj.size();
    }

}//namespace
#endif
