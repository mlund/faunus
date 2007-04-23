#ifndef _average_h
#define _average_h

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <cmath>

using std::ostream;
using namespace std;

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

    /*! Example of average class
     * \example average-test.C
     */
};

/*!
 *  \brief Statistical functions
 *  \author Mikael Lund
 *  \date Prague, 2007
 */
template<class T>
class statistics {
  private:
    vector<T> x;
  public:
    void add(T);        //!< Add a data point to set
    T avg();            //!< Return average
    T stddev();         //!< Return standard deviation
};

template<class T>
void statistics<T>::add(T val) { x.push_back(val); }

template<class T>
T statistics<T>::avg() {
  int N=x.size();
  T s=0;
  for (int i=0; i<N; i++) s+=x[i];
  return s/N;
};

/*!
 * \f$ s_{N-1} = \sqrt{\frac{1}{N-1} \sum \left ( x_i - <x> \right )^2}  \f$
 */
template<class T>
T statistics<T>::stddev() {
  int N=x.size();
  if (N<2) return 0;
  T xm=avg();
  T sum=0;
  for (int i=0; i<N; i++)
    sum+=pow( x[i]-xm, 2 );
  return sqrt( sum / (N-1) );
}

template<class T>
average<T>::average() {
  reset();
}

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

template<class T>
void average<T>::operator+=(T x) { add(x); }

template<class T>
void average<T>::add(T x) {
  sum += x;
  cnt++;
}

template<class T>
T average<T>::avg() { return sum/cnt; }

#endif
