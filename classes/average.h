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
 *
 * Example:\n
 * \code
 *   average x;
 *   x+=2.0;
 *   x+=6.0;
 *   cout << x.avg() ;   // --> 4.0
 * \endcode
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
class avgmatrix {
 private:
  int nx,ny;
  int x2i(T val) { return int((val-xbeg)/steps+.5); };
  int y2i(T val) { return int((val-ybeg)/steps+.5); };
 public:
  T xbeg,xend,ybeg,yend,steps;
  vector< vector< average<T> > > d;
  avgmatrix(T x1, T x2, T y1, T y2, T s=0.5) {
    xbeg=x1;
    xend=x2;
    ybeg=y1;
    yend=y2;
    steps=s;
    nx=x2i(xend);
    ny=y2i(yend);
    d.resize(nx);
    for (unsigned int i=0; i<nx; i++)
      d[i].resize(ny);
  }
  void add(T x, T y, T val) { d[ x2i(x) ][ y2i(y) ].add(val); }
  bool save(string file) {
    ofstream f(file.c_str());
    if (f) {
      f.precision(30);
      for (T x=xbeg; x<xend; x+=steps) {
	for (T y=ybeg; y<yend; y+=steps)
	  f << x << " " << y << " "
	       << d[x2i(x)][y2i(y)].avg() << endl;
	f << endl;
      };
      f.close();
      return true;
    };
    return false;
  }
};

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

/*!
 * \brief A class to handle a vector of averages. Restriction: x>=0.
 * \author Mikael Lund
 * \date Prague 2007.
 * \warning Try not to use this -- works, but xytable<float,average> would be better!
 */
template<class T>
class distribution {
  private:
    unsigned int x2i(T x) { return static_cast<unsigned int>(x/res+.5); }
    T res;                      // data resolution
    vector<average <T> > y;     // vector of averages
    string note;                // arbitrary comment
    unsigned long long cnt;     // counter for ALL
    bool hist;                  // true if a histogram
  public:
    distribution(T steps=.5, string comment="", bool histogram=false) {
      set(steps,comment,histogram);
    }
    void set(T steps, string comment, bool histogram=false) {
      cnt=0;
      res=steps;
      note=comment;
      hist=histogram;
    }
    T get(T x) { return y[ x2i(x) ].avg(); }
    void add(T x) {             // add to histogram
      add(x, 1);
      cnt++;
    }
    void add(T x, T val) {
      unsigned int i=x2i(x);    // add to average
      if (i>=y.size())
        y.resize(i+1);
      y[i]+=val;
    }
    string show(T beg=0,T end=-1) {
      ostringstream o;
      o.precision(4);
      o << "# " << note << endl << scientific;
      unsigned int i,n=(end>-1) ? x2i(end) : y.size();
      for (i=x2i(beg); i<n; i++)
        if (y[i].cnt>0)
          o << T(i)*res << " " << ( (hist==false) ?
            y[i].avg() : y[i].sum/T(cnt) ) << endl ;
      return o.str();
    }
};
#endif
