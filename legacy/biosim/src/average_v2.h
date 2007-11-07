/*
 * Average class
 *
 * Convenient handling of average
 * values.
 *
 * M. Lund, 2002-2007
 *
 */
#ifndef _average_h
#define _average_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <cmath>

using std::ostream;
using namespace std;

//! Class to collect average values.
/*!
 * This will automatically count the number of
 * values you add and can even return the standard
 * deviation (assuming the values are independent).
 * Example:\n
 *   average x;\n
 *   x+=2.0;\n
 *   x+=6.0;\n
 *   cout << x ;   // --> 4.0\n
 */

template<class T>
class average {
  friend ostream &operator<<(ostream &, average &);
  public:
    vector<T> v;                    ///< Vector to store sets of length average::maxcnt
    unsigned int maxcnt;            ///< Number of values in each set.
    unsigned long long int cnt;     ///< Number of values in current set.
    average(int=10000);             ///< Constructor. Specify average::maxcnt
    T sum;                          ///< Sum of current set.
    T operator+(average &);      
    T operator-(average &);      
    void operator+=(T);             ///< Add value to current set. Same is average::add
    void add(T);                    ///< Add value to current set.
    void info();
    void reset();                   ///< Clear all data
    T avg();                        ///< Return average of all sets in average::v
    T stddev();                     ///< If more than one set, return standard deviation.
};

template<class T>
class avgmatrix {
 private:
  int nx,ny;
  int x2i(T val) { return int((val-xbeg)/steps); };
  int y2i(T val) { return int((val-ybeg)/steps); };
  
 public:
  T xbeg,xend,ybeg,yend,steps;
  vector< vector<average> > d;
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
  };

  void add(T x, T y, T val) {
    d[ x2i(x) ][ y2i(y) ].add(val);
  };

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
  };

};

#endif

template<class T>
average<T>::average(int max) {
  maxcnt=max;
  sum=cnt=0;
  v.clear();
};

template<class T>
void average<T>::reset() {
  sum=cnt=0;
  v.clear();
};

template<class T>
T average<T>::operator+(average &a) { return avg() + a.avg(); };

template<class T>
T average<T>::operator-(average &a) { return avg() - a.avg(); };

template<class T>
void average<T>::operator+=(T x) { add(x); };

template<class T>
void average<T>::add(T x) {
  sum+=x;
  cnt++;
  if (cnt==maxcnt) {
    v.push_back( sum/cnt );
    sum=0;
    cnt=0;
  };
};

template<class T>
T average<T>::avg() {
  return (v.size()==0) ? sum/cnt :
    accumulate(v.begin(),v.end(),0.)/v.size();
};

template<class T>
T average<T>::stddev() {
  T vav=avg(), sum=0;
  for (unsigned int i=0; i<v.size(); i++)
    sum+=pow(v[i]-vav, 2);
  return sqrt( sum / ( v.size() * (v.size()-1) ) );  
};

//! Print average of all sets and more than one set, also the std. deviation.
ostream &operator<<( ostream &output, average &a) {
  output << a.avg();
  if (a.v.size()>1) output << " " << a.stddev();
  return output;
};

template<class T>
void average<T>::info() {
  cout << "# " << v.size() << " sets each containing " << maxcnt << " points." << endl;
};

