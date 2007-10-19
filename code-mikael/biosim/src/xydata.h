#ifndef _xydata_h
#define _xydata_h
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include "average.h"

using namespace std;

// Class to handle a dataset of xy values with
// equidistant x-values.
// M. Lund, Prague 2007.

template <class T>
class xydata {
 private:
  T invres;
  vector<T> y;                  // y value vector
  unsigned short x2i(T);        // xvalue->index
  T i2x(int);                   // index->xvalue
  T i2y(int);                   // index->yvalue
 public:
  string comment;               // arbitrary comment
  T xmin,xmax,res;              // boundaries and resolution
  void setup(T,T,T);            // re-set boundaries
  xydata(T=0.0,T=0.0,T=0.1);    // constructor
  void add(T,T);                // add a point
  bool add(vector<T> &, vector<T> &,bool=false); // add a high res dataset
  inline T x2y(T);              // xvalue->yvalue
  void list();                  // print info
};

// I M P L E M E N T A T I O N

template <class T>
unsigned short xydata<T>::x2i(T xin) {
  return static_cast<unsigned short>((xin-xmin)*invres+0.5);
};
template <class T> T xydata<T>::i2x(int i) { return i*res+xmin; }
template <class T> T xydata<T>::i2y(int i) { return y[i]; }
template <class T> T xydata<T>::x2y(T x) { return y[ x2i(x) ]; }

// re-set xmin, xmax and data resolution
template <class T>
void xydata<T>::setup(T min, T max, T resolution) {
  xmin=min; //minimum xvalue
  xmax=max; //maximum xvalue
  res=resolution; //xsteps
  invres=1./res;   //for speed
  y.resize( int( (xmax-xmin)/res+0.5) ); //adjust y vector
};

// constructor
template <class T>
xydata<T>::xydata(T min, T max, T resolution) {
  setup(min,max,resolution);
};

// convert x-value to index
template <class T>
void xydata<T>::add(T xin, T yin) {
  int i=x2i(xin);
  y[i] = yin;
};

// Load a high resolution dataset and average it
// to the current resolution.
// Note the 'sqr' bool that can be used to square the
// loaded xvalues. Useful for distances to avoid sqrt().
template <class T>
bool xydata<T>::add(vector<T> &xin, vector<T> &yin, bool sqr) {
  T d0,d=0;
  unsigned int i=0, n=xin.size()-1;
  average meanx, meany;
  setup(
      (sqr) ? xin[0]*xin[0] : xin[0],
      (sqr) ? xin[n]*xin[n] : xin[n], res); 
  d0=abs( xin[1] - xin[0] );
  if (d0>res)
    return false;
  while (i<xin.size()) {
    meanx+=(sqr) ? xin[i]*xin[i] : xin[i];
    meany+=yin[i];
    i++;
    d+=d0;
    if (d>=d0) {
      add( meanx.avg(), meany.avg() );
      meanx.reset();
      meany.reset();
      d=0;
    }
  };
  xmax=i2x( y.size()-1 ) ;
  return true;
};

template <class T>
void xydata<T>::list() {
  for (unsigned int i=0; i<y.size(); i++)
    cout << i2x(i) << " " << i2y(i) << endl;
};

#endif
