#ifndef _xydata_h
#define _xydata_h
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include "average.h"

using namespace std;

/*!
 * \brief Class to handle a dataset of xy values with equidistant x-values.
 * \author Mikael Lund
 * \date March 2007
 * \warning This class will be replaced by \link xytable \endlink
 */
template <class T=double>
class xydata {
 private:
  T invres;
  vector<T> y;                  //!< y value vector
  unsigned short x2i(T);        //!< xvalue->index
  T i2x(int);                   //!< index->xvalue
  T i2y(int);                   //!< index->yvalue
 public:
  string comment;               //!< Arbitrary comment
  T xmin;                       //!< Minimum x-value
  T xmax;                       //!< Maximum x-value
  T res;                        //!< Data resolution
  void setup(T,T,T);            //!< Reset boundaries
  xydata(T=0.,T=0.,T=0.1);      //!< Constructor
  void add(T,T);                //!< Add a point
  bool add(vector<T> &, vector<T> &); //!< Add a high res. dataset
  inline T x2y(T);              //!< xvalue->yvalue
  void list();                  //!< Print info
};

template <class T>
unsigned short xydata<T>::x2i(T xin) {
  return static_cast<unsigned short>((xin-xmin)*invres+.5);
}
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
}

/*!
 *  \param min Minimum x-value
 *  \param max Maximum x-value
 *  \param resolution Data resolution
 */
template <class T>
xydata<T>::xydata(T min, T max, T resolution) {
  setup(min,max,resolution);
}

template <class T>
void xydata<T>::add(T xin, T yin) {
  int i=x2i(xin);
  y[i] = yin;
}

/*! Load a high resolution dataset and average it
 * to the current resolution.
 */
template <class T>
bool xydata<T>::add(vector<T> &xin, vector<T> &yin) {
  T d0,d=0;
  unsigned int i=0, n=xin.size()-1;
  average<T> meanx, meany;
  setup( xin[0], xin[n], res); 
  d0=abs( xin[1] - xin[0] );
  if (d0>res)
    return false;
  while (i<xin.size()) {
    meanx+=xin[i];
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
}

template <class T>
void xydata<T>::list() {
  for (unsigned int i=0; i<y.size(); i++)
    cout << i2x(i) << " " << i2y(i) << endl;
}
#endif
