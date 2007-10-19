/*
 * Average class
 *
 * Convenient handling of average
 * values.
 *
 * M. Lund, 2002
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
class average {
  friend ostream &operator<<(ostream &, average &);
  public:
    vector<double> v;                    ///< Vector to store sets of length average::maxcnt
    unsigned int maxcnt;                          ///< Number of values in each set.
    double sum;                          ///< Sum of current set.
    unsigned long long int cnt;                   ///< Number of values in current set.
   
    average(int=10000);                  ///< Constructor. Specify average::maxcnt

    double operator+(average &);      
    double operator-(average &);      
    void operator+=(double);             ///< Add value to current set. Same is average::add
    void add(double);                    ///< Add value to current set.
    void info();
    void reset();                        ///< Clear all data
    double avg();                        ///< Return average of all sets in average::v
    double stddev();                     ///< If more than one set, return standard deviation.
};

class avgmatrix {
 private:
  int nx,ny;
  int x2i(double val) { return int((val-xbeg)/steps); };
  int y2i(double val) { return int((val-ybeg)/steps); };
  
 public:
  double xbeg,xend,ybeg,yend,steps;
  vector< vector<average> > d;
  avgmatrix(double x1, double x2, double y1, double y2, double s=0.5) {
    xbeg=x1;
    xend=x2;
    ybeg=y1;
    yend=y2;
    steps=s;
    nx=x2i(xend);
    ny=y2i(yend);
    d.resize(nx);
    for (int i=0; i<nx; i++)
      d[i].resize(ny);
  };

  void add(double x, double y, double val) {
    d[ x2i(x) ][ y2i(y) ].add(val);
  };

  bool save(string file) {
    ofstream f(file.c_str());
    if (f) {
      f.precision(30);
      for (double x=xbeg; x<xend; x+=steps) {
	for (double y=ybeg; y<yend; y+=steps)
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
