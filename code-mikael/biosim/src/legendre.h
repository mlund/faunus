#include<vector>
#include<iostream>

using namespace std;

// Class to handle a 2D table.
class table{
private:
  vector<double> y;
  double min,max,steps;
public:
  table( double beg, double end, double resolution) {
    min=beg;
    max=end;
    steps=resolution;
    y.resize(int(max/steps));
  };
  void set(double x, double val) {
    if (x>=min && x<=max)
      y[int(val/steps)]=val;
  };
  double get(double x) {
    if (x>=min && x<=max)
      return y[int(x/steps)];
  };
};

//Evaluate n'th degree Legendre polynomium at x.
//Result placed in vector p[0..n]
class legendre{
public:
  int n;
  vector<double> p;
  legendre(int order=-1) { resize(order); };
  void resize(int order) {
    n=order;
    if (n!=-1) {
      p.resize(n+1);
      p[0]=1.;
    };
  };
  inline void eval(double x) {
    if (n == 0) return;
    p[1] = x;
    double di;
    for( int i=1; i<n; i++ ) {
      di=double(i);
      p[i+1] = (2*di + 1) / (di+1) * x * p[i] - di/(di+1)*p[i-1];
    };
  }
};
