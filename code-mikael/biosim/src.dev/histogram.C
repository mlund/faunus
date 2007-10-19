#include <valarray>
#include <iostream>

using namespace std;

/*!
 * \brief Class for a 2D table, XY data etc.
 * \author Mikael Lund
 *
 * The types of X and Y are arbitrarily chosen.\n
 * X values:
 *  - can be smaller than zero (see constructor)
 *  - is assumed equidistant
 *
 *  Example:\n
 *  \code
 *  xy<float,int> y(0.5, -10);
 *  y.add( -3.1, 512);
 *  cout << y(-3.1); --> 512
 *  \endcode
 */
template <class TX, class TY>
class xy {
  private:
    TX xres, xmin;
    valarray<TY> y; 
    int x2i(TX x) { return int( (x-xmin)/xres+0.5); }
  public:
    //!< \param resolution Distance between x values
    //!< \param xminimum Minimum x value
    xy(TX resolution=0.5, TX xminimum=0) {
      xres=resolution;
      xmin=xminimum;
    }
    void add(TX xval, TY yval) {
      int i=x2i(xval);
      if (i>=y.size())
        y.resize(i+1);
      y[i]=yval;
    }
    TY operator()(TX val) { return y[ x2i(val) ]; }
};

int main() {
  xy<double,int> y(0.5,-10);
  y.add(-3.2, 10);
  cout << y(-3.2) << endl;
};
