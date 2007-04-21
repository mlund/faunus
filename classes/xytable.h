#include <valarray>
#include <iostream>

/*!
 * \brief Class for a 2D table, XY data etc.
 * \author Mikael Lund
 *
 * The types of X and Y are arbitrarily chosen.\n
 * X values:
 *  - can be smaller than zero (see constructor)
 *  - are assumed equidistant
 *
 *  Example:\n
 *  \code
 *  xytable<float,int> y(0.5, -10);
 *  y.add( -3.1, 512);
 *  cout << y(-3.1); --> 512
 *  \endcode
 */
template <class TX, class TY>
class xytable {
  private:
    TX xres, xmin;
    std::valarray<TY> y; 
    int x2i(TX x) { return int( (x-xmin)/xres+0.5); }
  public:
    //! \param resolution Distance between x values
    //! \param xminimum Minimum x value (can be smaller than zero)
    xytable(TX resolution=0.5, TX xminimum=0) {
      xres=resolution;
      xmin=xminimum;
    }

    /*!
     * \brief Add a data point
     * \note Memory is automatically allocated.
     */
    void add(TX xval, TY yval) {
      int i=x2i(xval);
      if (i>=y.size())
        y.resize(i+1);
      y[i]=yval;
    }
    TY operator()(TX val) { return y[ x2i(val) ]; }
};

