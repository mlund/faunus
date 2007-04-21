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
 *  xytable<float,int> f(0.1, -10);
 *  float x=-5.6;
 *  f(x)=10;
 *  cout << f(x); --> 10
 *
 *  xytable<float, average<float> > h(0.5);
 *  float r=7.5;
 *  h(r)+=2;
 *  h(r)+=4;
 *  cout << h(r).avg(); --> 3
 *  \endcode
 *
 *  \todo Check if the STL resize command assign zero to new data
 */
template <class TX, class TY>
class xytable {
  private:
    TX xres, xmin;
    int x2i(TX x) {
      int i=int( (x-xmin)/xres+0.5);
      if (i>=y.size())
        y.resize(i+1); //hmm..automatically filled w. zero?
      return i;
    }
    std::valarray<TY> y; // actual data
  public:
    //! \param resolution Distance between x values
    //! \param xminimum Minimum x value (can be smaller than zero)
    xytable(TX resolution=0.5, TX xminimum=0) {
      xres=resolution;
      xmin=xminimum;
    }

    //! Convenient data access
    TY& operator()(TX val) { return y[ x2i(val) ]; }

    //! Max x-value in set
    TX xmax() { return y.size()*xres-xmin; }
};

