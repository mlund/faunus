#include <vector>
#include <iostream>

/*!
 * \brief Class for a 2D table, XY data etc.
 * \author Mikael Lund
 * \todo Check if the STL resize command assign zero to new data
 * 
 * The types of X and Y are arbitrarily chosen.\n
 * X values:
 *  - can be smaller than zero (see constructor)
 *  - are assumed equidistant
 *
 */
template <class TX, class TY>
class xytable {
  private:
    int x2i(TX x) {
      int i=int( (x-xmin)/xres+0.5);
      if (i>=y.size())
        y.resize(i+1);
      return i;
    }
    std::vector<TY> y; // actual data
  public:
    TX xres,    //!< Distance between x values
       xmin;    //!< Minimum x value
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


/*! 
 * \example xytable-test.C
 * Example of the xytable class
 */
