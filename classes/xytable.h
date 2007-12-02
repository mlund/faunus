#ifndef _XYTABLE.H
#define _XYTABLE.H

#include <vector>
#include <iostream>
#include <string>

/*!
 * \brief Class for a 2D table, XY data etc.
* 
 * The types of X and Y are arbitrarily chosen.\n
 * X values:
 *  - can be smaller than zero (see constructor)
 *  - are assumed equidistant
 *
 * The in ternal data vector will be automatically resized to fit
 * the added data. However, for better memory utilization it is
 * possible to specify a maximum x value. (if this is exceeded the
 * vector will be resized anyway).
 *  
 * \author Mikael Lund
 * \todo Check if the STL resize command assign zero to new data
 *
 */
template <class TX, class TY>
class xytable {
  friend class histogram;
  friend class rdf;
  friend class rdfP3;
  private:
    int x2i(TX x) {
      int i=int( (x-xmin)/xres+0.5);
      if (i>=y.size())
        y.resize(i+1);
      return i;
    }
    std::vector<TY> y; 
  public:
    TX xres,            //!< Distance between x values
       xmin;            //!< Minimum x value

    //! \param resolution Distance between x values
    //! \param xminimum Minimum x value (can be smaller than zero)
    //! \param xmaximum Maximum x value (for better memory optimization)
    xytable(TX resolution=0.5, TX xminimum=0, TX xmaximum=0) {
      xres=resolution;
      xmin=xminimum;
      if (xmaximum>0)
        y.resize(xmaximum/resolution);
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

#endif
