#ifndef FAU_XYTABLE_H
#define FAU_XYTABLE_H

#include <vector>
#include <iostream>
#include <string>

#include "faunus/common.h"

namespace Faunus {
/*!
 * \brief Class for a 2D table, XY data etc.
* 
 * The types of X and Y are arbitrarily chosen.\n
 * X values:
 *  - can be smaller than zero (see constructor)
 *  - are assumed equidistant
 *
 * The internal data vector will be automatically resized to fit
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
  public:
    std::vector<TY> y; 
    TX xres,            //!< Distance between x values
       xmin;            //!< Minimum x value
    //! \param resolution Distance between x values
    //! \param xminimum Minimum x value (can be smaller than zero)
    //! \param xmaximum Maximum x value (for better memory optimization)
    xytable(TX resolution=.5, TX xminimum=0, TX xmaximum=0) {
      init(resolution, xminimum, xmaximum);
    }
    void init(TX resolution=.5, TX xminimum=0, TX xmaximum=0) {
      xres=resolution;
      xmin=xminimum;
      if (xmaximum>0)
        y.resize( int(std::abs(xmaximum-xminimum)/resolution) );
    }
    //! Convenient data access
    TY& operator()(TX val) { return y[ x2i(val) ]; }
    //! Max x-value in set
    TX xmax() { return y.size()*xres+xmin; }
};
/*! 
 * \example xytable-test.C
 * Example of the xytable class
 */
}//namespace
#endif
