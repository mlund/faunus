#ifndef FAU_XYTABLE_H
#define FAU_XYTABLE_H

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
 * Examples\n
 * \code
 * #include "average.h"
 * #include "xytable.h"
 * ...
 * xytable<float,int> f(0.1, -10);
 * float x=-5.6;
 * f(x)=10;
 * cout << f(x); // --> 10
 * 
 * xytable<float, average<float> > h(0.5);
 * float r=7.5;
 * h(r)+=2;
 * h(r)+=4;
 * cout << h(r).avg(); // --> 3
 * \endcode
 *  
 * \author Mikael Lund
 * \todo Check if the STL resize command assign zero to new data
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
    typedef typename std::vector<TY>::iterator yiter;
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
      y.clear();
      xres=resolution;
      xmin=xminimum;
      if (xmaximum>0)
        y.resize( int(std::abs(xmaximum-xminimum)/resolution) );
    }

    //! Returns the x value with the highest y value 
    TX x_atmax_y() {
      yiter iter = std::max_element( y.begin(), y.end() );
      int i=iter-y.begin();  // convert iterator to index integer
      return i*xres + xmin;
    }
    
    //! Convenient data access
    TY& operator()(TX val) { return y.at( x2i(val) ); }

    //! Max x-value in set
    TX xmax() { return y.size()*xres+xmin; }
};
/*! 
 * \example xytable-test.C
 * Example of the xytable class
 */
}//namespace
#endif
