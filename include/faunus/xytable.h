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

      //! \todo Test for negative x-values
      inline int x2ifloor(TX x) {
        return int( invxres*(x-xmin) );
      }

      public:
      typedef typename std::vector<TY>::iterator yiter;
      std::vector<TY> y; 
      TX xres,            //!< Distance between x values
         invxres,         //!< Inverse distance between x values
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

      TX x(int i) { 
        return i*xres+xmin; 
      }

      /*! \brief Interpolation
       *  \todo Test for negative x-values
       */
      TY interpolate(TX val) {
        int i = x2ifloor(val);
        TY  y_low  = y[i];
        TX  x_low  = x(i);
        TY slope = (y[i+1]-y_low) / (x(i+1)-x_low);
        return y_low + slope * (val-x_low); 
      }

      //! Max x-value in set
      TX xmax() { return y.size()*xres+xmin; }

      std::vector<TY> gety() { return y; }

      std::vector<TX> getx() {
        std::vector<TX> v;
        for (int i=0; i<y.size; ++i)
          v.push_back( i*xres + xmin );
        return v;
      }

      //! Dump table to disk
      bool dumptodisk(string filename) {
        std::ofstream f(filename.c_str());
        if (f) {
          f.precision(6);
          f << "% xytable dump: tablesize, xmin, xres; ydata" << endl;
          f << y.size() << " " << xmin << " " << xres << endl;
          for (int i=0; i<y.size(); i+=1) {
            f << y[i] << endl; 
          }
          f.close();
          return true;
        }
        return false;
      }

      bool dumptodisk(string file, int num, string ext) {
        string filename;
        std::stringstream buffer;
        buffer << file << num << "." << ext;
        filename = buffer.str();
        return dumptodisk(filename);
      }

      //! Load table from disk
      bool loadfromdisk(string filename) {
        int n;
        string s;
        std::ifstream f(filename.c_str());
        if (f) {
          getline(f,s);
          f >> n >> xmin >> xres;
          y.resize(n);
          for (int i=0; i<y.size(); i+=1) {
            f >> y[i];
          }
          f.close();
          return true;
        }         
        return false;
      }
    };
}//namespace
#endif
