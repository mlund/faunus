#ifndef FAU_XYTABLE_H
#define FAU_XYTABLE_H

#ifndef SWIG
#include <faunus/common.h>
#endif

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
    class xytable_center {
      protected:
        int x2i(TX x) {
          int i=int( invxres*(x-xmin)+0.5);
          if (i>=y.size())
            y.resize(i+1);
          return i;
        }

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
        xytable_center(TX xminimum, TX xmaximum, TX resolution) {
          init(resolution, xminimum, xmaximum);
        }

        void init(TX xminimum, TX xmaximum, TX resolution) {
          y.clear();
          xres=resolution;
          invxres=1/xres;
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

        /*! \brief Linear interpolation
         *  \warning If one of the reference values is Inf, nan is returned
         */
        TY interpolate(TX val) {
          int i = x2ifloor(val);
          TY  y_low  = y[i];
          TY slope = (y[i+1]-y_low) * invxres;
          return y_low + slope * (val-x(i)); 
        }

        //! Max x-value in set
        TX xmax() { return (y.size()-1)*xres+xmin; }

        std::vector<TY> gety() { return y; }

        std::vector<TX> getx() {
          std::vector<TX> v;
          for (int i=0; i<y.size; ++i)
            v.push_back( i*xres + xmin );
          return v;
        }

        //! Dump table to disk
        bool save(string filename) {
          std::ofstream f(filename.c_str());
          if (f) {
            f.precision(6);
            f << "% xytable dump: tablesize, xmin, xres, xmax; ydata" << endl;
            f << y.size() << " " << xmin << " " << xres << " " << xmax() << endl;
            for (int i=0; i<y.size(); i+=1) {
              f << y[i] << endl; 
            }
            f.close();
            return true;
          }
          return false;
        }

        bool save(string file, int num, string ext) {
          string filename;
          std::stringstream buffer;
          buffer << file << num << "." << ext;
          filename = buffer.str();
          return save(filename);
        }

        //! Load table from disk
        bool load(string filename) {
          int n;
          string s;
          std::ifstream f(filename.c_str());
          if (f) {
            getline(f,s);             // ignore first input line
            f >> n >> xmin >> xres;   // extract number of points, xmin, xres
            getline(f,s);             // ignore aditional columns
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

  /*!
   * \brief Class for a 2D table, XY data, etc. #2
   * Based on xytable_center with some modifications, 
   * here x-values correspond to the lower bound of the slice
   * for which the value is stored instead of the middle 
   * of the slice
   * 
   * \warning Does not automatically increase table size
   *          but points to the last slice instead
   * \warning Function fails if x < xmin
   * \author Chris Evers, Lund 2011-05-13
   */
  template <class TX, class TY>
    class xytable {
      protected:
        inline int x2i(TX x) { return floor( invxres*(x-xmin)); }

      public:
        typedef typename std::vector<TY>::iterator yiter;
        std::vector<TY> y; 
        TX xres,            //!< Distance between x values
           invxres,         //!< Inverse distance between x values
           xmin;            //!< Minimum x value
        int imax;           //!< Index of maximum x value

        //! \param resolution Distance between x values
        //! \param xminimum Minimum x value (can be smaller than zero)
        //! \param xmaximum Maximum x value (for better memory optimization)
        xytable(TX xminimum, TX xmaximum, TX resolution) {
          init(resolution, xminimum, xmaximum);
        }

        void init(TX xminimum, TX xmaximum, TX resolution) {
          y.clear();
          xres=resolution;
          invxres=1/resolution;
          xmin=xminimum;
          imax=ceil(std::abs(xmaximum-xminimum)/resolution)-1;
          if (xmaximum>0)
            y.resize(imax+1); // last one is used if x=xmaximum
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

        /*! \brief Linear interpolation
         *  \warning If one of the reference values is Inf, nan is returned
         */
        TY interpolate(TX val) {
          int i = x2i(val);
          TY  y_low  = y[i];
          TY slope = (y[i+1]-y_low) * invxres;
          return y_low + slope * (val-x(i)); 
        }

        //! Max x-value in set
        TX xmax() { return imax*xres+xmin; }

        std::vector<TY> gety() { return y; }

        std::vector<TX> getx() {
          std::vector<TX> v;
          for (int i=0; i<y.size; ++i)
            v.push_back( i*xres + xmin );
          return v;
        }

        //! Dump table to disk
        bool save(string filename) {
          std::ofstream f(filename.c_str());
          if (f) {
            f.precision(6);
            f << "% xytable2 dump: tablesize, xmin, xres, xmax; ydata" << endl
              << y.size() << " " << xmin << " " << xres << " " << xmax() << endl;
            for (auto y_i : y)
              f << y_i << endl; 
            f.close();
            return true;
          }
          return false;
        }

        bool save(string file, int num, string ext) {
          string filename;
          std::stringstream buffer;
          buffer << file << num << "." << ext;
          filename = buffer.str();
          return save(filename);
        }

        bool load(string filename) {
          int n;
          string s;
          std::ifstream f(filename.c_str());
          if (f) {
            getline(f,s);
            f >> n >> xmin >> xres;
            y.resize(n);
            for (auto y_i : y)
              f >> y_i;
            f.close();
            return true;
          }         
          return false;
        }
    };
}//namespace
#endif
