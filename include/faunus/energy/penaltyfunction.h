#ifndef faunus_penaltyfunction_h
#define faunus_penaltyfunction_h

/*!
 * \brief General class for penalty functions along a coordinate
 * \author Mikael Lund
 * \date Malmo, 2011
 *
 * This class stores a penalty function, v(x), along a given coordinate, x,
 * which could be a distance, angle, volume etc. Initially v(x) is zero for
 * all x. Each time the system visits x the update(x) function should be called
 * so as to add the penalty energy, du. In the energy evaluation, the coordinate
 * x should be associated with the extra energy v(x) as returned by energy(x).
 */

#include "faunus/common.h"
#include "faunus/histogram.h"

namespace Faunus {
  class penaltyfunction {
    private:
      unsigned int number_of_updates; //!< Number of times the penalty function has been updated
      double du;                  //!< Energy to update penalty function with
      double sf;                  //!< Factor to scale du with when duscale is called
      xytable<double,double> v;   //!< The penalty function
      xytable<double,int> cnt;    //!< Number of times the system is observed at x
    public:
      penaltyfunction(double, double, double, double, double); //!< Constructor (xmin,xmax,xresolution,penalty energy, scaling factor)
      double update(double);      //!< Update penalty function at x, return energy change
      double energy(double);      //!< Get penalty energy at x
      void reset();               //!< Reset penalty function (set to zero)
      void write(string);         //!< Write penalty function to disk
      void dump(string);          //!< Dump penalty function to disk
      void dump(string, int, string);  //!< Dump penalty function to disk
      void gofrdump(string);      //!< Dump gofr to disk
      void gofrdump(string, int, string); //!< Dump gofr to disk
      bool load(string);          //!< Load penalty function from disk
      bool gofrload(string);      //!< Load gofr from disk
      double scaledu();           //!< Scale the penalty energy
      string info();              //!< Information string
  };
} // namespace
#endif
