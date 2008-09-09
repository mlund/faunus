/*     Simbox class (for the use of box geometry)
 *     Martin Trulsson (2007) 
 *     A rebuild from Mikaels Lund's Cell class
 *     Periodic conditions lies in the modified point and particle classes
 */
#ifndef _simbox_h
#define _simbox_h
#include "faunus/physconst.h"

namespace Faunus {
  class Simbox : private Physconst {
    public:
      double vol_A, vol_L;                          // Volume in cubic angstrom and liter
      Simbox(double boxlen, double boxsep);         // Constructor, gives the dimensions
      double len, sep, area;                        // Dimensions of the box
      double len_half, sep_half;                    
      double sqlen_half;                            
      double conc(int);                             // Returns concentrations
      //  double ionStr( vector<particle> & );
      //  inline double segmentVolume(double,double);
  };
}
#endif
