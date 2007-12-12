#ifndef _POT_HYDROPHOBIC_H
#define _POT_HYDROPHOBIC_H

#include "potentials.h"

/*! \brief Double vdW interaction for iodide and hydrophobic  groups 
 *  \author Mikael Lund
 */
class pot_hydrophobic : private pot_lj {
  private:
    double scale;
  public:
    double f;             //!< Factor to convert to kT
    pot_hydrophobic( pot_setup &pot ) : pot_lj( pot.eps/pot.lB ) {
      f=pot.lB;
      scale=pot.hydroscale;
    }
    inline double pairpot(particle &p1, particle &p2) {
      register double r2=p1.sqdist(p2), u=lj(p1,p2,r2);
      if (p1.id==particle::I && p2.hydrophobic==true)
        u=scale*u;
      if (p2.id==particle::I && p1.hydrophobic==true)
        u=scale*u;
      return u + p1.charge*p2.charge/sqrt(r2);
    }
    string info() {
      string s="Hydrophobic\n";
      return s;
    }
};

#endif
