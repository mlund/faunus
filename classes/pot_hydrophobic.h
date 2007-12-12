#ifndef _POT_HYDROPHOBIC_H
#define _POT_HYDROPHOBIC_H

#include "potentials.h"

/*! \brief Double vdW interaction for iodide and hydrophobic  groups 
 *  \author Mikael Lund
 */
class pot_hydrophobic : private pot_lj {
  public:
    double f;             //!< Factor to convert to kT
    pot_hydrophobic( pot_setup &pot ) : pot_lj( pot.eps/pot.lB ) {
      f=pot.lB;
    }
    inline double pairpot(particle &p1, particle &p2) {
      register double r2=p1.sqdist(p2), u=lj(p1,p2,r2);
      if (p1.id==particle::I && p2.id==particle::HYDROPHOBIC)
        u=4.*u;
      if (p2.id==particle::I && p1.id==particle::HYDROPHOBIC)
        u=4.*u;
      return u + p1.charge*p2.charge/sqrt(r2);
    }
    string info() {
      string s="Hydrophobic\n";
      return s;
    }
    short renamehydrophobic(vector<particle> &, species &);
};

short pot_hydrophobic::renamehydrophobic( vector<particle> &p, species &s) {
  unsigned short i,cnt=0;
  for (i=0; i<p.size(); i++)
    if (s.d[p[i].id].hydrophobic==true) {
      p[i].id=particle::HYDROPHOBIC;
      cnt++;
    }
  return cnt;
}

#endif
