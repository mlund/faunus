#ifndef _MATUBAYASHI
#define _MATUBAYASHI

#include "histogram.h"
#include "potentials.h"

//comment out this when done testing:
//typedef pot_coulomb T_pairpot;

class matubayashi {
  public:
    histogram h;
    matubayashi() : h(0.1, -20., 20.) {
      h.comment = "Matubayashi energy histogram";
    }
    void add(interaction<T_pairpot> &pot, vector<particle> &p, group &g) {
      for (int i=0; i<g.beg; i++) 
        h.add ( pot.energy(p,g,i) ); // particles before the group
      for (int i=g.end+1; i<p.size(); i++)
        h.add( pot.energy(p,g,i) );  // particles after the group
    }
};
#endif
