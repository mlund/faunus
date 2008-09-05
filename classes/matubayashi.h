#ifndef FAU_MATUBAYASHI_H
#define FAU_MATUBAYASHI_H
#include "histogram.h"
#include "analysis.h"

namespace Faunus {
  /*!
   * Samples a histogram of pair energies of one
   * group with all other single particles.
   *
   * \author Mikael Lund
   * \date Prague, 2008
   */
  class matubayashi : public analysis {
    public:
      histogram hist; //!< The histogram
      matubayashi() : hist(0.1, -20., 20.) {
        hist.comment="Matubayashi energy histogram";
      }
      string info() { return "(info to come)"; }
      void add(interaction<T_pairpot> &pot, vector<particle> &p, group &g) {
        for (int i=0; i<g.beg; i++)           // particles before the group
          hist.add ( pot.energy(p,g,i) );     // add energy to histogram      
        for (int i=g.end+1; i<p.size(); i++)  // ...and after the group
          hist.add( pot.energy(p,g,i) );
      }
  };
}
#endif

