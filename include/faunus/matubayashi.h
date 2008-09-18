#ifndef FAU_MATUBAYASHI_H
#define FAU_MATUBAYASHI_H
#include "faunus/histogram.h"
#include "faunus/analysis.h"

namespace Faunus {
  /*!
   * Samples a histogram of pair energies of one
   * group with all other single particles.
   *
   * \author mikaek lund
   * \date Prague, 2008
   */
  class matubayashi : public analysis {
    public:
      histogram hist; //!< The histogram
      matubayashi(); 
      string info(); 
      void add(interaction<T_pairpot> &, vector<particle> &, group &);
  };
  matubayashi::matubayashi() : hist(0.02, -10., 10.) {
    runfraction=0.3;
    hist.comment="Matubayashi energy histogram";
  }
  void matubayashi::add(interaction<T_pairpot> &pot, vector<particle> &p, group &g) {
    for (int i=0; i<g.beg; i++)           // particles before the group
      hist.add ( pot.energy(p,g,i) );     // add energy to histogram      
    for (int i=g.end+1; i<p.size(); i++)  // ...and after the group
      hist.add( pot.energy(p,g,i) );
  }
  string matubayashi::info() {
    std::ostringstream o;                                                        
    o << endl << "# MAUBAYASHI HISTOGRAM:" << endl
      << "#   More information:  Paper ref?" << endl;
    return o.str();
  }
}
#endif
