#ifndef FAU_ENERGY_SPRINGINTERACTIONEXFIELD_H
#define FAU_ENERGY_SPRINGINTERACTIONEXFIELD_H

#include "faunus/energy/springinteraction.h"
#include "faunus/potentials/pot_debyehuckel.h"

namespace Faunus {
  /*!
   * \brief  
   * \author Bjorn Persson
   * \date Lund 2009
   *
   * some description
   */
  template<class pot_debyehuckelXYcc> class springinteractionfield : public springinteraction<pot_debyehuckelXYcc> {
    private:
    public:
      using springinteraction<pot_debyehuckelXYcc>::pair;
      using springinteraction<pot_debyehuckelXYcc>::energy;
      springinteractionfield(inputfile &in) : springinteraction<pot_debyehuckelXYcc>(in) {
        springinteractionfield::name+="SOMETHING!";
      }
      // ENERGY WITH EXTERNAL CORRECTION
      double energy(const vector<particle> &p, const group &g ) {
        double u=interaction<pot_debyehuckelXYcc>::energy(p,g);
        for (int i=g.beg; i<=g.end; i++) {
          u += pair.expot(p[i]);
          u += pair.hphobpot(p[i]);
        }
        return u;
      }
/*      double energy(const vector<particle> &p, const group &g1 , const group &g2) {
        double u=interaction<pot_debyehyckelXYcc>energy(p,g1,g2);
        for (int i=g1.beg; i<=g1.end; i++)
          u += pair.expot(p[i]);
        for (int i=g2.beg; i<=g2.end; i++)
          u += pair.expot(p[i]);
        return u;
      }
*/
     // HYDROPHOBIC POTENTIAL
      double hydrophobic(const vector<particle> &p, const group &g) {
        double u=0;
        for (int i=g.beg; i<=g.end; i++)
          u+=pair.hphobpot(p[i]);
        return u;
      }
      // LJ-ENERGY
      double ljenergy(const vector<particle> &p, const group &g) {
        int s=p.size();
        double u=0;
        for (int i =g.beg; i<=g.end; i++) {
          for (int j=0; j<g.beg; j++)
            u+=pair.lj(p[i], p[j], pair.sqdist(p[i], p[j]));
          for (int j=g.end+1; j<s; j++)
            u+=pair.lj(p[i], p[j], pair.sqdist(p[i], p[j]));
        }
        return u*pair.f;
      }
      // EXTERNAL CORRECTION
      double extpot(const vector<particle> &p, const group &g) {
      double u=0;
      for (int i=g.beg; i<=g.end; i++) 
        u+= pair.expot(p[i]);
      return u;
      }
      string info() {
        std::ostringstream o;
        o << interaction<pot_debyehuckelXYcc>::info()<<std::endl;
        return o.str();
      }
  };
}//namespace
#endif
