#ifndef FAU_ENERGY_RF_H
#define FAU_ENERGY_RF_H

#include "faunus/energy.h"
#include "faunus/potentials/pot_test_rf.h"

namespace Faunus {
  /*!
   * \brief .... 
   * \todo Not finished
   * \warning Untested!
   * \author Mikael Lund and Bjorn Persson
   * \date Lund 2009
   *
   * some description
   */
  class reactionfield : public interaction<pot_test_rf> {
    private:
    public:
      reactionfield(inputfile &in) : interaction<pot_test_rf>(in) {
        reactionfield::name+="SOMETHING!";
      }
      // TOTAL ENERGY
      double potential(const vector<particle> &p, int i) {
        double phi=0;
        for (int j=0; j<i; j++)
          phi += pair.rf.phi(p[j], p[i]);
        for (int j=0; j<i; j++)
          phi += pair.rf.phi(p[j], p[i]);
        return pair.rf.f*phi;  //*lB ??
      }
      double energy(const vector<particle> &p ) {
        double u=interaction<pot_test_rf>::energy(p);
        for (int i=0; i<p.size(); i++)
          u+=pair.rf.f*pair.rf.selfenergy(p[i]);
        return u;
      }
      double energy(const vector<particle> &p, int i) {
        return pair.rf.f*pair.rf.selfenergy(p[i]) + interaction<pot_test_rf>::energy(p,i);
      }
      double energy(const vector<particle> &p, const group &g ) {
        double u=interaction<pot_test_rf>::energy(p,g);
        for (int i=g.beg; i<=g.end; i++)
          u += pair.rf.selfenergy(p[i]);
        return pair.rf.f*u;
      }
      string info() {
        std::ostringstream o;
        o << interaction<pot_test_rf>::info()
          << pair.info();
        return o.str();
      }
  };
}//namespace
#endif
