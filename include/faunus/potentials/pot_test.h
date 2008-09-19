#ifndef FAU_POT_TEST_H
#define FAU_POT_TEST_H
#include "faunus/potentials/base.h"
namespace Faunus {
  class pot_test : public pot_lj {
    public:
      pot_test( pot_setup &pot ) : pot_lj(pot.eps/pot.lB) { f=pot.lB; }
      string info();
      inline double pairpot(const particle &p1, const particle &p2) {
        register double r2=p1.sqdist(p2);
        register double a=p1.radius+p2.radius;
        a=a*a;
        if (r2<4*a) {
          a=a/r2;
          a=a*a*a;
          a=a*a/f;
        } else a=0;
        register double qq=p1.charge*p2.charge;
        return (qq!=0) ? qq/sqrt(r2)+a : a;
      }
  };
  typedef Faunus::pot_test T_pairpot;
}
#endif
