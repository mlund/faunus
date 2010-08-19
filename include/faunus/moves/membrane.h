#ifndef FAU_MEMBRANE_H
#define FAU_MEMBRANE_H

#include "faunus/moves/base.h"
#include "faunus/histogram.h"

namespace Faunus {
  class popscmembrane;

  class membranemove : public markovmove {
    public:
    int steps;
    double rfend,rflatt,rfpopc,rfpops;
    double endacc,popcacc, popsacc, lattacc;
    int    endcnt,popccnt, popscnt, lattcnt;
    double dpmm, dpgp, dplatt, utemp;
    point dpv, dpvg;
    membranemove( ensemble&, container&, energybase&, popscmembrane&);
    double move(popscmembrane &);
    string info();
  };
}//namespace
#endif
