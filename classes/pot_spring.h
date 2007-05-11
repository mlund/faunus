#ifndef _POT_SPRING_H
#define _POT_SPRING_H

#include "point.h"

class pot_spring {
  public:
    double f,k,r0;
    pot_spring() {f=1;};
    inline double pairpot(particle &p1, particle &p2) {
      register double r=p1.dist(p2)-r0;
      return k*r*r;
    }
};

#endif
