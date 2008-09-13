#ifndef FAU_GEOMETRY_H
#define FAU_GEOMETRY_H

#include "faunus/point.h"

namespace Faunus {
  class shape {
    public:
      double volume, area;
      inline int anint(double x) { return int(x>0 ? x+.5 : x-.5); }
      virtual void recalc()=0;
  };

  class cube : public shape {
    public:
      double len, len_inv, len_half; 
      void setvolume(double v) {
        len=pow(v,1/3.);
        recalc();
      }
      void setlen(double l) {
        len=l;
        recalc();
      }
      void recalc() {
        len_inv=1./len;
        len_half=len/2.;
        volume=len*len*len;
      }
      inline void boundary(point &p) {
        p.x=p.x-len*anint(p.x*len_inv);
        p.y=p.y-len*anint(p.y*len_inv);
        p.z=p.z-len*anint(p.z*len_inv);
      }
  };
}//namespace
#endif
