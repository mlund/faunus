#ifndef FAUNUS_FIELD_H
#define FAUNUS_FIELD_H

#include <faunus/potentials.h>

namespace Faunus {

  template<typename Tfield>
    class field {
      
    };

    MatrixXd getField(Energy::NonbondedVector<Potential::DipoleDipoleHS,Geometry::Cuboid> &pot, Space &spc) {
      int size = spc.p.size();
      double R1 = 0.0;          double R2 = 0.0;   double R3 = 0.0;
      Point E(0,0,0);           Point r(0,0,0);
      MatrixXd field(size,3);   field.setZero();

      for(int I =0; I < size; I ++) {
        for(int i = 0; i < size; i ++) {
          if(i == I) continue;
          r = pot.getGeometry().vdist(spc.p[i],spc.p[I]);
          R1 = 1.0/r.norm();   r = r*R1;   R2 = R1*R1;   R3 = R2*R1;
          E = E + spc.p[i].charge*R2*r;                                             // From charges
          E = E + (3.0*spc.p[i].mu.dot(r)*r - spc.p[i].mu)*spc.p[i].muscalar*R3;    // From dipoles
          field.row(I) = E;                     
        }
        E.setZero();
      }
      return field;
    }

} //namespace

#endif
