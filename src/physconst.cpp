#include <cmath>
#include "physconst.h"
namespace Faunus {
//! \param kelvin Temperature in Kelvin
//! \param dielec Relative dielectric constant
physconst::physconst(double kelvin, double dielec) {
  T   = kelvin; 
  e_r = dielec; 
  pi  = acos(-1.);
  k   = 1.380658e-23; 
  e   = 1.602177e-19; 
  e_0 = 8.85419e-12; 
  Na  = 6.022137e23; 
  beta= 1./(k*T); 
  lB  = e*e / (4.*pi*e_0*e_r*1e-10) * beta;
}
void physconst::lB_TO_T(double bjerr) {       //Calculate the absolute T (K) from lB
  T=e*e*1e10/(4.*pi*e_0*e_r*k*bjerr);    //assuming e_r has been set
  beta=1./(k*T);
  lB=bjerr;
}
}
