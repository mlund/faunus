#include <cmath>
#include "physconst.h"

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
  ecf = e*e / (4.*pi*e_0*e_r*1e-10);
  beta_ecf = beta * ecf; 
}
