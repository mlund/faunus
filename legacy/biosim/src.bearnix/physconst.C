#include <cmath>
#include "physconst.h"

physconst::physconst(double kelvin, double dielec) {
  T   = kelvin;       ///! Temperature in kelvin
  e_r = dielec;       ///! Relative permittivity
  pi  = acos(-1.);
  k   = 1.380658e-23; ///! Boltzmann's constant
  e   = 1.602177e-19; ///! Electron charge
  e_0 = 8.85419e-12;  ///! Vacuum permittivity
  Na  = 6.022137e23;  ///! Avogadros number

  beta= 1./(k*T);     ///! Thermal energy
  ecf = e*e / (4.*pi*e_0*e_r*1e-10);
  beta_ecf = beta * ecf; ///! The Bjerrum length
};
