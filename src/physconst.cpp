#include "faunus/physconst.h"
namespace Faunus {
  const double pc::pi=std::acos(-1.);
  const double pc::e0=8.85419e-12;
  const double pc::e=1.602177e-19;
  const double pc::kB=1.380658e-23;
  const double pc::Nav=6.022137e23;
  const double pc::R=kB*Nav;

  pc::pc(double temp) { T=temp; }
  double pc::lB(double e_r) const { return e*e / (4*pi*e0*e_r*1e-10*kB*T); }
  double pc::kT2kJ(double u) const { return u*kB*T*Nav*1e-3; }
 
  // GLOBAL
  pc pyc;
}
