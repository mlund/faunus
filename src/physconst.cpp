#include "faunus/physconst.h"
namespace Faunus {
  const double physical_const::pi=std::acos(-1.);
  const double physical_const::e0=8.85419e-12;
  const double physical_const::e=1.602177e-19;
  const double physical_const::kB=1.380658e-23;
  const double physical_const::Nav=6.022137e23;
  const double physical_const::R=kB*Nav;

  physical_const::physical_const(double temp) { T=temp; }
  double physical_const::T=298.15;
  double physical_const::lB(double e_r) { return e*e / (4*pi*e0*e_r*1e-10*kB*T); }
  double physical_const::kT2kJ(double u) { return u*kB*T*Nav*1e-3; }
}
