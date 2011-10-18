#include "faunus/physconst.h"
namespace Faunus {
  const double PhysicalConstants::pi=std::acos(-1.);
  const double PhysicalConstants::e0=8.85419e-12;
  const double PhysicalConstants::e=1.602177e-19;
  const double PhysicalConstants::kB=1.380658e-23;
  const double PhysicalConstants::Nav=6.022137e23;
  const double PhysicalConstants::R=kB*Nav;

  PhysicalConstants::PhysicalConstants(double temp) { T=temp; }
  double PhysicalConstants::T=298.15;
  double PhysicalConstants::lB(double e_r) { return e*e / (4*pi*e0*e_r*1e-10*kB*T); }
  double PhysicalConstants::kT2kJ(double u) { return u*kB*T*Nav*1e-3; }
}
