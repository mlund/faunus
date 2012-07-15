#include "faunus/physconst.h"
namespace Faunus {
  const double PhysicalConstants::pi=std::acos(-1.);
  const double PhysicalConstants::infty=std::numeric_limits<double>::infinity();
  const double PhysicalConstants::e0=8.85419e-12;
  const double PhysicalConstants::e=1.602177e-19;
  const double PhysicalConstants::kB=1.380658e-23;
  const double PhysicalConstants::Nav=6.022137e23;
  const double PhysicalConstants::R=kB*Nav;

  double PhysicalConstants::_T=298.15;

  PhysicalConstants::PhysicalConstants(double temp) {
    setT(temp);
  }

  void PhysicalConstants::setT(double temp) { _T=temp; }

  double PhysicalConstants::T() { return _T; }

  double PhysicalConstants::lB(double e_r) { return e*e / (4*pi*e0*e_r*1e-10*kB*_T); }

  double PhysicalConstants::kT2kJ(double u) { return u*kB*_T*Nav*1e-3; }
}
