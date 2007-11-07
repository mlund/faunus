#ifndef _physconst_h
#define _physconst_h

#include <iostream>

//! Placeholder for physical constants and parameters.
class physconst {
  public:
    double k,     //!< Boltzmann's constant
           e,     //!< Electronic charge
           Na,    //!< Avogadro's number
           pi,    //!< pi
           e_0,   //!< Permittivity of vacuum
           T,     //!< Temperature
           e_r,   //!< Relative dielectric constant
           beta,  //!< 1/kT
           ecf,
           beta_ecf; //!< Bjerrum length [AA]
    physconst(double=298.15, double=80.);
};
#endif
