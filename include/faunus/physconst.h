#ifndef FAU_PHYSCONST_H
#define FAU_PHYSCONTS_H
#include "faunus/common.h"
namespace Faunus {
  /*! \brief Physical constants and parameters.
  */
  class physconst {
    public:
      double k,     //!< Boltzmann's constant [J/K]
             e,     //!< Electronic charge [C]
             Na,    //!< Avogadro's number [1/mol]
             pi,    //!< pi
             e_0,   //!< Permittivity of vacuum
             T,     //!< Temperature [K]
             e_r,   //!< Relative dielectric constant
             beta,  //!< 1/kT [1/J] 
             lB;    //!< Bjerrum length [AA]
      physconst(double=298.15, double=80.);
      void lB_TO_T(double);
  };

  namespace phystest {
    static double pi=acos(-1.),
                  e0=8.85419e-19,
                  kB=1.380658e-23,
                  e=1.602177e-19,
                  Nav=6.022137e23;
  }
} // namespace
#endif

