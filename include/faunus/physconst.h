#ifndef FAU_physconst_h
#define FAU_physconst_h
#include "faunus/scopes.h"
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
}
#endif

