#ifndef FAU_PHYSCONST_H
#define FAU_PHYSCONST_H
#include "faunus/common.h"
namespace Faunus {
  /*! \brief Physical constants and parameters.
   */
  class pc {
    public:
      pc(double temp=298.15) { T=temp; }
      static const double
        pi,
        e0,                 //!< Permittivity of vacuum [C^2/(J*m)]
        kB,                 //!< Boltzmann's constant [J/K]
        e,                  //!< Electronic charge [C] 
        R,                  //!< Molar gas constant [J/(K*mol)]
        Nav;                //!< Avogadro's number [1/mol]
      double T;             //!< Temperature [K]
      double lB(double=80); //!< Bjerrum length [Aangstrom]
      double kT2kJ(double); //!< kT/molecule -> kJ/mol
  };
  
  extern pc pyc; // Global instance
  
} // namespace
#endif

