#ifndef FAU_PHYSCONST_H
#define FAU_PHYSCONST_H
#include "faunus/common.h"
namespace Faunus {
  /*! \brief Physical constants and parameters.
   */
  class PhysicalConstants {
    public:
      PhysicalConstants(double temp=298.15);
      static const double
        pi,
        infty,                       //!< Infinity
        e0,                          //!< Permittivity of vacuum [C^2/(J*m)]
        kB,                          //!< Boltzmann's constant [J/K]
        e,                           //!< Electronic charge [C] 
        R,                           //!< Molar gas constant [J/(K*mol)]
        Nav;                         //!< Avogadro's number [1/mol]
      static double T;               //!< Temperature [K]
      static double lB(double=78.5); //!< Bjerrum length [Aangstrom]
      static double kT2kJ(double);   //!< kT/molecule -> kJ/mol
  };

  typedef PhysicalConstants pc;
  
} // namespace
#endif

