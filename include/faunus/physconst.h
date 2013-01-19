#ifndef FAU_PHYSCONST_H
#define FAU_PHYSCONST_H

#ifndef SWIG
#include <faunus/common.h>
#endif

namespace Faunus {
  /**
   * @brief Physical constants and parameters.
   *
   * This class stores physical constants, temperature
   * as well as various conversion functions. All data
   * and functions are `static` and can thus be called
   * anywhere without construction.
   */
  class PhysicalConstants {
    private:
      static double _T;               //!< Temperature [K]
    public:
      PhysicalConstants(double=298.15);
      static const double
        pi,                          //!< The number pi
        infty,                       //!< Infinity
        e0,                          //!< Permittivity of vacuum [C^2/(J*m)]
        kB,                          //!< Boltzmann's constant [J/K]
        e,                           //!< Absolute electronic unit charge [C] 
        R,                           //!< Molar gas constant [J/(K*mol)]
        Nav;                         //!< Avogadro's number [1/mol]
      static double lB(double);      //!< Bjerrum length [Aangstrom]
      static double kT2kJ(double);   //!< kT/molecule -> kJ/mol
      static double kJ2kT(double);   //!< kJ/mol -> kT/molecule
      static double T();             //!< Return temperature [K]
      static void setT(double);      //!< Set temperature [K]
      static double D2eA(double);    //!< Debye to electron Angstrom
  };

  typedef PhysicalConstants pc;      //!< Typedef for PhysicalConstants

} // namespace
#endif

