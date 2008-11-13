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

  class pc {
    public:
      pc(double temp=298.15) { T=temp; }
      const static double
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
} // namespace
#endif

