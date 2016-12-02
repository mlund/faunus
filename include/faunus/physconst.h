#ifndef FAU_PHYSCONST_H
#define FAU_PHYSCONST_H

#ifndef SWIG
#include <cmath>
#include <climits>

#endif

namespace Faunus
{
  /**
   * @brief Physical constants and parameters.
   *
   * This class stores physical constants, temperature
   * as well as various conversion functions. All data
   * and functions are `static` and can thus be called
   * anywhere without construction.
   */
  template<class Td>
  class PhysicalConstants
  {
  private:
      static Td _T;          //!< Temperature [K]
  public:
      PhysicalConstants( Td= 298.15 );
      static const Td
          pi,                  //!< The number pi
          infty,               //!< Infinity
          e0,                  //!< Permittivity of vacuum [C^2/(J*m)]
          kB,                  //!< Boltzmann's constant [J/K]
          e,                   //!< Absolute electronic unit charge [C] 
          R,                   //!< Molar gas constant [J/(K*mol)]
          Nav;                 //!< Avogadro's number [1/mol]
      static Td lB( Td );      //!< Bjerrum length [Aangstrom]
      static Td T();         //!< Return temperature [K]
      static Td kT();        //!< Returns k_bT [J]
      static void setT( Td );  //!< Set temperature [K]
  };

#ifdef __INTEL_COMPILER
  // needed due to constexpr bug in intel13 compiler. Fixed?
  template<class Td>
    const Td PhysicalConstants<Td>::infty=-std::log(0.);
#else
  template<class Td>
  const Td PhysicalConstants<Td>::infty = std::numeric_limits<Td>::infinity();
#endif

  template<class Td>
  const Td PhysicalConstants<Td>::pi = std::acos(-1.);

  template<class Td>
  const Td PhysicalConstants<Td>::e0 = 8.85419e-12;

  template<class Td>
  const Td PhysicalConstants<Td>::e = 1.602177e-19;

  template<class Td>
  const Td PhysicalConstants<Td>::kB = 1.380658e-23;

  template<class Td>
  const Td PhysicalConstants<Td>::Nav = 6.022137e23;

  template<class Td>
  const Td PhysicalConstants<Td>::R = kB * Nav;

  template<class Td>
  Td PhysicalConstants<Td>::_T = 298.15;

  template<class Td>
  PhysicalConstants<Td>::PhysicalConstants( Td temp ) { setT(temp); }

  template<class Td>
  void PhysicalConstants<Td>::setT( Td temp ) { _T = temp; }

  template<class Td>
  Td PhysicalConstants<Td>::T() { return _T; }

  template<class Td>
  Td PhysicalConstants<Td>::kT() { return kB * _T; }

  template<class Td>
  Td PhysicalConstants<Td>::lB( Td e_r )
  {
      return e * e / (4 * pi * e0 * e_r * 1e-10 * kB * _T);
  }

  typedef PhysicalConstants<double> pc;      //!< Typedef for PhysicalConstants

  /**
   * @brief Chemistry units
   *
   * This is the default and currently only unit system in Faunus.
   * By using floating-point literals one may specify properties in
   * arbitrary units that will automatically be converted to the
   * following:
   *
   * Property           | Unit
   * :----------------- | :--------------------------
   * Energy             | Thermal energy (kT)
   * Temperature        | Kelvin (K)
   * Length             | Angstrom (A) 
   * Charge             | Electron unit charge (e)
   * Dipole moment      | Electron angstrom (eA)
   * Concentration      | Particles / angstrom^3
   * Pressure           | Particles / angstrom^3
   * Angle              | Radians
   *
   * Example:
   *
   * ~~~~
   * std::cout << 1.0_nm;       // 10
   * std::cout << 180.0_deg;    // 3.1415...
   * std::cout << 25.0_C;       // 298.15
   * std::cout << 50.0_K;       // 50
   * ~~~~
   */
  namespace ChemistryUnits
  {
    /// Temperature in Kelvin
    constexpr long double operator "" _K( long double T ) { return T; }

    /// Temperature in Celcius
    constexpr long double operator "" _C( long double T ) { return 273.15 + T; }

    /// Dipole moment in Debye
    constexpr long double operator "" _Debye( long double mu ) { return mu * 0.20819434; }

    /// Dipole moment in electron angstrom
    constexpr long double operator "" _eA( long double mu ) { return mu; }

    /// Dipole moment in Coulomb meter
    constexpr long double operator "" _Cm( long double mu ) { return mu * 1.0_Debye / 3.33564e-30; }

    /// Length in Angstrom
    constexpr long double operator "" _angstrom( long double l ) { return l; }

    /// Length in meters
    constexpr long double operator "" _m( long double l ) { return l * 1e10; }

    /// Length in Bohr
    constexpr long double operator "" _bohr( long double l ) { return l * 0.52917721092; }

    /// Length in nanometers
    constexpr long double operator "" _nm( long double l ) { return l * 10; }

    /// Volume in litres
    constexpr long double operator "" _liter( long double v ) { return v * 1e27; } // -> angstrom^3

    constexpr long double operator "" _m3( long double v ) { return v * 1e30; } // -> angstrom^3

    /// Number of molecules in moles
    inline long double operator "" _mol( long double n ) { return n * PhysicalConstants<double>::Nav; } // -> particles

    /// Concentration in moles per liter
    inline long double operator "" _molar( long double c )
    {
        return c * 1.0_mol / 1.0_liter;
    } // -> particle / angstrom^3

    /// Concentration in millimoles per liter
    inline long double operator "" _mM( long double c )
    {
        return c * 1.0e-3_mol / 1.0_liter;
    } // -> particle / angstrom^3

    /// Pressure in Pascal
    inline long double operator "" _Pa( long double p )
    {
        return p / PhysicalConstants<double>::kT() / 1.0_m3;
    } // -> particle / angstrom^3

    /// Pressure in atmosphere
    inline long double operator "" _atm( long double p ) { return p * 101325.0_Pa; } // -> particle / angstrom^3

    /// Pressure in bar
    inline long double operator "" _bar( long double p ) { return p * 100000.0_Pa; } // -> particle / angstrom^3

    /// Angle in radians
    constexpr long double operator "" _rad( long double a ) { return a; }

    /// Angle in degrees
    constexpr long double operator "" _deg( long double a ) { return a * 3.14159265358979323846 / 180; }

    /// Energy in kT
    constexpr long double operator "" _kT( long double u ) { return u; }

    /// Energy in Joule
    inline long double operator "" _J( long double u ) { return u / PhysicalConstants<double>::kT(); }

    /// Energy in kJ/mol
    inline long double operator "" _kJmol( long double u )
    {
        return u / PhysicalConstants<double>::kT() / PhysicalConstants<double>::Nav * 1e3;
    }

    /// Energy in kcal/mol
    inline long double operator "" _kcalmol( long double u ) { return u * 4.1868_kJmol; }

    /// Energy in hartree
    inline long double operator "" _hartree( long double u ) { return u * 4.35974434e-18_J; }

  }//namespace

  // set default unit system
  using namespace ChemistryUnits;

}// namespace
#endif
