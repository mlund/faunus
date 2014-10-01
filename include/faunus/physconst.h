#ifndef FAU_PHYSCONST_H
#define FAU_PHYSCONST_H

#ifndef SWIG
#include <cmath>
#include <climits>
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
  template<class Td>
    class PhysicalConstants {
      private:
        static Td _T;          //!< Temperature [K]
      public:
        PhysicalConstants(Td=298.15);
        static const Td
          pi,                  //!< The number pi
          infty,               //!< Infinity
          e0,                  //!< Permittivity of vacuum [C^2/(J*m)]
          kB,                  //!< Boltzmann's constant [J/K]
          e,                   //!< Absolute electronic unit charge [C] 
          R,                   //!< Molar gas constant [J/(K*mol)]
          Nav;                 //!< Avogadro's number [1/mol]
        static Td lB(Td);      //!< Bjerrum length [Aangstrom]
        static Td kT2kJ(Td=1); //!< kT/molecule -> kJ/mol
        static Td kJ2kT(Td=1); //!< kJ/mol -> kT/molecule
        static Td T();         //!< Return temperature [K]
        static Td kT();        //!< Returns k_bT [J]
        static void setT(Td);  //!< Set temperature [K]
        static Td D2eA(Td=1);  //!< Debye to electron Angstrom
        static Td eA2Cm(Td=1); //!< Converts eÅ to SI-units Cm
        
        //!< Convert Å to atomic unit(Bohr). Convert length scales. E.g. Length -> dim=1, Area -> dim=2, Volume -> dim=3
        template<class T>
          static T Ang2Bohr(const T &Ang, int dim=1) { return (Ang*pow(1.88971616463,dim)); }
        //!< Electron Angstrom to Debye. E.g. @f$ \mu_{e\AA} -> \mu_{Debye} @f$ -> dim=1, @f$ \mu^2_{e\AA} -> \mu^2_{Debye} @f$ -> dim=2
        template<class T>
        static T eA2D(const T &eA, int dim=1) { return eA*pow(4.803204544369458,dim); }

        static Td kT2Hartree(Td=kT());    //!< Convert energy in kT to atomic unit(Hartree Energy)
    };


#ifdef __INTEL_COMPILER
  // needed due to constexpr bug in intel13 compiler. Fixed?
  template<class Td>
    const Td PhysicalConstants<Td>::infty=-std::log(0.);
#else
  template<class Td>
    const Td PhysicalConstants<Td>::infty=std::numeric_limits<Td>::infinity();
#endif

  template<class Td>
    const Td PhysicalConstants<Td>::pi=std::acos(-1.);

  template<class Td>
    const Td PhysicalConstants<Td>::e0=8.85419e-12;

  template<class Td>
    const Td PhysicalConstants<Td>::e=1.602177e-19;

  template<class Td>
    const Td PhysicalConstants<Td>::kB=1.380658e-23;

  template<class Td>
    const Td PhysicalConstants<Td>::Nav=6.022137e23;

  template<class Td>
    const Td PhysicalConstants<Td>::R=kB*Nav;

  template<class Td>
    Td PhysicalConstants<Td>::_T=298.15;

  template<class Td>
    PhysicalConstants<Td>::PhysicalConstants(Td temp) { setT(temp); }

  template<class Td>
    void PhysicalConstants<Td>::setT(Td temp) { _T=temp; }

  template<class Td>
    Td PhysicalConstants<Td>::T() { return _T; }

  template<class Td>
    Td PhysicalConstants<Td>::kT() { return kB*_T; }

  template<class Td>
    Td PhysicalConstants<Td>::lB(Td e_r) {
      return e*e / (4*pi*e0*e_r*1e-10*kB*_T);
    }

  template<class Td>
    Td PhysicalConstants<Td>::kT2kJ(Td u) { return u*kB*_T*Nav*1e-3; }

  template<class Td>
    Td PhysicalConstants<Td>::kJ2kT(Td u) { return u/kT2kJ(1); }

  template<class Td>
    Td PhysicalConstants<Td>::D2eA(Td D) { return 0.20819434*D; }
    
  template<class Td>
    Td PhysicalConstants<Td>::eA2Cm(Td eA) { return (eA*3.33564*(1e-30)/0.20819434); }
    
  // Converts to atomic units
  template<class Td>
    Td PhysicalConstants<Td>::kT2Hartree(Td E_kT) { return (E_kT/(4.3597441775*pow(10,-18))); }

  typedef PhysicalConstants<double> pc;      //!< Typedef for PhysicalConstants

  /**
   * @brief Chemistry units
   *
   * This is the default and currently only unit system in Faunus.
   * By using string literals one may specify properties in
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
   * std::cout << 25_C;         // 298.15
   * std::cout << 50_K;         // 50
   * ~~~~
   */
  namespace ChemistryUnits {
    /// Temperature in kelvin
    constexpr long double operator "" _K(long double T)
    { return T; }

    /// Temperature in celcius
    constexpr long double operator "" _C(long double T)
    { return 273.15+T; }

    /// Dipole moment in Debye
    constexpr long double operator "" _Debye(long double mu)
    { return mu * 0.20819434; }

    /// Length in angstrom
    constexpr long double operator "" _angstrom(long double l)
    { return l; }

    /// Length in meters
    constexpr long double operator "" _m(long double l)
    { return l * 1e10; }

    /// Length in nanometers
    constexpr long double operator "" _nm(long double l)
    { return l*10; }

    /// Volume in litres
    constexpr long double operator "" _liter(long double v)
    { return v*1e27; }

    /// Number of molecules in moles
    inline long double operator "" _mol(long double n)
    { return n*PhysicalConstants<double>::Nav; } // -> particles

    /// Concentration in moles per liter
    inline long double operator "" _molar(long double c)
    { return c * 1.0_mol / 1.0_liter; } // -> particle / angstrom^3

    /// Pressure in Pascal
    inline long double operator "" _Pa(long double p)
    { return p / PhysicalConstants<double>::kT() / 1.0_liter; }

    /// Pressure in atmosphere
    inline long double operator "" _atm(long double p)
    { return p*101325.0_Pa; }

    /// Pressure in bar
    inline long double operator "" _bar(long double p)
    { return p*100000.0_Pa; }

    /// Angle in radians
    constexpr long double operator"" _rad ( long double a )
    { return a; }

    /// Angle in degrees
    inline long double operator"" _deg ( long double a )
    { return a*PhysicalConstants<double>::pi/180; }

    /// Energy in kT
    constexpr long double operator"" _kT (long double u)
    { return u; }

    /// Energy in Joule
    inline long double operator"" _J (long double u)
    { return u / PhysicalConstants<double>::kT(); }

    /// Energy in kJ/mol
    inline long double operator"" _kJmol (long double u)
    { return u / PhysicalConstants<double>::kT() / PhysicalConstants<double>::Nav * 1e3; }

    /// Energy in kcal/mol
    inline long double operator"" _kcalmol (long double u)
    { return u * 4.1868_kJmol; }

    /// Energy in hartree
    inline long double operator"" _hartree (long double u)
    { return u * 4.35974434e-18_J; }

  }//namespace

  // set default unit system
  using namespace ChemistryUnits;

}// namespace
#endif


