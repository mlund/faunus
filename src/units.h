#pragma once
#include <limits>
#include <string>
#include <cmath>
#include <memory>

namespace Faunus {

    /** @brief Physical constants */
    namespace PhysicalConstants {
        typedef double T; //!< Float size
        constexpr T infty = std::numeric_limits<T>::infinity(), //!< Numerical infinity
                  pi = 3.141592653589793, //!< Pi
                  e0 = 8.85419e-12,  //!< Permittivity of vacuum [C^2/(J*m)]
                  e = 1.602177e-19,  //!< Absolute electronic unit charge [C]
                  kB = 1.380658e-23, //!< Boltzmann's constant [J/K]
                  Nav = 6.022137e23, //!< Avogadro's number [1/mol]
                  c = 299792458.0,   //!< Speed of light [m/s]
                  R = kB * Nav;      //!< Molar gas constant [J/(K*mol)]
        extern T temperature;        //!< Temperature (Kelvin)
        static inline T kT() { return temperature*kB; } //!< Thermal energy (Joule)
        static inline T lB( T epsilon_r ) {
            return e*e/(4*pi*e0*epsilon_r*1e-10*kT());
        } //!< Bjerrum length (angstrom)
        static inline T lB2epsr( T bjerrumlength ) {
            return lB(bjerrumlength);
        } //!< lB --> relative dielectric constant
    }

    namespace pc = PhysicalConstants;

    /**
     * @brief Chemistry units
     *
     * String literals to convert to the following internal units:
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
     */
    namespace ChemistryUnits
    {
        typedef long double T; //!< Floating point size
        constexpr T operator "" _K( T temp ) { return temp; } //!< Kelvin to Kelvin
        constexpr T operator "" _C( T temp ) { return 273.15 + temp; } //!< Celcius to Kelvin
        constexpr T operator "" _Debye( T mu ) { return mu * 0.208194334424626; }  //!< Debye to eA
        constexpr T operator "" _eA( T mu ) { return mu; } //!< eA to eA
        constexpr T operator "" _Cm( T mu ) { return mu * 1.0_Debye / 3.335640951981520e-30; } //!< Cm to eA
        constexpr T operator "" _angstrom( T l ) { return l; } //!< Angstrom to Angstrom
        constexpr T operator "" _angstrom3( T l ) { return l; } //!< Angstrom^3 to Angstrom^3
        constexpr T operator "" _m( T l ) { return l * 1e10; } //!< Meter to angstrom
        constexpr T operator "" _bohr( T l ) { return l * 0.52917721092; } //!< Bohr to angstrom
        constexpr T operator "" _nm( T l ) { return l * 10; } //!< nanometers to angstrom
        constexpr T operator "" _liter( T v ) { return v * 1e27; } //!< liter to angstrom cubed
        constexpr T operator "" _m3( T v ) { return v * 1e30; } //!< -> cubic meter to angstrom cubed
        constexpr T operator "" _mol( T n ) { return n * pc::Nav; } //!< moles to particles
        constexpr T operator "" _molar( T c ) { return c * 1.0_mol / 1.0_liter; } //!< molar to particle / angstrom^3
        constexpr T operator "" _mM( T c ) { return c * 1.0e-3_mol / 1.0_liter; } //!< millimolar to particle / angstrom^3
        constexpr T operator "" _rad( T a ) { return a; } //!< Radians to radians
        constexpr T operator "" _deg( T a ) { return a * pc::pi / 180; } //!< Degrees to radians
        inline T operator "" _Pa( T p ) { return p / pc::kT() / 1.0_m3; } //!< Pascal to particle / angstrom^3
        inline T operator "" _atm( T p ) { return p * 101325.0_Pa; } //!< Atmosphere to particle / angstrom^3
        inline T operator "" _bar( T p ) { return p * 100000.0_Pa; } //!< Bar to particle / angstrom^3
        constexpr T operator "" _kT( T u ) { return u; } //!< kT to kT (thermal energy)
        inline T operator "" _J( T u ) { return u / pc::kT(); } //!< Joule to kT
        inline T operator "" _kJmol( T u ) { return u / pc::kT() / pc::Nav * 1e3; } //!< kJ/mol to kT/particle
        inline T operator "" _kcalmol( T u ) { return u * 4.1868_kJmol; } //!< kcal/mol to kT/particle
        inline T operator "" _hartree( T u ) { return u * 4.35974434e-18_J; } //!< Hartree to kT
    }
    using namespace ChemistryUnits;

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Units and string literals")
    {
        using doctest::Approx;
        pc::temperature = 298.15_K;
        CHECK( 1.0e-10_m == 1 );
        CHECK( (1/1.0_Debye) == Approx(4.8032) );
        CHECK( 1.0_Debye == Approx( 3.33564e-30_Cm ) );
        CHECK( 1.0_Debye == Approx( 0.20819434_eA ) );
        CHECK( 360.0_deg == Approx( 2*std::acos(-1) ) );
        CHECK( (1.0_mol / 1.0_liter) == Approx( 1.0_molar ) );
        CHECK( 1.0_bar == Approx( 0.987_atm ) );
        CHECK( 1.0_atm == Approx( 101325._Pa) );
        CHECK( 1.0_kT == Approx( 2.47897_kJmol ) );
        CHECK( 1.0_hartree == Approx( 2625.499_kJmol ) );
    }
#endif

    namespace u8 {
        const std::string angstrom = "\u00c5";   //!< Angstrom symbol
        const std::string beta = "\u03b2";       //!< Greek beta
        const std::string cubed = "\u00b3";      //!< Superscript 3
        const std::string cuberoot = "\u221b";   //!< Cubic root
        const std::string degrees = "\u00b0";    //!< Degrees
        const std::string Delta = "\u0394";      //!< Greek Delta
        const std::string epsilon = "\u03f5";    //!< Greek epsilon
        const std::string epsilon_m = "\u03b5";  //!< Greek epsilon (minuscule)
        const std::string gamma = "\u0263";      //!< Greek gamma
        const std::string Gamma = "\u0393";      //!< Greek capital gamma
        const std::string infinity="\u221E";     //!< Infinity
        const std::string kappa = "\u03ba";      //!< Greek kappa
        const std::string mu = "\u03bc";         //!< Greek mu
        const std::string partial = "\u2202";    //!< Partial derivative
        const std::string percent = "\ufe6a";    //!< Percent sign
        const std::string pm = "\u00b1";         //!< Plus minus sign
        const std::string rho = "\u03C1";        //!< Greek rho
        const std::string rootof = "\u221a";     //!< Square root sign
        const std::string squared = "\u00b2";    //!< Superscript 2
        const std::string sigma = "\u03c3";      //!< Greek sigma
        const std::string superminus = "\u207b"; //!< Superscript minus (-)
        const std::string subr = "\u1D63";       //!< Subscript "r"
        const std::string theta = "\u03b8";      //!< Greek theta

        std::string bracket( const std::string &s );
    } //!< Unicode

} // end of Faunus namespace
