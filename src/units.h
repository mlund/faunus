#pragma once
#include <limits>
#include <string>
#include <cmath>
#include <memory>
#include <numeric>
#include <array>

namespace Faunus {

    /** @brief Physical constants */
    namespace PhysicalConstants {
        typedef double T; //!< Float size
        constexpr T infty = std::numeric_limits<T>::infinity(), //!< Numerical infinity
            epsilon_dbl = std::numeric_limits<T>::epsilon(),    //!< Numerical precision
            max_value = std::numeric_limits<T>::max(),          //!< Maximal (finite) representable value
            max_exp_argument =
                709.782712893384,    //!< Largest value exp() can take before overflow (hard-coded for double)
            pi = 3.141592653589793,  //!< Pi
            e0 = 8.85419e-12,        //!< Permittivity of vacuum [C^2/(J*m)]
            e = 1.602177e-19,        //!< Absolute electronic unit charge [C]
            kB = 1.380658e-23,       //!< Boltzmann's constant [J/K]
            Nav = 6.022137e23,       //!< Avogadro's number [1/mol]
            c = 299792458.0,         //!< Speed of light [m/s]
            R = kB * Nav;            //!< Molar gas constant [J/(K*mol)]
        extern T temperature;        //!< Temperature [K]
        static inline T kT() { return temperature * kB; } //!< Thermal energy [J]
        static inline T RT() { return temperature * R; }  //!< Thermal energy [J/mol]
        static inline T bjerrumLength(T epsilon_r) {
            return e * e / ( 4 * pi * e0 * epsilon_r * 1e-10 * kT());
        } //!< Bjerrum length [Å]

        static inline T relativeDielectricFromBjerrumLength(T bjerrumlength) {
            return bjerrumLength(bjerrumlength);
        } //!< Bjerrum length to a relative dielectric constant
    }

    namespace pc = PhysicalConstants;

    //! Calculate the ionic strength of a binary salt
    double ionicStrength(double, const std::array<int, 2> &);

    //! Calculate the Debye screening length for a binary salt
    double debyeLength(double, const std::array<int, 2> &, double);

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] ionicStrength") {
        using doctest::Approx;
        CHECK(ionicStrength(0.1, {1, 1}) == Approx(0.1));
        CHECK(ionicStrength(0.1, {2, 2}) == Approx(0.5 * (0.1 * 4 + 0.1 * 4)));
        CHECK(ionicStrength(0.1, {2, 1}) == Approx(0.5 * (0.1 * 4 + 0.2)));
        CHECK(ionicStrength(0.1, {1, 2}) == Approx(0.5 * (0.2 + 0.1 * 4)));
        CHECK(ionicStrength(0.1, {1, 3}) == Approx(0.5 * (0.3 + 0.1 * 9)));
        CHECK(ionicStrength(0.1, {3, 1}) == Approx(0.5 * (0.3 + 0.1 * 9)));
        CHECK(ionicStrength(0.1, {1, 3}) == Approx(0.5 * (0.3 + 0.1 * 9)));
        CHECK(ionicStrength(0.1, {2, 3}) == Approx(0.5 * (0.3 * 4 + 0.2 * 9)));

        SUBCASE("debyeLength") { CHECK(debyeLength(0.03, {1, 1}, 7) == Approx(17.7376102214)); }
    }
#endif

    /**
     * @brief Chemistry units
     *
     * String literals to convert to the following internal units:
     *
     * Property           | Unit
     * :----------------- | :--------------------------
     * energy             | thermal energy unit (kT)
     * temperature        | kelvin (K)
     * time               | picosecond (ps)
     * length             | ångström (Å)
     * charge             | electron unit charge (e)
     * dipole moment      | electron-ångström (eÅ)
     * concentration      | number of particles per cubic ångström (Å⁻³)
     * pressure           | number of particles per cubic ångström (Å⁻³)
     * mass               | dalton (Da)
     * angle              | radians (rad)
     *
     * The value of pressure in particle densities assumes the ideal gas law.
     */
    namespace ChemistryUnits {
        typedef long double T; //!< floating point size
        constexpr T operator"" _rad(T a) { return a; }                  //!< angle in radians
        constexpr T operator"" _deg(T a) { return a * pc::pi / 180; }   //!< angle in degrees (to radians)
        constexpr T operator"" _K(T temp) { return temp; }              //!< temperature in kelvins
        constexpr T operator"" _C(T temp) { return 273.15 + temp; }     //!< temperature in degrees Celcius (to kelvins)
        constexpr T operator"" _ps(T tau) { return tau; }               //!< time in picoseconds
        constexpr T operator"" _s(T tau) { return tau * 1e12; }         //!< time in seconds (to picoseconds)
        constexpr T operator"" _angstrom(T l) { return l; }             //!< length in ångströms
        constexpr T operator"" _nm(T l) { return l * 10; }              //!< length in nanometers (to ångströms)
        constexpr T operator"" _m(T l) { return l * 1e10; }             //!< length in meters (to ångströms)
        constexpr T operator"" _bohr(T l) { return l * 0.52917721092; } //!< length in Bohr radii (to ångströms)
        constexpr T operator"" _angstrom3(T l) { return l; }            //!< volume in cubic ångströms
        constexpr T operator"" _m3(T v) { return v * 1e30; }            //!< volume in cubic meters (to cubic ångströms)
        constexpr T operator"" _liter(T v) { return v * 1e27; }         //!< volume in liters (to cubic ångströms)
        constexpr T operator"" _Da(T m) { return m; }                   //!< mass in daltons
        constexpr T operator"" _kg(T m) { return m * 1e3 * pc::Nav; }   //!< mass in kilograms (to daltons)
        //! amount of substance in moles (to number of particles)
        constexpr T operator"" _mol(T n) { return n * pc::Nav; }
        //! molar concentration (to particle density in number of particles per cubic ångström)
        constexpr T operator"" _molar(T c) { return c * 1.0_mol / 1.0_liter; }
        //! millimolar concentration (to particle density in number of particles per cubic ångström)
        constexpr T operator"" _millimolar(T c) { return c * 1.0e-3_mol / 1.0_liter; }
        //! pressure in pascals (to particle density in number of particles per cubic ångström – assuming the ideal gas law)
        inline T operator"" _Pa(T p) { return p / pc::kT() / 1.0_m3; }
        //! pressure in atmospheres (to particle density in number of particles per cubic ångström – assuming the ideal gas law)
        inline T operator"" _atm(T p) { return p * 101325.0_Pa; }
        //! pressure in bars (to particle density in number of particles per cubic ångström – assuming the ideal gas law)
        inline T operator"" _bar(T p) { return p * 100000.0_Pa; }
        //! dipole moment in electron—ångströms
        constexpr T operator"" _eA(T mu) { return mu; }
        //! dipole moment in debyes (to electron-ångströms)
        constexpr T operator"" _debye(T mu) { return mu * 0.208194334424626; }
        //! dipole moment in coulomb-meters (to electron ångströms)
        constexpr T operator"" _Cm(T mu) { return mu * 1.0_debye / 3.335640951981520e-30; }
        //! energy in thermal energy units kT
        constexpr T operator"" _kT(T u) { return u; }
        //! energy in joules (to thermal energy units kT)
        inline T operator"" _J(T u) { return u / pc::kT(); }
        //! energy in hartrees (to thermal energy units kT)
        inline T operator"" _hartree(T u) { return u * 4.35974434e-18_J; }
        //! energy in kilojoules per mole (to thermal energy units kT per particle)
        inline T operator"" _kJmol(T u) { return u / pc::kT() / pc::Nav * 1e3; }
        //! energy in kilocalories per mole (to thermal energy units kT per particle)
        inline T operator"" _kcalmol(T u) { return u * 4.1868_kJmol; }
   }  // namespace ChemistryUnits
    using namespace ChemistryUnits;

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Units and string literals")
    {
        using doctest::Approx;
        pc::temperature = 298.15_K;
        CHECK( 1.0e-10_m == 1 );
        CHECK( (1/1.0_debye) == Approx(4.8032) );
        CHECK( 1.0_debye == Approx( 3.33564e-30_Cm ) );
        CHECK( 1.0_debye == Approx( 0.20819434_eA ) );
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

        std::string bracket( const std::string & );
    } //!< Unicode

} // end of Faunus namespace
