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
        constexpr T
            infty = std::numeric_limits<T>::infinity(),      //!< Numerical infinity
            neg_infty = -std::numeric_limits<T>::infinity(), //!< Numerical negative infinity
            epsilon_dbl = std::numeric_limits<T>::epsilon(), //!< Numerical precision
            max_value = std::numeric_limits<T>::max(),       //!< Maximal (finite) representable value
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
     * mass               | grams per mole (g/mol)
     * angle              | radians (rad)
     *
     * The value of pressure in particle densities assumes the ideal gas law.
     */
    namespace ChemistryUnits {

    typedef double T;                                                       //!< floating point size
    constexpr T operator"" _rad(long double a) { return a; }                //!< angle in radians
    constexpr T operator"" _deg(long double a) { return a * pc::pi / 180; } //!< angle in degrees (to radians)
    constexpr T operator"" _K(long double temp) { return temp; }            //!< temperature in kelvins
    constexpr T operator"" _C(long double temp) {
        return 273.15 + temp;
    }                                                                 //!< temperature in degrees Celcius (to kelvins)
    constexpr T operator"" _ps(long double tau) { return tau; }       //!< time in picoseconds
    constexpr T operator"" _s(long double tau) { return tau * 1e12; } //!< time in seconds (to picoseconds)
    constexpr T operator"" _angstrom(long double l) { return l; }     //!< length in ångströms
    constexpr T operator"" _nm(long double l) { return l * 10; }      //!< length in nanometers (to ångströms)
    constexpr T operator"" _m(long double l) { return l * 1e10; }     //!< length in meters (to ångströms)
    constexpr T operator"" _bohr(long double l) { return l * 0.52917721092; } //!< length in Bohr radii (to ångströms)
    constexpr T operator"" _angstrom3(long double l) { return l; }            //!< volume in cubic ångströms
    constexpr T operator"" _m3(long double v) { return v * 1e30; }    //!< volume in cubic meters (to cubic ångströms)
    constexpr T operator"" _liter(long double v) { return v * 1e27; } //!< volume in liters (to cubic ångströms)
    constexpr T operator"" _gmol(long double m) { return m; }         //!< mass in grams per mole
    constexpr T operator"" _kg(long double m) {
        return m * 1e3 * pc::Nav;
    }                                                       //!< mass in kilograms per particle (to grams per mole)
    constexpr T operator"" _Da(long double m) { return m; } //!< mass in daltons per particle (to grams per mole)

    //! amount of substance in moles (to number of particles)
    constexpr T operator"" _mol(long double n) { return n * pc::Nav; }

    //! molar concentration (to particle density in number of particles per cubic ångström)
    constexpr T operator"" _molar(long double c) { return c * 1.0_mol / 1.0_liter; }

    //! millimolar concentration (to particle density in number of particles per cubic ångström)
    constexpr T operator"" _millimolar(long double c) { return c * 1.0e-3_mol / 1.0_liter; }

    //! pressure in pascals (to particle density in number of particles per cubic ångström – assuming the ideal gas law)
    inline T operator"" _Pa(long double p) { return p / pc::kT() / 1.0_m3; }

    //! pressure in atmospheres (to particle density in number of particles per cubic ångström – assuming the ideal gas
    //! law)
    inline T operator"" _atm(long double p) { return p * 101325.0_Pa; }

    //! pressure in bars (to particles per cubic ångström – assuming the ideal gas law)
    inline T operator"" _bar(long double p) { return p * 100000.0_Pa; }

    //! dipole moment in electron—ångströms
    constexpr T operator"" _eA(long double mu) { return mu; }

    //! dipole moment in debyes (to electron-ångströms)
    constexpr T operator"" _debye(long double mu) { return mu * 0.208194334424626; }

    //! dipole moment in coulomb-meters (to electron ångströms)
    constexpr T operator"" _Cm(long double mu) { return mu * 1.0_debye / 3.335640951981520e-30; }

    //! energy in kT (to kT)
    constexpr T operator"" _kT(long double u) { return u; }

    //! energy in joules (to thermal energy units kT)
    inline T operator"" _J(long double u) { return u / pc::kT(); }

    //! energy in hartrees (to kT)
    inline T operator"" _hartree(long double u) { return u * 4.35974434e-18_J; }

    //! energy in kilojoules per mole (to kT per particle)
    inline T operator"" _kJmol(long double u) { return u / pc::kT() / pc::Nav * 1e3; }

    //! energy in kilocalories per mole (to kT per particle)
    inline T operator"" _kcalmol(long double u) { return u * 4.1868_kJmol; }

    } // namespace ChemistryUnits

    using namespace ChemistryUnits;

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

        std::string bracket(std::string_view sv);
    } //!< Unicode

} // end of Faunus namespace
