#pragma once
#include <limits>
#include <string>
#include <cmath>
#include <numbers>

namespace Faunus {

/** @brief Physical constants */
namespace PhysicalConstants {
using T = double;                                       //!< floating point type
constexpr T infty = std::numeric_limits<T>::infinity(), //!< Numerical infinity
    neg_infty = -std::numeric_limits<T>::infinity(),    //!< Numerical negative infinity
    epsilon_dbl = std::numeric_limits<T>::epsilon(),    //!< Numerical precision
    max_value = std::numeric_limits<T>::max(),          //!< Maximal (finite) representable value
    max_exp_argument =
        709.782712893384, //!< Largest value exp() can take before overflow (hard-coded for double)
    pi = std::numbers::pi,
            vacuum_permittivity = 8.85419e-12,          //!< Permittivity of vacuum [C^2/(J*m)]
    elementary_charge = 1.602177e-19,                   //!< Absolute electronic unit charge [C]
    boltzmann_constant = 1.380658e-23,                  //!< Boltzmann's constant [J/K]
    avogadro = 6.022137e23,                             //!< Avogadro's number [1/mol]
    speed_of_light = 299792458.0,                       //!< Speed of light [m/s]
    molar_gas_constant = boltzmann_constant * avogadro; //!< Molar gas constant [J/(K*mol)]
extern T temperature;                                   //!< Temperature [K]

static inline T kT()
{
    return temperature * boltzmann_constant;
} //!< Thermal energy [J]

static inline T RT()
{
    return temperature * molar_gas_constant;
} //!< Thermal energy [J/mol]

static inline T bjerrumLength(T relative_dielectric_constant)
{
    return elementary_charge * elementary_charge /
           (4 * pi * vacuum_permittivity * relative_dielectric_constant * 1e-10 * kT());
} //!< Bjerrum length [Å]

static inline T relativeDielectricFromBjerrumLength(T bjerrumlength)
{
    return bjerrumLength(bjerrumlength);
} //!< Bjerrum length to a relative dielectric constant
} // namespace PhysicalConstants

namespace pc = PhysicalConstants;

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

using T = double; //!< floating point size

constexpr T operator"" _rad(long double radians)
{
    return radians;
} //!< angle in radians

constexpr T operator"" _deg(long double degrees)
{
    return degrees * pc::pi / 180.0;
} //!< angle in degrees (to radians)

constexpr T operator"" _K(long double kelvin)
{
    return T(kelvin);
} //!< temperature in kelvins

constexpr T operator"" _C(long double celsius)
{
    return T(celsius) + 273.15;
} //!< temperature in degrees Celsius (to kelvins)

constexpr T operator"" _ps(long double picoseconds)
{
    return picoseconds;
} //!< time in picoseconds

constexpr T operator"" _s(long double secs)
{
    return T(secs) * 1e12;
} //!< time in seconds (to picoseconds)

constexpr T operator"" _angstrom(long double angstroms)
{
    return angstroms;
} //!< length in ångströms

constexpr T operator"" _nm(long double nanometers)
{
    return T(nanometers * 10.0);
} //!< length in nanometers (to ångströms)

constexpr T operator"" _m(long double meters)
{
    return T(meters * 1e10);
} //!< length in meters (to ångströms)

constexpr T operator"" _bohr(long double bohrs)
{
    return T(bohrs * 0.52917721092);
} //!< length in Bohr radii (to ångströms)

constexpr T operator"" _angstrom3(long double cubic_angstroms)
{
    return T(cubic_angstroms);
} //!< volume in cubic ångströms

constexpr T operator"" _m3(long double cubic_meters)
{
    return T(cubic_meters * 1e30);
} //!< volume in cubic meters (to cubic ångströms)

constexpr T operator"" _liter(long double liters)
{
    return T(liters * 1e27);
} //!< volume in liters (to cubic ångströms)

constexpr T operator"" _gmol(long double grams_per_mole)
{
    return grams_per_mole;
} //!< mass in grams per mole

constexpr T operator"" _kg(long double kilograms)
{
    return T(kilograms) * 1e3 * pc::avogadro;
} //!< mass in kilograms per particle (to grams per mole)

constexpr T operator"" _Da(long double daltons)
{
    return T(daltons);
} //!< mass in daltons per particle (to grams per mole)

//! amount of substance in moles (to number of particles)
constexpr T operator"" _mol(long double moles)
{
    return T(moles) * pc::avogadro;
}

//! molar concentration (to particle density in number of particles per cubic ångström)
constexpr T operator"" _molar(long double molarity)
{
    return T(molarity) * 1.0_mol / 1.0_liter;
}

//! millimolar concentration (to particle density in number of particles per cubic ångström)
constexpr T operator"" _millimolar(long double millimolar)
{
    return T(millimolar) * 1.0e-3_mol / 1.0_liter;
}

//! pressure in pascals (to particle density in number of particles per cubic ångström – assuming
//! the ideal gas law)
inline T operator"" _Pa(long double pascals)
{
    return T(pascals) / pc::kT() / 1.0_m3;
}

//! pressure in atmospheres (to particle density in number of particles per cubic ångström –
//! assuming the ideal gas law)
inline T operator"" _atm(long double atmospheres)
{
    return T(atmospheres) * 101325.0_Pa;
}

//! pressure in bars (to particles per cubic ångström – assuming the ideal gas law)
inline T operator"" _bar(long double bars)
{
    return T(bars) * 100000.0_Pa;
}

//! dipole moment in electron—ångströms
constexpr T operator"" _eA(long double electron_angstroms)
{
    return T(electron_angstroms);
}

//! dipole moment in debyes (to electron-ångströms)
constexpr T operator"" _debye(long double debyes)
{
    return T(debyes) * 0.208194334424626;
}

//! dipole moment in coulomb-meters (to electron ångströms)
constexpr T operator"" _Cm(long double coulomb_meters)
{
    return T(coulomb_meters) * 1.0_debye / 3.335640951981520e-30;
}

//! energy in kT (to kT)
constexpr T operator"" _kT(long double kT)
{
    return T(kT);
}

//! energy in joules (to thermal energy units kT)
inline T operator"" _J(long double joules)
{
    return T(joules) / pc::kT();
}

//! energy in hartrees (to kT)
inline T operator"" _hartree(long double hartrees)
{
    return T(hartrees) * 4.35974434e-18_J;
}

//! energy in kilojoules per mole (to kT per particle)
inline T operator"" _kJmol(long double kJ_per_mole)
{
    return T(kJ_per_mole) / pc::kT() / pc::avogadro * 1e3;
}

//! energy in kilocalories per mole (to kT per particle)
inline T operator"" _kcalmol(long double kcal_per_mole)
{
    return T(kcal_per_mole) * 4.1868_kJmol;
}

} // namespace ChemistryUnits

using namespace ChemistryUnits;

namespace unicode {
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
const std::string infinity = "\u221E";   //!< Infinity
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
} // namespace unicode

} // namespace Faunus
