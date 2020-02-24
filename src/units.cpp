#include "units.h"

double Faunus::PhysicalConstants::temperature = 298.15;

std::string Faunus::u8::bracket(const std::string &s) {
    return "\u27e8" + s + "\u27e9";
}

/**
 * @param concentration Salt concentration (arbitrary units)
 * @param valency Absolute valencies, i.e. 1:1, 2:1, 2:2 salt
 * @return Ionic strength (same unit as `conc`)
 */
double Faunus::ionicStrength(double concentration, const std::array<unsigned char, 2> &valency) {
    double gcd = std::gcd(valency[0], valency[1]); // greatest common divider
    double n0 = valency[1] / gcd;                  // stoichiometric coefficient
    double n1 = valency[0] / gcd;                  // stoichiometric coefficient
    return 0.5 * concentration * (n0 * valency[0] * valency[0] + n1 * valency[1] * valency[1]);
}

/**
 * @param molarity Molar salt concentration
 * @param valency Absolute valencies, i.e. 1:1, 2:1, 2:2 salt
 * @param bjerrum_length Bjerrum length in Angstrom
 * @return Debye screening length in Angstrom
 */
double Faunus::debyeLength(double molarity, const std::array<unsigned char, 2> &valency, double bjerrum_length) {
    double ionic_strength = ionicStrength(molarity, valency) * 1.0_molar;
    return 1.0 / std::sqrt(8.0 * pc::pi * bjerrum_length * ionic_strength);
}
