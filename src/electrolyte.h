#pragma once

#include <optional>
#include <vector>
#include <nlohmann/json.hpp>

namespace Faunus {

using json = nlohmann::json;

/**
 * @brief Stores information about salts for calculation of Debye screening length etc.
 *
 * In this context a "salt" is an arbitrary set of cations and anions, combined to form
 * a net-neutral compound. The object state is _temperature independent_.
 */
class Electrolyte
{
  private:
    double ionic_strength;
    double molarity;
    std::vector<int> valencies;

  public:
    Electrolyte(double molarity, const std::vector<int>& valencies);
    Electrolyte(double debye_length,
                double bjerrum_length); //!< Initialize from existing Debye and Bjerrum length
    [[nodiscard]] double ionicStrength() const; //!< Molar ionic strength (mol/l)
    [[nodiscard]] double
    debyeLength(double bjerrum_length) const; //!< Debye screening length in Angstrom
    [[nodiscard]] double getMolarity() const; //!< Input salt molarity (mol/l)
    [[nodiscard]] const std::vector<int>&
    getValencies() const; //!< Charges of each participating ion in the salt
};

void to_json(json& j, const Electrolyte& electrolyte);
std::optional<Electrolyte> makeElectrolyte(const json& j); //!< Create ionic salt object from json

} // namespace Faunus
