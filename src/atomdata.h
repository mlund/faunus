#pragma once
#include "core.h"

namespace Faunus {

/**
 * @brief A stub to hold various parameters of interactions.
 *
 * The class serves as a container for parameters of MixerPairPotentialBase. Currently only raw numbers
 * without units are stored. In the future, it is supposed to keep a track which parameters are needed
 * by potentials employed, as well as to convert units when initialized from JSON.
 */
class InteractionData {
    typedef std::string Tkey;
    std::map<Tkey, double> data; //!< arbitrary additional properties
    friend void to_json(json &, const InteractionData &);
  public:
    bool has(const Tkey name) const;               // like C++20 map::contains
    double get(const Tkey name) const;             // like map::at()
    double &get(const Tkey name);                  // like map::at()
    void set(const Tkey name, const double value); // like C++17 map::insert_or_assign
};

void to_json(json &j, const InteractionData &a);
void from_json(const json &j, InteractionData &a);
void from_single_use_json(SingleUseJSON &j, InteractionData &a);

/**
 * @brief General properties for atoms
 */
class AtomData { // has to be a class when a constant reference is used
  public:
    typedef int Tid;

  private:
    Tid _id = -1;
    friend void to_json(json &, const AtomData &);
    friend void from_json(const json &, AtomData &);

  public:
    std::string name;         //!< Name
    double charge = 0;        //!< Particle charge [e]
    double mw = 1;            //!< Weight
    double sigma = 0;         //!< Diameter for e.g Lennard-Jones etc. [angstrom]
                              //!< Do not set! Only a temporal class member during the refactorization
    double activity = 0;      //!< Chemical activity [mol/l]
    double alphax = 0;        //!< Excess polarisability (unit-less)
    double dp = 0;            //!< Translational displacement parameter [angstrom]
    double dprot = 0;         //!< Rotational displacement parameter [degrees]
    double mulen = 0;         //!< Dipole moment scalar [eÃ]
    double sclen = 0;         //!< Sphere-cylinder length [angstrom]
    double tension = 0;       //!< Surface tension [kT/Å^2]
    double tfe = 0;           //!< Transfer free energy [J/mol/angstrom^2/M]
    Point mu = {0, 0, 0};     //!< Dipole moment unit vector
    Point scdir = {1, 0, 0};  //!< Sphero-cylinder direction
    bool hydrophobic = false; //!< Is the particle hydrophobic?
    bool implicit = false;    //!< Is the particle implicit (e.g. proton)?
    InteractionData interaction; //!< Arbitrary interaction parameters, e.g., epsilons in various potentials

    Tid &id();                //!< Type id
    const Tid &id() const;    //!< Type id
};

void to_json(json &j, const AtomData &a);
void from_json(const json &j, AtomData &a);

/**
 * @brief Construct vector of atoms from json array
 *
 * Items are added to existing items while if an item
 * already exists, it will be overwritten.
 */
void from_json(const json &j, std::vector<AtomData> &v);

extern std::vector<AtomData> atoms; //!< Global instance of atom list

/**
 * @brief Finds the first element with a member attribute `name` matching the input.
 *
 * @param rng  a range of elements
 * @param name  a name to look for
 * @return an iterator to the first element, or `last` if not found
 * @see obtainName()
 */
template <class Trange> auto findName(Trange &rng, const std::string &name) {
    return std::find_if(rng.begin(), rng.end(), [&name](auto &i) { return i.name == name; });
}

/**
 * @brief Finds the first element with a member attribute `name` matching the input. Throws an exception if not found.
 *
 * @param rng  a range of elements
 * @param name  a name to look for
 * @return an iterator to the first element
 * @throw std::range_error if such an element not found
 * @see findName()
 */
template <class Trange> auto obtainName(Trange &rng, const std::string &name) {
    const auto result = findName(rng, name);
    if(result == rng.end()) {
        throw std::out_of_range(name + " not found");
    }
    return result;
}

/**
 * @brief Search for `name` in `database` and return `id()`
 * @tparam Trange Container of object having `.name` and `.id()` data members
 * @param database Iterable range having `.name` and `.id()` members
 * @param names Container with names to convert to id
 * @return Vector of ids matching `names`
 *
 * This is typically used with `Faunus::atoms` or `Faunus::molecules`
 * to lookup atom or molecule names and return them as id numbers.
 * If the string `*` occurs in `names`, the returned vector will be
 * a sequence containing all id's of the database, i.e.
 * `0, ..., database.size()-1`.
 */
template <class Trange> std::vector<int> names2ids(const Trange &database, const std::vector<std::string> &names) {
    std::vector<int> index;
    index.reserve(names.size());
    for (const auto &name : names) {
        if (name == "*") { // wildcard selecting all id's
            index.resize(database.size());
            std::iota(index.begin(), index.end(), 0);
            return index;
        } else if (auto it = findName(database, name); it != database.end())
            index.push_back(it->id());
        else {
            throw std::out_of_range("name '" + name + "' not found");
        }
    }
    return index;
}

} // namespace Faunus
