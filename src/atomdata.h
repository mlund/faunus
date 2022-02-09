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
    using key_type = std::string;
    std::map<key_type, double> data; //!< arbitrary additional properties
    friend void to_json(json&, const InteractionData&);

  public:
    bool contains(const key_type& name) const;                 // like C++20 map::contains
    double at(const key_type& name) const;                     // like map::at()
    double& at(const key_type& name);                          // like map::at()
    void insert_or_assign(const key_type& name, double value); // like C++17 map::insert_or_assign
};

void to_json(json& j, const InteractionData& a);
void from_json(const json& j, InteractionData& a);
void from_single_use_json(SingleUseJSON& j, InteractionData& a);

/**
 * @brief Static properties for patchy sphero cylinders (PSC)
 */
class SpheroCylinderData {
  protected:
    friend void from_json(const json&, SpheroCylinderData&);
    friend void to_json(json&, const SpheroCylinderData&);

  public:
    enum PatchType {
        None = 0,                      //!< No patch
        Full = 1,                      //!< Patch runs the full length of the SC
        Capped = 2                     //!< Patch stops before the end caps
    };                                 //!< Not yet in use @todo replace `patch_type` integer
    double chiral_angle = 0;           //!< Chiral angle (rotation of patch) on PSC
    double length = 0;                 //!< Sphere-cylinder length
    double patch_angle = 0;            //!< Angle of attrative patch on PSC
    double patch_angle_switch = 0;     //!< Angle of angular switch on sides of patch on PSC
    double patch_attraction_range = 0; //!< Dist. to which attraction is flat on attrative patch on PSC
    double patch_cutoff = 0;           //!< Distance on which attraction switches to zero on PSC
    PatchType type = None;       //!< Patch type of spherocylinder
};

void from_json(const json& j, SpheroCylinderData& psc);
void to_json(json& j, const SpheroCylinderData& psc);

/**
 * @brief General properties for atoms
 */
class AtomData { // has to be a class when a constant reference is used
  public:
    using index_type = std::size_t; //!< Unsigned int used for atom id and atom indexing

  private:
    index_type _id = 0;
    friend void to_json(json&, const AtomData&);
    friend void from_json(const json&, AtomData&);

  public:
    std::string name;                   //!< Name
    double charge = 0;                  //!< Particle charge [e]
    double mw = 1;                      //!< Weight
    double sigma = 0;                   //!< Diameter for e.g Lennard-Jones etc. [angstrom]
                                        //!< Do not set! Only a temporal class member during the refactorization
    double activity = 0;                //!< Chemical activity [mol/l]
    double alphax = 0;                  //!< Excess polarisability (unit-less)
    double dp = 0;                      //!< Translational displacement parameter [angstrom]
    double dprot = 0;                   //!< Rotational displacement parameter [degrees]
    double tension = 0;                 //!< Surface tension [kT/Å^2]
    double tfe = 0;                     //!< Transfer free energy [J/mol/angstrom^2/M]
    Point mu = {0, 0, 0};               //!< Dipole moment unit vector
    double mulen = 0;                   //!< Dipole moment length
    bool hydrophobic = false;           //!< Is the particle hydrophobic?
    bool implicit = false;              //!< Is the particle implicit (e.g. proton)?
    InteractionData interaction;        //!< Arbitrary interaction parameters, e.g., epsilons in various potentials
    SpheroCylinderData sphero_cylinder; //!< Data for patchy sphero cylinders (PSCs)

    index_type& id();             //!< Type id
    const index_type& id() const; //!< Type id
};

void to_json(json& j, const AtomData& a);
void from_json(const json& j, AtomData& a);

/**
 * @brief Construct vector of atoms from json array
 *
 * Items are added to existing items while if an item
 * already exists, it will be overwritten.
 */
void from_json(const json& j, std::vector<AtomData>& atom_vector);

extern std::vector<AtomData> atoms; //!< Global instance of atom list

/**
 * @brief Finds the first element with a member attribute `name` matching the input.
 *
 * @param rng  a range of elements
 * @param name  a name to look for
 * @return an iterator to the first element, or `last` if not found
 * @see findAtomByName(), findMoleculeByName()
 */
template <class Trange> auto findName(Trange& rng, std::string_view name) {
    return std::find_if(rng.begin(), rng.end(), [&name](auto& i) { return i.name == name; });
}

/**
 * @brief An exception to indicate an unknown atom name in the input.
 */
struct UnknownAtomError : public GenericError {
    explicit UnknownAtomError(std::string_view atom_name);
};

/**
 * @brief Finds an atom by its name in the global Faunus atoms lexicon.
 *
 * The first matching atom is returned, or an UnknownAtomError is thrown when not found.
 *
 * @param name  an atom name to look for
 * @return an atom found
 * @throw UnknownAtomError  when no atom found
 */
AtomData& findAtomByName(std::string_view name);

/**
 * @brief Search for `name` in `database` and return `id()`
 * @tparam Trange Container of object having `.name` and `.id()` data members
 * @param database Iterable range having `.name` and `.id()` members
 * @param names Container with names to convert to id
 * @throw if names not found
 * @return Vector of ids matching `names`
 *
 * This is typically used with `Faunus::atoms` or `Faunus::molecules`
 * to lookup atom or molecule names and return them as id numbers.
 * If the string `*` occurs in `names`, the returned vector will be
 * a sequence containing all id's of the database, i.e.
 * `0, ..., database.size()-1`.
 */
template <class Trange> auto names2ids(Trange& database, const std::vector<std::string>& names) {
    using id_type = typename Trange::value_type::index_type;
    std::vector<id_type> index;
    index.reserve(names.size());
    for (const auto& name : names) {
        if (name == "*") { // wildcard selecting all id's
            index.resize(database.size());
            std::iota(index.begin(), index.end(), id_type(0));
            return index;
        }
        if (auto it = findName(database, name); it != database.end()) {
            index.template emplace_back(it->id());
        } else {
            throw std::out_of_range("name '" + name + "' not found");
        }
    }
    return index;
}

} // namespace Faunus
