#pragma once
#include "core.h"
#include "geometry.h"
#include "group.h"
#include "molecule.h"
#include <range/v3/view/join.hpp>

namespace Faunus {

/**
 * @brief Specify changes made to a system
 *
 * This class is used to describe a change to the system. It carries
 * information about changed group index, changed atom index, if the
 * volume has been changed, if there has been a particle insertion etc.
 * It does NOT describe the extent of the change, but simply points to
 * the affected groups, particles, etc.
 * The class is set when performing a system perturbation like a Monte Carlo
 * move and then passed on to the energy function where it is used to establish
 * which energies to calculate.
 *
 * Some notes:
 *
 * - If `GroupChange::all==true` then `relative_atom_indices` may be left empty. This is to
 *   avoid constucting a large size N vector of indices.
 */
struct Change {
    using index_type = std::size_t;
    bool everything = false;                 //!< Everything has changed (particles, groups, volume)
    bool volume_change = false;              //!< The volume has changed
    bool matter_change = false;              //!< The number of atomic or molecular species has changed
    bool moved_to_moved_interactions = true; //!< If several groups are moved, should they interact with each other?

    //! Properties of changed groups
    struct GroupChange {
        index_type group_index; //!< Touched group index
        bool dNatomic = false;  //!< The number of atomic molecules has changed
        bool dNswap = false;    //!< The number of atoms has changed as a result of a swap move
        bool internal = false;  //!< The internal energy or configuration has changed
        bool all = false;       //!< All particles in the group have changed (leave relative_atom_indices empty)
        std::vector<index_type> relative_atom_indices; //!< A subset of particles changed (sorted; empty if `all`=true)

        bool operator<(const GroupChange& other) const; //!< Comparison operator based on `group_index`
    };

    std::vector<GroupChange> groups; //!< Touched groups by index in group vector

    //! List of moved groups (index)
    inline auto touchedGroupIndex() const { return ranges::cpp20::views::transform(groups, &GroupChange::group_index); }

    //! List of changed atom index relative to first particle in system
    std::vector<index_type> touchedParticleIndex(const std::vector<Group>&) const;

    void clear();                                                             //!< Clear all change data
    bool empty() const;                                                       //!< Check if change object is empty
    explicit operator bool() const;                                           //!< True if object is not empty
    void sanityCheck(const std::vector<Group>& group_vector) const;           //!< Sanity check on contained object data
};

void to_json(json& j, const Change::GroupChange& group_change); //!< Serialize Change data to json
void to_json(json& j, const Change& change);                    //!< Serialise Change object to json

/**
 * @brief Placeholder for atoms and molecules
 *
 * Space is pervasive in the code as it stores the state of
 * - particles
 * - groups
 * - simulation container (`geometry`)
 * - implicit particles and groups
 *
 * It has methods to insert, find, and probe particles, groups, as
 * well as scaling the volume of the system
 */
class Space {
  public:
    using GeometryType = Geometry::Chameleon;
    using GroupType = Group; //!< Continuous range of particles defining molecules
    using GroupVector = std::vector<GroupType>;
    using ScaleVolumeTrigger = std::function<void(Space&, double, double)>;
    using ChangeTrigger = std::function<void(Space&, const Change&)>;
    using SyncTrigger = std::function<void(Space&, const Space&, const Change&)>;

  private:
    /**
     * @brief Stores implicit molecules
     *
     * This stores the number (map value) of each implicit
     * molecule in the system, identified by the molecular id (map key)
     * The reservoir is used to emulate a canonical system where a finite
     * number of implicit molecules can participate in equilibrium reactions.
     */
    std::map<MoleculeData::index_type, std::size_t> implicit_reservoir;

    std::vector<ChangeTrigger> changeTriggers; //!< Call when a Change object is applied (unused)
    std::vector<SyncTrigger> onSyncTriggers;   //!< Call when two Space objects are synched (unused)

  public:
    ParticleVector particles;                            //!< All particles are stored here!
    GroupVector groups;                                  //!< All groups are stored here (i.e. molecules)
    GeometryType geometry;                               //!< Container geometry (boundaries, shape, volume)
    std::vector<ScaleVolumeTrigger> scaleVolumeTriggers; //!< Functions triggered whenever the volume is scaled

    const std::map<MoleculeData::index_type, std::size_t>& getImplicitReservoir() const; //!< Implicit molecules
    std::map<MoleculeData::index_type, std::size_t>& getImplicitReservoir();             //!< Implicit molecules

    //!< Keywords to select particles based on the their active/inactive state and charge neutrality
    enum class Selection { ALL, ACTIVE, INACTIVE, ALL_NEUTRAL, ACTIVE_NEUTRAL, INACTIVE_NEUTRAL };

    void clear(); //!< Clears particle and molecule list
    GroupType& addGroup(MoleculeData::index_type molid, const ParticleVector& particles); //!< Append a group
    GroupVector::iterator findGroupContaining(const Particle& particle); //!< Finds the group containing the given atom
    GroupVector::iterator findGroupContaining(AtomData::index_type atom_index); //!< Find group containing atom index
    size_t numParticles(Selection selection = Selection::ACTIVE) const;         //!< Number of (active) particles

    Point scaleVolume(
        double new_volume,
        Geometry::VolumeMethod method = Geometry::VolumeMethod::ISOTROPIC); //!< Scales atoms, molecules, container

    GroupVector::iterator randomMolecule(MoleculeData::index_type molid, Random& rand,
                                         Selection selection = Selection::ACTIVE); //!< Random group matching molid

    json info();

    std::size_t getGroupIndex(const GroupType& group) const;         //!< Get index of given group in the group vector
    std::size_t getFirstParticleIndex(const GroupType& group) const; //!< Index of first particle in group
    std::size_t getFirstActiveParticleIndex(
        const GroupType& group) const; //!< Index of first particle w. respect to active particles

    /**
     * @brief Update particles in Space from a source range
     *
     * @tparam iterator Iterator for source range
     * @tparam copy_operation Functor used to copy data from an element in the source range to element in `Space:p`
     * @param begin Begin of source particle range
     * @param end End of source particle range
     * @param destination Iterator to target particle vector (should be in `Space::p`)
     * @param copy_function Function used to copy from source range to destination particle
     *
     * By default a source range of particles is expected and all content is copied.
     * This can be customised to e.g. a position range by giving a `copy_function` such
     * as `[](const auto &pos, auto &particle){particle.pos = pos;}`.
     *
     * The following is updated:
     *
     * - particles
     * - molecular mass centers of affected groups
     * - future: update cell list?
     *
     * @todo Since Space::groups is ordered, binary search could be used to filter
     */
    template <class iterator, class copy_operation = std::function<void(const Particle &, Particle &)>>
    void updateParticles(
        const iterator begin, const iterator end, ParticleVector::iterator destination,
        copy_operation copy_function = [](const Particle &src, Particle &dst) { dst = src; }) {

        const auto size = std::distance(begin, end); // number of affected particles

        assert(destination >= particles.begin() && destination < particles.end());
        assert(size <= std::distance(destination, particles.end()));

        auto affected_groups = groups | ranges::cpp20::views::filter([=](auto &group) {
                                   return (group.begin() < destination + size) && (group.end() > destination);
                               }); // filtered group with affected groups only. Note we copy in org. `destination`

        // copy data from source range (this modifies `destination`)
        std::for_each(begin, end, [&](const auto &source) { copy_function(source, *destination++); });

        for (auto& group : affected_groups) { // update affected mass centers
            group.updateMassCenter(geometry.getBoundaryFunc(), group.begin()->pos);
        }
    }

    //! Iterable range of all particle positions
    auto positions() const {
        return ranges::cpp20::views::transform(particles, [](auto& particle) -> const Point& { return particle.pos; });
    }

    //! Mutable iterable range of all particle positions
    auto positions() {
        return ranges::cpp20::views::transform(particles, [](auto& particle) -> Point& { return particle.pos; });
    }

    /**
     * @brief Finds all groups of type `molid` (complexity: order N)
     * @param molid Molecular id to look for
     * @param selection Selection
     * @return range with all groups of molid
     */
    auto findMolecules(MoleculeData::index_type molid, Selection selection = Selection::ACTIVE) {

        auto is_active = [](const GroupType& group) { return group.size() == group.capacity(); };

        auto is_neutral = [](auto begin, auto end) {
            auto charge =
                std::accumulate(begin, end, 0.0, [](auto sum, auto& particle) { return sum + particle.charge; });
            return (std::fabs(charge) < 1e-6);
        }; //!< determines if range of particles is neutral

        std::function<bool(const GroupType&)> f; //!< Lambda to filter groups according to selection

        switch (selection) {
        case (Selection::ALL):
            f = []([[maybe_unused]] auto& group) { return true; };
            break;
        case (Selection::INACTIVE):
            f = [=](auto& group) { return !is_active(group); };
            break;
        case (Selection::ACTIVE):
            f = [=](auto& group) { return is_active(group); };
            break;
        case (Selection::ALL_NEUTRAL):
            f = [=](auto& group) { return is_neutral(group.begin(), group.trueend()); };
            break;
        case (Selection::INACTIVE_NEUTRAL):
            f = [=](auto& group) { return !is_active(group) && is_neutral(group.begin(), group.trueend()); };
            break;
        case (Selection::ACTIVE_NEUTRAL):
            f = [=](auto& group) { return is_active(group) && is_neutral(group.begin(), group.end()); };
            break;
        }
        f = [f, molid](auto& group) { return group.id == molid && f(group); };
        return groups | ranges::cpp20::views::filter(f);
    }

    auto activeParticles() { return groups | ranges::cpp20::views::join; }       //!< Range with all active particles
    auto activeParticles() const { return groups | ranges::cpp20::views::join; } //!< Range with all active particles

    /**
     * @brief Find active atoms of type `atomid` (complexity: order N)
     * @param atomid Atom id to look for
     * @return Range of filtered particles
     */
    auto findAtoms(AtomData::index_type atomid) {
        return activeParticles() |
               ranges::cpp20::views::filter([atomid](const Particle& particle) { return particle.id == atomid; });
    }

    /**
     * @brief Find active atoms of type `atomid` (complexity: order N)
     * @param atomid Atom id to look for
     * @return Range of filtered particles
     */
    auto findAtoms(AtomData::index_type atomid) const {
        return activeParticles() |
        ranges::cpp20::views::filter([atomid](const Particle& particle) { return particle.id == atomid; });
    }

    size_t countAtoms(AtomData::index_type atomid) const; //!< Count active particles

    /**
     * @brief Count number of molecules matching criteria
     * @param molid Molecule id to match
     * @tparam mask Selection mask based on `Group::Selectors`
     * @return Number of molecules matching molid and mask
     */
    template <unsigned int mask> auto numMolecules(MoleculeData::index_type molid) const {
        auto filter = [&](const GroupType& group) {
            return (group.id == molid) ? group.template match<mask>() : false;
        };
        return std::count_if(groups.begin(), groups.end(), filter);
    }

    void sync(const Space& other, const Change& change); //!< Copy differing data from other Space using Change object

}; // end of space

void to_json(json& j, const Space& spc);   //!< Serialize Space to json object
void from_json(const json &j, Space &spc); //!< Deserialize json object to Space

/**
 * @brief Insert molecules into space
 *
 * Stateless class that helps inserting molecules into Space, based on json input.
 * json input should be an array of molecules, for example:
 *
 * ~~~ yaml
 *     - salt:  { molarity: 0.1 }          # number calc. from volume
 *     - water: { N: 256, inactive: 6 }    # 250 active, 6 inactive
 *     - dummy: { N: 100, inactive: true } # all 100 inactive
 * ~~~
 */
class InsertMoleculesInSpace {
  private:
    static void insertAtomicGroups(MoleculeData& moldata, Space& spc, size_t num_molecules,
                                   size_t num_inactive_molecules);

    static void insertMolecularGroups(MoleculeData& moldata, Space& spc, size_t num_molecules, size_t num_inactive);

    static void setPositionsForTrailingGroups(Space& spc, size_t num_molecules, const Faunus::ParticleVector&,
                                              const Point&);

    static void insertImplicitGroups(const MoleculeData& moldata, Space& spc, size_t num_molecules);

    //! Get number of molecules to insert from json object
    static size_t getNumberOfMolecules(const json& j, double volume, const std::string& molecule_name);

    //! Get number of inactive molecules from json object
    static size_t getNumberOfInactiveMolecules(const json& j, size_t number_of_molecules);

    //!< Get position vector using json object
    static ParticleVector getExternalPositions(const json &j, const std::string &molname);

    //!< Aggregated version of the above, called on each item in json array
    static void insertItem(const std::string &molname, const json &properties, Space &spc);

    static void reserveMemory(const json& j, Space& spc);

  public:
    static void insertMolecules(const json& j, Space& spc);
}; // end of insertMolecules class

namespace SpaceFactory {

void makeNaCl(Space& space, size_t num_particles, const Geometry::Chameleon& geometry);  //!< Create a salt system
void makeWater(Space& space, size_t num_particles, const Geometry::Chameleon& geometry); //!< Create a water system

} // namespace SpaceFactory

} // namespace Faunus
