#pragma once
#include "core.h"
#include "geometry.h"
#include "group.h"
#include "molecule.h"
#include <range/v3/view/join.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/sample.hpp>
#include <unordered_set>

namespace Faunus {

/**
 * @brief
 *
 * holds indices of either particles or groups and remembers which are active/inactive
 *
 *
 *
 */
class Tracker {
  public:
    using id_type = size_t;
    using index_type = size_t;

  private:
    std::vector<std::unordered_map<index_type, index_type>> active, inactive; //
    std::vector<std::vector<index_type>> index_holder;

    std::vector<size_t> number_of_active;

    void activate(id_type, size_t index);
    void deactivate(id_type, size_t);
    void insert(id_type, size_t);

    void removeInactive(id_type, size_t);
    void swap(id_type, id_type, size_t);

    size_t countActive(id_type) const;
    size_t countInactive(id_type) const;

    auto getActive(const id_type id) const {
        return ranges::subrange(index_holder.at(id).begin(), index_holder.at(id).begin() + countActive(id));
    }
    auto getInactive(const id_type id) const {
        return ranges::subrange(index_holder.at(id).begin() + countActive(id), index_holder.at(id).end());
    }
    void init(const id_type max_number_of_ids) {
        active.resize(max_number_of_ids);
        inactive.resize(max_number_of_ids);
        index_holder.resize(max_number_of_ids);
        number_of_active.resize(max_number_of_ids);
    }
    friend class Space; // I think there is no need for this class to be used outside space
};

template <typename reference_type> class Tracker2 {
  public:
    using id_type = size_t;

  private:
    std::unordered_map<reference_type, size_t> active, inactive; //
    std::vector<reference_type> index_holder;

    size_t number_of_active;

    void activate(id_type, reference_type);
    void deactivate(id_type, reference_type);
    void insert(id_type, reference_type);
    void removeInactive(id_type, reference_type);
    void swap(id_type, id_type, reference_type);

    size_t countActive(id_type) const;
    size_t countInactive(id_type) const;

    auto getActive(const id_type id) const {
        return ranges::subrange(index_holder.at(id).begin(), index_holder.at(id).begin() + countActive(id));
    }
    auto getInactive(const id_type id) const {
        return ranges::subrange(index_holder.at(id).begin() + countActive(id), index_holder.at(id).end());
    }
    void init(const id_type max_number_of_ids) {
        active.resize(max_number_of_ids);
        inactive.resize(max_number_of_ids);
        index_holder.resize(max_number_of_ids);
    }
    friend class Space; // I think there is no need for this class to be used outside space
};

template <typename reference_type>
void Tracker2<reference_type>::activate(const id_type i3d, reference_type reference) {

    active.insert({index, number_of_active});
    assert(inactive.count(index) != 0);             // deactivating reference must be among inactive
    assert(number_of_active < index_holder.size()); // I cannot activate anything when everything is active

    auto old_reference_to_holder = inactive.at(index);
    assert(index_holder.size() > old_reference_to_holder);

    auto first_inactive = index_holder.at(number_of_active);
    index_holder.at(old_reference_to_holder) = first_inactive;
    index_holder.at(number_of_active) = index;

    assert(inactive.count(first_inactive) != 0);
    inactive.at(first_inactive) = old_reference_to_holder;

    inactive.erase(index);

    number_of_active++;
}
template <typename reference_type>
void Tracker2<reference_type>::deactivate(const id_type id2, const reference_type index) {

    inactive.insert({index, number_of_active - 1});
    assert(active.count(index) != 0); // deactivating reference must be among inactive
    assert(number_of_active > 0);     // I cannot deactivate anything when everything is inactive
    auto old_reference_to_holder = active.at(index);
    assert(index_holder.size() > old_reference_to_holder);

    auto last_active = index_holder.at(number_of_active - 1);
    index_holder.at(old_reference_to_holder) = last_active;
    index_holder.at(number_of_active - 1) = index;

    assert(active.count(last_active) != 0);
    active.at(last_active) = old_reference_to_holder;
    active.erase(index);

    number_of_active--;
}
template <typename reference_type> void Tracker2<reference_type>::insert(const id_type id, const reference_type index) {
    inactive.insert({index, index_holder.size()});
    index_holder.push_back(index);
}

template <typename reference_type>
void Tracker2<reference_type>::removeInactive(const id_type id, const reference_type index) {
    auto old_reference_to_holder = inactive.at(index); // swap deactivated index to end of original index holder
    auto last_inactive = index_holder.at(index_holder.at(id).size() - 1);

    assert(index_holder.size() > old_reference_to_holder);
    index_holder.at(old_reference_to_holder) = last_inactive;
    index_holder.pop_back();
    inactive.erase(index);
}

/**
 * @brief swaps index from one holder into another
 *
 * @param old_id of bucket where @param index originally resides
 * @param new_id of bucket where @param index will be placed
 * @param index is a reference to particle/group, which we are swapping
 *
 * @note needed for speciation only when species is Atomic.
 */
template <typename reference_type>
void Tracker2<reference_type>::swap(const id_type old_id, const id_type new_id, reference_type index) {

    insert(new_id, index);   // push to top of inactive indices of new holder
    activate(new_id, index); // activate the newly inserted inactive index

    deactivate(old_id, index);     // deactivate index in original holder
    removeInactive(old_id, index); // remove deactivated index from original holder
}
template <typename reference_type> size_t Tracker2<reference_type>::countActive(const id_type id) const {
    return number_of_active;
}
template <typename reference_type> size_t Tracker2<reference_type>::countInactive(const id_type id) const {
    return index_holder.size() - countActive(id);
}

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
    bool accepted = false;

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

    std::unique_ptr<Tracker> atom_tracker;
    std::unique_ptr<Tracker> molecule_tracker;
    std::unique_ptr<std::vector<Tracker2<size_t>>> molecule_tracker2;

    std::unique_ptr<Tracker> neutral_atom_tracker;
    std::unique_ptr<Tracker> neutral_molecule_tracker;

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

    std::vector<std::reference_wrapper<Group>>
    randomMolecules2(MoleculeData::index_type molid, Random& rand, size_t number_of_samples = 1,
                     Selection selection = Selection::ACTIVE); //!< Random group matching molid

    Space::GroupVector::iterator randomMolecule(MoleculeData::index_type molid, Random& rand,
                                                Space::Selection selection = Selection::ACTIVE);

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
     * @brief returns absolute index of particle in ParticleVector
     * @param particle
     */
    inline auto indexOf(const Particle& particle) const {
        return static_cast<size_t>(std::addressof(particle) - std::addressof(particles.at(0)));
    }
    inline auto indexOf(const GroupType& group) const {
        return static_cast<size_t>(std::addressof(group) - std::addressof(groups.at(0)));
    }

    void swapInTrackers(const size_t old_atom_id, const size_t new_atom_id, Particle& particle) {
        atom_tracker->swap(old_atom_id, new_atom_id, indexOf(particle));
    }

    void activate(const GroupType& group);
    void deactivate(const GroupType& group);
    void activateParticle(const Particle& particle);
    void deactivateParticle(const Particle& particle);

    Tracker& getMolecularTracker() { return *molecule_tracker; }
    Tracker& getAtomicTracker() { return *atom_tracker; }

    /**
     * @brief Returns a range of particle indices according to selection and 'atom'
     * @param atom_id Molecule id to look for
     * @param selection
     * @return unordered_set of particle indices
     */
    auto getAtoms(const AtomData::index_type atom_id, Selection selection = Selection::ACTIVE) const {
        switch (selection) {
        case (Selection::ACTIVE):
            return atom_tracker->getActive(atom_id);
        case (Selection::INACTIVE):
            return molecule_tracker->getInactive(atom_id);
        case (Selection::ACTIVE_NEUTRAL):
            return neutral_molecule_tracker->getActive(atom_id);
        case (Selection::INACTIVE_NEUTRAL):
            return neutral_molecule_tracker->getInactive(atom_id);
        default:
            throw ConfigurationError("invalid selection!");
        }
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
     * @brief Returns a set of molecule indices according to selection and 'molid'
     * @param molid Molecule id to look for
     * @param selection
     * @return unordered_set of molecule indices
     */
    auto getMoleculesIndices(const MoleculeData::index_type molid, Selection selection = Selection::ACTIVE) {
        switch (selection) {
        case (Selection::ACTIVE):
            return molecule_tracker->getActive(molid);
        case (Selection::INACTIVE):
            return molecule_tracker->getInactive(molid);
        case (Selection::ACTIVE_NEUTRAL):
            return neutral_molecule_tracker->getActive(molid);
        case (Selection::INACTIVE_NEUTRAL):
            return neutral_molecule_tracker->getInactive(molid);
        default:
            throw ConfigurationError("invalid selection!");
        }
    }

    auto getMolecules(const MoleculeData::index_type molid, Selection selection = Selection::ACTIVE) {
        return getMoleculesIndices(molid, selection) |
               ranges::views::transform([this](auto index) { return groups[index]; });
    }

    /**
     * @brief Find active atoms of type `atomid` (complexity: order N)
     * @param atomid Atom id to look for
     * @return Range of filtered particles
     */
    auto findAtoms(AtomData::index_type atomid) {
        return activeParticles() |
               ranges::cpp20::views::filter([atomid](const Particle& particle) { return particle.id == atomid; });
    }

    std::vector<size_t> randomMolecules(const size_t molid, const size_t number_to_sample, Random& slump,
                                        Selection selection = Selection::ACTIVE) {
        auto molecule_indices = getMoleculesIndices(molid, selection);
        if (countMolecules(molid, selection) == 0) {
            return {};
        }
        if (number_to_sample == 1) {
            auto random_index = slump.range<size_t>(0, countMolecules(molid, selection) - 1);
            return {molecule_indices[random_index]};
        }
        return molecule_indices | ranges::views::sample(number_to_sample, slump.engine) | ranges::to<std::vector>;
    }

    std::vector<size_t> randomAtoms(const size_t atom_id, const size_t number_to_sample, Random& slump,
                                    Selection selection = Selection::ACTIVE) {
        auto atoms = getAtoms(atom_id, selection);
        if (countAtoms(atom_id, selection) == 0) {
            return {};
        }
        if (number_to_sample == 1) {
            auto random_index = slump.range<size_t>(0, countAtoms(atom_id, selection) - 1);
            return {atoms[random_index]};
        }
        return atoms | ranges::views::sample(number_to_sample, slump.engine) | ranges::to<std::vector>;
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

    //!< Count particles specified by selection
    size_t countAtoms(AtomData::index_type, Selection selection = Selection::ACTIVE) const;

    //!< Count molecules specified by selection
    size_t countMolecules(MoleculeData::index_type, Selection selection = Selection::ACTIVE) const;

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

    void sync(Space& other, const Change& change); //!< Copy differing data from other Space using Change object

    void reserveTrackers();
    void initTrackers();

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