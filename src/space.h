#pragma once
#include "core.h"
#include "geometry.h"
#include "group.h"
#include "molecule.h"
#include <range/v3/view/join.hpp>

namespace Faunus {

class Space;

/**
 * @brief Specify change to a new state
 *
 * - If `moved` or `removed` are defined for a group, but are
 *   empty, it is assumed that *all* particles in the group are affected.
 */
struct Change {
    bool dV = false;         //!< Set to true if there's a volume change
    bool all = false;        //!< Set to true if *everything* has changed
    bool dN = false;         //!< True if the number of atomic or molecular species has changed
    bool moved2moved = true; //!< If several groups are moved, should they interact with each other?

    struct data {
        bool dNatomic = false;  //!< True if the number of atomic molecules has changed
        bool dNswap = false;    //!< True if the number of atoms has changed as a result of a swap move
        int index;              //!< Touched group index
        bool internal = false;  //!< True if the internal energy/config has changed
        bool all = false;       //!< True if all particles in group have been updated
        std::vector<int> atoms; //!< Touched atom index w. respect to `Group::begin()`

        bool operator<(const data &other) const;
    }; //!< Properties of changed groups

    std::vector<data> groups; //!< Touched groups by index in group vector

    //! List of moved groups (index)
    inline auto touchedGroupIndex() {
        return ranges::cpp20::views::transform(groups, [](data &i) -> int { return i.index; });
    }

    //! List of changed atom index relative to first particle in system
    std::vector<int> touchedParticleIndex(const std::vector<Group<Particle>> &);

    void clear();                                                 //!< Clear all change data
    bool empty() const;                                           //!< Check if change object is empty
    explicit operator bool() const;                               //!< True if object is not empty
    void sanityCheck(const std::vector<Group<Particle>> &) const; //!< Sanity check on contained object data
};

void to_json(json &, const Change::data &); //!< Serialize Change data to json
void to_json(json &, const Change &);       //!< Serialise Change object to json

/**
 * @brief Placeholder for atoms and molecules
 */
class Space {
  public:
    typedef Geometry::Chameleon Tgeometry;
    typedef Particle Tparticle; // remove
    typedef Faunus::ParticleVector Tpvec;
    typedef Group<Particle> Tgroup;
    typedef std::vector<Tgroup> Tgvec;
    typedef Change Tchange;
    typedef std::function<void(Space &, double, double)> ScaleVolumeTrigger;
    typedef std::function<void(Space &, const Tchange &)> ChangeTrigger;
    typedef std::function<void(Space &, const Space &, const Tchange &)> SyncTrigger;

  private:
    /**
     * @brief Stores implicit molecules
     *
     * This stores the number (map value) of each implicit
     * molecule in the system, identified by the molecular id (map key)
     * The reservoir is used to emulate a canonical system where a finite
     * number of implicit molecules can participate in equilibrium reactions.
     */
    std::map<int, int> implicit_reservoir;
    std::vector<ChangeTrigger> changeTriggers; //!< Call when a Change object is applied (unused)
    std::vector<SyncTrigger> onSyncTriggers;   //!< Call when two Space objects are synched (unused)

  public:
    ParticleVector p;                                       //!< Particle vector storing all particles in system
    Tgvec groups;                                           //!< Group vector storing all molecules in system
    Tgeometry geo;                                          //!< Container geometry (boundaries, shape, volume)
    std::vector<ScaleVolumeTrigger> scaleVolumeTriggers;    //!< Called whenever the volume is scaled
    const std::map<int, int> &getImplicitReservoir() const; //!< Storage for implicit molecules
    std::map<int, int> &getImplicitReservoir();             //!< Storage for implicit molecules

    //!< Keywords to select particles based on the their active/inactive state and charge neutrality
    enum Selection { ALL, ACTIVE, INACTIVE, ALL_NEUTRAL, ACTIVE_NEUTRAL, INACTIVE_NEUTRAL };

    void clear();                                           //!< Clears particle and molecule list
    void push_back(int, const ParticleVector &);            //!< Safely add particles and corresponding group to back
    Tgvec::iterator findGroupContaining(const Particle &i); //!< Finds the group containing the given atom
    Tgvec::iterator findGroupContaining(size_t atom_index); //!< Finds the group containing given atom index
    size_t numParticles(Selection selection = ACTIVE) const; //!< Number of particles, all or active (default)
    Point scaleVolume(double, Geometry::VolumeMethod = Geometry::ISOTROPIC); //!< Scales atoms, molecules, container
    Tgvec::iterator randomMolecule(int, Random &, Selection = ACTIVE);       //!< Random group matching molid
    json info();

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

        assert(destination >= p.begin() && destination < p.end());
        assert(size <= std::distance(destination, p.end()));

        auto affected_groups = groups | ranges::cpp20::views::filter([=](auto &group) {
                                   return (group.begin() < destination + size) && (group.end() > destination);
                               }); // filtered group with affected groups only. Note we copy in org. `destination`

        // copy data from source range (this modifies `destination`)
        std::for_each(begin, end, [&](const auto &source) { copy_function(source, *destination++); });

        for (auto &group : affected_groups) { // update affected mass centers
            if (!group.empty()) {
                group.updateMassCenter(geo.getBoundaryFunc(), group.begin()->pos);
            }
        };
    }

    //! Iterable range of all particle positions
    auto positions() const {
        return ranges::cpp20::views::transform(p, [](auto &i) -> const Point & { return i.pos; });
    }

    //! Mutable iterable range of all particle positions
    auto positions() {
        return ranges::cpp20::views::transform(p, [](auto &i) -> Point & { return i.pos; });
    }

    /**
     * @brief Finds all groups of type `molid` (complexity: order N)
     * @param molid Molecular id to look for
     * @param sel Selection
     * @return range with all groups of molid
     */
    auto findMolecules(int molid, Selection sel = ACTIVE) {
        std::function<bool(Tgroup &)> f;
        switch (sel) {
        case (ALL):
            f = [molid](Tgroup &i) { return i.id == molid; };
            break;
        case (INACTIVE):
            f = [molid](Tgroup &i) { return (i.id == molid) && (i.size() != i.capacity()); };
            break;
        case (ACTIVE):
            f = [molid](Tgroup &i) { return (i.id == molid) && (i.size() == i.capacity()); };
            break;
        case (ALL_NEUTRAL):
            f = [molid](Tgroup &group) {
                if (group.id != molid)
                    return false;
                else {
                    double charge = std::accumulate(group.begin(), group.trueend(), 0.0,
                                                    [](double sum, auto &particle) { return sum + particle.charge; });
                    return (std::fabs(charge) < 1e-6);
                }
            };
            break;
        case (INACTIVE_NEUTRAL):
            f = [molid](Tgroup &group) {
                if ((group.id == molid) && (group.size() != group.capacity())) {
                    double charge = 0.0;
                    for (auto it = group.begin(); it != group.trueend(); ++it) {
                        charge += it->charge;
                    }
                    return (std::fabs(charge) <= pc::epsilon_dbl);
                } else {
                    return false;
                }
            };
            break;
        case (ACTIVE_NEUTRAL):
            f = [molid](Tgroup &group) {
                if ((group.id == molid) && (group.size() == group.capacity())) {
                    double charge = 0.0;
                    for (auto it = group.begin(); it != group.end(); ++it) {
                        charge += it->charge;
                    }
                    return (std::fabs(charge) <= pc::epsilon_dbl);
                } else {
                    return false;
                }
            };
            break;
        }
        return groups | ranges::cpp20::views::filter(f);
    }

    auto findAtoms(int atomid) {
        return p | ranges::cpp20::views::filter([&, atomid](const Particle &i) {
                   if (i.id == atomid)
                       for (const auto &g : groups)
                           if (g.contains(i, false))
                               return true;
                   return false;
               });
    } //!< Range with all active atoms of type `atomid` (complexity: order N)

    auto activeParticles() {
        return groups | ranges::cpp20::views::join;
    } //!< Returns range with all *active* particles in space

    /**
     * @brief Count number of molecules matching criteria
     * @param molid Molecule id to match
     * @tparam mask Selection mask based on `Group::Selectors`
     * @return Number of molecules matching molid and mask
     */
    template <unsigned int mask> auto numMolecules(int molid) const {
        auto filter = [&](const Tgroup &g) { return (g.id == molid) ? g.template match<mask>() : false; };
        return std::count_if(groups.begin(), groups.end(), filter);
    }

    void sync(const Space &other,
              const Tchange &change); //!< Copy differing data from other (o) Space using Change object

}; // end of space

void to_json(json &j, Space &spc);         //!< Serialize Space to json object
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
    static void insertAtomicGroups(MoleculeData &, Space &, int num_molecules, int num_inactive_molecules);
    static void insertMolecularGroups(MoleculeData &, Space &, int num_molecules, int num_inactive);
    static void setPositionsForTrailingGroups(Space &, int, const Faunus::ParticleVector &, const Point &);
    static void insertImplicitGroups(const MoleculeData &, Space &, int);

    //! Get number of molecules to insert from json object
    static int getNumberOfMolecules(const json &j, double volume, const std::string &molecule_name);

    //! Get number of inactive molecules from json object
    static int getNumberOfInactiveMolecules(const json &j, int number_of_molecules);

    //!< Get position vector using json object
    static ParticleVector getExternalPositions(const json &j, const std::string &molname);

    //!< Aggregated version of the above, called on each item in json array
    static void insertItem(const std::string &molname, const json &properties, Space &spc);

  public:
    static void insertMolecules(const json &, Space &);
}; // end of insertMolecules class

/**
 * @brief Helper class for range-based for-loops over *active* particles
 *
 * This class is currently not used as `Space::activeParticles` achieves the
 * same with much reduced code. However, this approach is expected to be
 * more efficient and is left here as an example of implementing a custom,
 * non-linear iterator.
 */
struct getActiveParticles {
    const Space &spc;
    class const_iterator {
      private:
        typedef typename Space::Tgvec::const_iterator Tgroups_iter;
        typedef typename Space::Tpvec::const_iterator Tparticle_iter;
        const Space &spc;
        Tparticle_iter particle_iter;
        Tgroups_iter groups_iter;

      public:
        const_iterator(const Space &spc, Tparticle_iter it);
        const_iterator operator++();
        bool operator!=(const const_iterator &other) const { return particle_iter != other.particle_iter; }
        auto operator*() const { return *particle_iter; }
    }; // enable range-based for loops

    const_iterator begin() const;
    const_iterator end() const;
    size_t size() const;
    getActiveParticles(const Space &spc);
};

// Make a global alias to easy transition to non-templated code
using Tspace = Space;

namespace SpaceFactory {

void makeNaCl(Space &, int, const Geometry::Chameleon &); //!< Create a simple salt system
void makeWater(Space &space, int num_particles, const Geometry::Chameleon &geometry);

} // namespace SpaceFactory

} // namespace Faunus
