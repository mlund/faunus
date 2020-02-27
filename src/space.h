#pragma once
#include "core.h"
#include "geometry.h"
#include "group.h"
#include "molecule.h"

namespace Faunus {

struct reservoir {
    std::string name;
    int N_reservoir;
    bool canonic;
};

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
    // bool chargeMove = false; //!< Not in use...

    struct data {
        bool dNatomic = false;  //!< True if the number of atomic molecules has changed
        bool dNswap = false;    //!< True if the number of atoms has changed as a result of a swap move
        int index;              //!< Touched group index
        bool internal = false;  //!< True if the internal energy/config has changed
        bool all = false;       //!< True if all particles in group have been updated
        std::vector<int> atoms; //!< Touched atom index w. respect to `Group::begin()`

        bool operator<(const data &a) const;
    }; //!< Properties of changed groups

    std::vector<data> groups; //!< Touched groups by index in group vector

    inline auto touchedGroupIndex() {
        return ranges::cpp20::views::transform(groups, [](data &i) -> int { return i.index; });
    } //!< List of moved groups (index)

    /** List of changed atom index relative to first particle in system) */
    std::vector<int> touchedParticleIndex(const std::vector<Group<Particle>> &);

    void clear();       //!< Clear all change data
    bool empty() const; //!< Check if change object is empty
    explicit operator bool() const;
    bool sanityCheck(Space &) const; //!< Performs a sanity check on contained object data
};

void to_json(json &, const Change::data &); //!< Serialize Change data to json
void to_json(json &, const Change &);       //!< Serialise Change object to json

/**
 * @brief Placeholder for atoms and molecules
 * @tparam Tparticletype Particle type for the space
 *
 * @todo
 * Split definitions into cpp file and forward instatiate
 * an instance of `Space<Tparticle>` where `Tparticle` has
 * been defined as `typedef Particle<Charge> Tparticle`.
 */
class Space {
  private:
    /**
     * Storage for the reservoir size (map value) of each implicit
     * molecule in the system, identified by the molecular id (map key)
     * The reservoir is used to emulate a canonical system where a finite
     * number of implicit molecules can participate in equilibrium reactions.
     */
    std::map<int, int> implicit_reservoir;

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

    std::vector<ScaleVolumeTrigger> scaleVolumeTriggers; //!< Call when volume is scaled
    std::vector<ChangeTrigger> changeTriggers;           //!< Call when a Change object is applied
    std::vector<SyncTrigger> onSyncTriggers;             //!< Call when two Space objects are synched

    Tpvec p;       //!< Particle vector
    Tgvec groups;  //!< Group vector
    Tgeometry geo; //!< Container geometry // TODO as a dependency injection in the constructor

    const std::map<int, int> &getImplicitReservoir() const; //!< Map of implicit molecule reservoirs
    std::map<int, int> &getImplicitReservoir();             //!< Map of implicit molecule reservoirs

    auto positions() const {
        return ranges::cpp20::views::transform(p, [](auto &i) -> const Point & { return i.pos; });
    } //!< Iterable range with positions

    //!< Keywords to select particles based on the their active/inactive state and charge neutrality
    enum Selection { ALL, ACTIVE, INACTIVE, ALL_NEUTRAL, ACTIVE_NEUTRAL, INACTIVE_NEUTRAL  };

    void clear(); //!< Clears particle and molecule list

    /*
     * The following is considered:
     *
     * - `groups` vector is expanded with a new group at the end
     * - if the particle vector is relocated, all existing group
     *   iterators are updated to reflect the new memory positions
     */
    void push_back(int molid, const Tpvec &in); //!< Safely add particles and corresponding group to back

    /**
     * @brief Find all groups of type `molid` (complexity: order N)
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
            f = [molid](Tgroup &i) {
                if (i.id != molid)
                    return false;
                else {
                    int charge = 0;
                    for (auto p = i.begin(); p != i.trueend(); ++p)
                        charge += p->charge;
                    return (charge == 0);
                }
            };
            break;
        case (INACTIVE_NEUTRAL):
            f = [molid](Tgroup &i) {
                if ((i.id != molid) && (i.size() == i.capacity()))
                    return false;
                else {
                    int charge = 0;
                    for (auto p = i.begin(); p != i.trueend(); ++p)
                        charge += p->charge;
                    return (charge == 0);
                }
            };
            break;
        case (ACTIVE_NEUTRAL):
            f = [molid](Tgroup &i) {
                if ((i.id != molid) && (i.size() != i.capacity()))
                    return false;
                else {
                    int charge = 0;
                    for (auto p = i.begin(); p != i.trueend(); ++p)
                        charge += p->charge;
                    return (charge == 0);
                }
            };
            break;

        }
        return groups | ranges::cpp20::views::filter(f);
    }

    typename Tgvec::iterator randomMolecule(int molid, Random &rand,
                                            Selection sel = ACTIVE); //!< Random group; groups.end() if not found

    auto findAtoms(int atomid) {
        return p | ranges::cpp20::views::filter([&, atomid](const Particle &i) {
                   if (i.id == atomid)
                       for (const auto &g : groups)
                           if (g.contains(i, false))
                               return true;
                   return false;
               });
    } //!< Range with all active atoms of type `atomid` (complexity: order N)

    auto findGroupContaining(const Particle &i) {
        return std::find_if(groups.begin(), groups.end(), [&i](auto &g) { return g.contains(i); });
    } //!< Finds the groups containing the given atom

    auto findGroupContaining(size_t index) {
        assert(index < p.size());
        return std::find_if(groups.begin(), groups.end(),
                            [&](auto &g) { return index < std::distance(p.begin(), g.end()); });
    } //!< Finds group containing given atom index

    auto activeParticles() {
        auto f = [&groups = groups](Particle &i) {
            for (auto &g : groups)
                if (g.contains(i)) // true if particle is within active part
                    return true;
            return false;
        };
        return p | ranges::cpp20::views::filter(f);
    } //!< Returns range with all *active* particles in space

    size_t numParticles(Selection sel = ACTIVE) const {
        size_t n = 0;
        if (sel == ALL)
            n = p.size();
        else if (sel == ACTIVE)
            for (auto &g : groups)
                n += g.size();
        else
            throw std::runtime_error("invalid selection");
        return n;
    } //!< Number of particles, all or active (default)

    void sync(Space &other, const Tchange &change); //!< Copy differing data from other (o) Space using Change object

    /*
     * Scales:
     * - positions of free atoms
     * - positions of molecular masscenters
     * - simulation container
     */
    void scaleVolume(double Vnew, Geometry::VolumeMethod method = Geometry::ISOTROPIC); //!< scale space to new volume

    json info();

}; // end of space

void to_json(json &j, Space &spc); //!< Serialize Space to json object

void from_json(const json &j, Space &spc); //!< Deserialize json object to Space

/**
 * @brief Insert molecules into space
 *
 * Expects an array of molecules, for example:
 *
 * ~~~ yaml
 *     - salt:  { N: 10 }
 *     - water: { N: 256 }
 *     - water: { N: 1, inactive: true }
 * ~~~
 */
void insertMolecules(const json &j, Space &spc); //!< Insert `N` molecules into space as defined in `insert`

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

} // namespace SpaceFactory

} // namespace Faunus
