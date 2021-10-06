#include "space.h"
#include "io.h"
#include "aux/iteratorsupport.h"
#include "spdlog/spdlog.h"
#include "aux/eigensupport.h"
#include <range/v3/algorithm/for_each.hpp>
#include <memory>
#include <stdexcept>

namespace Faunus {

bool Change::GroupChange::operator<(const Faunus::Change::GroupChange& other) const {
    return group_index < other.group_index;
}

void Change::clear() {
    *this = Change();
    assert(empty());
}
bool Change::empty() const {
    if (volume_change || everything || matter_change || !groups.empty()) {
        return false;
    }
    return true;
}
Change::operator bool() const { return !empty(); }

std::vector<Change::index_type> Change::touchedParticleIndex(const std::vector<Group<Particle>>& group_vector) const {
    std::vector<index_type> indices;                 // atom index rel. to first particle in system
    auto begin_first = group_vector.front().begin(); // first particle, first group
    for (const auto& changed : groups) {             // loop over changed groups
        auto begin_current = group_vector.at(changed.group_index).begin(); // first particle, current group
        auto offset = std::distance(begin_first, begin_current);           // abs. distance from first particle
        if (offset < 0) {
            throw std::runtime_error("negative index");
        }
        const auto offset_i = static_cast<index_type>(offset);
        indices.reserve(indices.size() + changed.relative_atom_indices.size());
        for (const auto index : changed.relative_atom_indices) { // atom index relative to group
            indices.push_back(index + offset_i);                 // atom index relative to first
        }
    }
    return indices;
}

/**
 * @param group_vector Vector of group connected to the change; typically `Space::groups`.
 * @throw If the atoms in the change object is outside range of given group index.
 */
void Change::sanityCheck(const std::vector<Group<Particle>>& group_vector) const {
    const auto first_particle = group_vector.at(0).begin(); // first particle in first group
    for (const auto& changed : groups) {
        const auto& group = group_vector.at(changed.group_index);
        for (auto index : changed.relative_atom_indices) {
            if (index >= group.capacity()) {
                auto first = std::distance(first_particle, group.begin());
                auto last = std::distance(first_particle, group.trueend()) - 1;
                faunus_logger->error("atom {} outside group capacity: '{}' ({}-{})", index + first,
                                     molecules.at(group.id).name, first, last);
                throw std::runtime_error("insane change object: atom outside group capacity");
            }
        }
    }
}

void to_json(json& j, const Change::GroupChange& group_change) {
    j = {{"all", group_change.all},           {"internal", group_change.internal},
         {"dNswap", group_change.dNswap},     {"dNatomic", group_change.dNatomic},
         {"index", group_change.group_index}, {"atoms", group_change.relative_atom_indices}};
}

void to_json(json& j, const Change& change) {
    j = {{"dV", change.volume_change},
         {"all", change.everything},
         {"dN", change.matter_change},
         {"moved2moved", change.moved_to_moved_interactions},
         {"groups", change.groups}};
}

TEST_CASE("[Faunus] Change") {
    Change change;
    CHECK(not change);
    change.volume_change = true;
    CHECK(not change.empty());
    CHECK(change);
    change.clear();
    CHECK(change.empty());
}

void Space::clear() {
    particles.clear();
    groups.clear();
    implicit_reservoir.clear();
}

/**
 * The following is considered:
 *
 * - `groups` vector is expanded with a new group at the end
 * - if the space particle vector is relocated, all existing group
 *   iterators are updated to reflect the new memory positions
 * - for molecular groups, the mass center is calculated and set
 *
 * @param molid Molecule id of inserted data
 * @param new_particles Particles to generate group from
 * @return Reference to the inserted group
 * @throw if particles is empty or if the given particles do not match the molecule id
 */
Space::GroupType& Space::addGroup(MoleculeData::index_type molid, const ParticleVector& new_particles) {
    if (new_particles.empty()) {
        throw std::runtime_error("cannot add empty molecule");
    }
    auto original_begin = particles.begin(); // used to detect if `particles` is relocated
    particles.insert(particles.end(), new_particles.begin(), new_particles.end()); // insert particle into space
    if (particles.begin() != original_begin) { // update group iterators if `particles` is relocated
        std::for_each(groups.begin(), groups.end(),
                      [&](auto& group) { group.relocate(original_begin, particles.begin()); });
    }
    GroupType group(particles.end() - new_particles.size(), particles.end()); // create a group
    const auto& moldata = Faunus::molecules.at(molid);
    group.id = molid;
    group.compressible = moldata.compressible;
    group.atomic = moldata.atomic;
    if (group.isMolecular()) {
        if (new_particles.size() != moldata.atoms.size()) {
            faunus_logger->error("{} requires {} atoms but {} were provided", moldata.name, moldata.atoms.size(),
                                 new_particles.size());
            throw std::runtime_error("particle size mismatch");
        }
        group.cm =
            Geometry::massCenter(group.begin(), group.end(), geometry.getBoundaryFunc(), -new_particles.begin()->pos);
    } else {
        if (new_particles.size() % moldata.atoms.size() != 0) {
            throw std::runtime_error("indivisible by atomic group size: "s + moldata.name);
        }
    }
    return groups.emplace_back(group);
}

/**
 * @param other Space to copy from
 * @param change Change object describing the changes beteeen the two Space objects
 *
 * Copy data from another Space according to Change object. This is typically done
 * after a Monte Carlo move has either been accepted or rejected.
 * The `other` space *must* be populated in the exact same way, i.e. must have the same
 * molecules and particles. In DEBUG mode, several assertions are included to ensure this is true.
 *
 * Copied data includes:
 *
 * - geometry
 * - groups
 * - particles
 * - implicit molecules
 */
void Space::sync(const Space &other, const Change &change) {
    if (&other != this && !change.empty()) {
        assert(!groups.empty());
        assert(particles.size() == other.particles.size());
        assert(groups.size() == other.groups.size());
        assert(implicit_reservoir.size() == other.implicit_reservoir.size());
        if (change.volume_change or change.everything) {
            geometry = other.geometry; // copy simulation geometry
        }
        if (change.everything) {                                            // deep copy *everything*
            implicit_reservoir = other.implicit_reservoir;                  // copy all implicit molecules
            particles = other.particles;                                    // copy all positions
            groups = other.groups;                                          // copy all groups
            assert(particles.begin() != other.particles.begin());           // check deep copy problem
            assert(groups.front().begin() != other.groups.front().begin()); // check deep copy problem
        } else {
            for (const auto &changed : change.groups) {             // look over changed groups
                auto& group = groups.at(changed.group_index);       // old group
                auto& other_group = other.groups.at(changed.group_index); // new group
                assert(group.id == other_group.id);
                if (group.traits().isImplicit()) { // the molecule is implicit
                    implicit_reservoir[group.id] = other.implicit_reservoir.at(group.id);
                } else if (changed.all) {
                    group = other_group;            // copy everything
                } else {                            // copy only a subset
                    group.shallowcopy(other_group); // copy group data but *not* particles
                    for (auto i : changed.relative_atom_indices) { // loop over atom index (rel. to group)
                        group.at(i) = other_group.at(i); // deep copy select particles
                    }
                }
            }
        }
    }
}

/**
 * @param Vnew New volume
 * @param method Scaling policy
 * @returns Scaling factors in each dimension
 * @warning Check new_volume/old_volume for ISOCHORIC in case of external triggers (see end of function)
 */
Point Space::scaleVolume(double Vnew, Geometry::VolumeMethod method) {
    for (auto &group : groups) {             // remove PBC on molecular groups ...
        group.unwrap(geometry.getDistanceFunc()); // ... before scaling the volume
    }
    auto Vold = geometry.getVolume();               // simulation volume before move
    Point scale = geometry.setVolume(Vnew, method); // scale volume of simulation container

    auto scale_position = [&](auto &particle) {
        particle.pos = particle.pos.cwiseProduct(scale);
        geometry.boundary(particle.pos);
    }; //!< unary helper function to scale position and apply PBC

    for (auto &group : groups) { // loop over all molecules in system
        if (group.empty()) {
            continue;
        } else if (group.isAtomic()) { // scale all particle positions
            std::for_each(group.begin(), group.end(), scale_position);
        } else {
            auto original_mass_center = group.cm;
            if (group.traits().compressible) { // scale positions; recalculate mass center
                std::for_each(group.begin(), group.end(), scale_position);
                group.cm =
                    Geometry::massCenter(group.begin(), group.end(), geometry.getBoundaryFunc(), -original_mass_center);
            } else { // scale mass center; translate positions
                group.cm = group.cm.cwiseProduct(scale);
                geometry.boundary(group.cm);
                auto mass_center_displacement = group.cm - original_mass_center;
                std::for_each(group.begin(), group.end(), [&](auto &particle) {
                    particle.pos += mass_center_displacement; // translate internal coordinates
                    geometry.boundary(particle.pos);          // apply PBC
                });
#ifndef NDEBUG
                auto recalc_cm = Geometry::massCenter(group.begin(), group.end(), geo.getBoundaryFunc(), -group.cm);
                if (double error = geo.sqdist(group.cm, recalc_cm); error > 1e-6) {
                    assert(false); // mass center mismatch
                }
#endif
            }
        }
    }
    if (method == Geometry::VolumeMethod::ISOCHORIC) { // ? not used for anything...
        Vold = std::pow(Vold, 1.0 / 3.0);              // ?
    }
    for (auto trigger_function : scaleVolumeTriggers) { // external clients may have added function
        trigger_function(*this, Vold, Vnew);            // to be triggered upon each volume change
    }
    return scale;
}

json Space::info() {
    json j = {{"number of particles", particles.size()}, {"number of groups", groups.size()}, {"geometry", geometry}};
    auto &j_groups = j["groups"];
    for (const auto& group : groups) {
        auto &molname = Faunus::molecules.at(group.id).name;
        json tmp, d = group;
        d.erase("cm");
        d.erase("id");
        d.erase("atomic");
        auto ndx = group.to_index(particles.begin()); // absolute index
        if (not group.empty()) {
            d["index"] = {ndx.first, ndx.second};
        }
        tmp[molname] = d;
        j_groups.push_back(tmp);
    }
    auto &j_reactionlist = j["reactionlist"];
    for (auto &reaction : Faunus::reactions) {
        j_reactionlist.push_back(reaction);
    }
    return j;
}

Space::GroupVector::iterator Space::randomMolecule(MoleculeData::index_type molid, Random& rand,
                                                   Space::Selection selection) {
    auto found_molecules = findMolecules(molid, selection);
    if (ranges::cpp20::empty(found_molecules)) {
        return groups.end();
    }
    return groups.begin() + (&*rand.sample(found_molecules.begin(), found_molecules.end()) - &*groups.begin());
}

const std::map<MoleculeData::index_type, std::size_t>& Space::getImplicitReservoir() const {
    return implicit_reservoir;
}

std::map<MoleculeData::index_type, std::size_t>& Space::getImplicitReservoir() { return implicit_reservoir; }

std::vector<Space::GroupType, std::allocator<Space::GroupType>>::iterator
Space::findGroupContaining(const Particle& particle) {
    return std::find_if(groups.begin(), groups.end(), [&particle](auto& group) { return group.contains(particle); });
}

std::vector<Space::GroupType, std::allocator<Space::GroupType>>::iterator
Space::findGroupContaining(AtomData::index_type atom_index) {
    assert(atom_index < particles.size());
    return std::find_if(groups.begin(), groups.end(),
                        [&](auto& g) { return atom_index < std::distance(particles.begin(), g.end()); });
}

std::size_t Space::numParticles(Space::Selection selection) const {
    if (selection == Selection::ALL) {
        return particles.size();
    }
    if (selection == Selection::ACTIVE) {
        return std::accumulate(groups.begin(), groups.end(), 0u,
                               [](auto sum, auto& group) { return sum + group.size(); });
    }
    throw std::runtime_error("invalid selection");
}

/**
 * @throw If group is not part of space
 */
std::size_t Space::getGroupIndex(const Space::GroupType& group) const {
    const auto distance = std::addressof(group) - std::addressof(groups.front()); // std::ptrdiff_t
    if (distance >= 0) {
        const auto index = static_cast<std::size_t>(distance);
        if (index < groups.size()) {
            return index;
        }
    }
    throw std::out_of_range("invalid group index");
}

std::size_t Space::getFirstParticleIndex(const GroupType& group) const {
    if (group.capacity() > 0 && !particles.empty()) {
        const auto distance = std::distance<ParticleVector::const_iterator>(particles.cbegin(), group.begin());
        if (distance >= 0) {
            const auto particle_index = static_cast<std::size_t>(distance);
            if (particle_index < particles.size()) {
                return particle_index;
            }
        }
    }
    throw std::overflow_error("group empty or not part of space");
}

/**
 * Returns the index of the first particle of the group in the range returned by `activeParticles()`
 */
std::size_t Space::getFirstActiveParticleIndex(const GroupType& group) const {
    const auto group_index = getGroupIndex(group);
    std::size_t index = 0u;
    std::for_each(groups.begin(), groups.begin() + group_index, [&](const auto& g) { index += g.size(); });
    return index;
}

TEST_CASE("Space::numParticles") {
    Space spc;
    spc.particles.resize(10);
    CHECK(spc.particles.size() == spc.numParticles(Space::Selection::ALL));
    CHECK(spc.numParticles(Space::Selection::ACTIVE) == 0);  // zero as there are still no groups
    spc.groups.emplace_back(spc.particles.begin(), spc.particles.end() - 2); // enclose first 8 particles in group
    CHECK(spc.numParticles(Space::Selection::ACTIVE) == 8);
    spc.groups.emplace_back(spc.particles.end() - 2, spc.particles.end()); // enclose last 2 particles in group
    CHECK(spc.numParticles(Space::Selection::ACTIVE) == 10);
    spc.groups.front().resize(0); // deactivate first group with 8 particles
    CHECK(spc.numParticles(Space::Selection::ACTIVE) == 2);
}

void to_json(json& j, const Space& spc) {
    j["geometry"] = spc.geometry;
    j["groups"] = spc.groups;
    j["particles"] = spc.particles;
    j["reactionlist"] = reactions;
    j["implicit_reservoir"] = spc.getImplicitReservoir();
}
void from_json(const json &j, Space &spc) {
    try {
        if (atoms.empty()) {
            atoms = j.at("atomlist").get<decltype(atoms)>();
        }
        if (molecules.empty()) {
            molecules = j.at("moleculelist").get<decltype(molecules)>();
        }
        if (reactions.empty()) {
            if (j.contains("reactionlist")) {
                reactions = j.at("reactionlist").get<decltype(reactions)>();
            }
        }

        spc.clear();
        spc.geometry = j.at("geometry");

        if (!j.contains("groups")) {
            InsertMoleculesInSpace::insertMolecules(j.at("insertmolecules"), spc);
        } else {
            spc.particles = j.at("particles").get<ParticleVector>();
            if (!spc.particles.empty()) {
                auto begin = spc.particles.begin();
                Space::GroupType g(begin, begin); // create new grou[
                for (auto &i : j.at("groups")) {
                    g.begin() = begin;
                    from_json(i, g);
                    spc.groups.push_back(g);
                    begin = g.trueend();
                }
                if (begin != spc.particles.end()) {
                    throw ConfigurationError("load error");
                }
            }
        }

        if (auto it = j.find("implicit_reservoir"); it != j.end() && it->is_array()) {
            for (auto vec : *it) {
                assert(vec.is_array() && vec.size() == 2);
                spc.getImplicitReservoir()[vec[0]] = vec[1];
            }
            faunus_logger->trace("{} implicit molecules loaded", it->size());
        }

        // check correctness of molecular mass centers
        auto check_mass_center = [&](const auto& group) {
            const auto should_be_small = spc.geometry.sqdist(
                group.cm, Geometry::massCenter(group.begin(), group.end(), spc.geometry.getBoundaryFunc(), -group.cm));
            if (should_be_small > 1e-9) {
                throw std::runtime_error(fmt::format(
                    "couldn't calculate mass center for {}; increase periodic box size?", group.traits().name));
            }
        };
        auto active_and_molecular = [](const auto& group) { return (!group.empty() && group.isMolecular()); };
        ranges::cpp20::for_each(spc.groups | ranges::cpp20::views::filter(active_and_molecular), check_mass_center);
    } catch (const std::exception& e) { throw std::runtime_error("error building space -> "s + e.what()); }
}

TEST_SUITE_BEGIN("Space");

TEST_CASE("[Faunus] Space") {
    using doctest::Approx;
    Space spc1;
    spc1.geometry = R"( {"type": "sphere", "radius": 1e9} )"_json;

    // check molecule insertion
    Faunus::molecules.at(0).atomic = false;
    REQUIRE_EQ(Faunus::molecules.at(0).atomic, false);
    Faunus::molecules.at(0).atoms.resize(2);
    CHECK(Faunus::molecules.at(0).atoms.size() == 2);

    Faunus::atoms.resize(2);
    CHECK(Faunus::atoms.at(0).mw == 1);

    Particle a;
    a.id = 0;
    a.pos.setZero();
    ParticleVector p(2, a);
    CHECK(p[0].traits().mw == 1);
    p[0].pos.x() = 2;
    p[1].pos.x() = 3;
    spc1.addGroup(0, p); // insert molecular group
    CHECK(spc1.particles.size() == 2);
    CHECK(spc1.groups.size() == 1);
    CHECK(spc1.groups.back().isMolecular());
    CHECK(spc1.groups.front().id == 0);
    CHECK(spc1.groups.front().cm.x() == doctest::Approx(2.5));

    // check `positions()`
    CHECK(&spc1.positions()[0] == &spc1.particles[0].pos);

    SUBCASE("ActiveParticles") {
        // add three groups to space
        Space spc;
        spc.geometry = R"( {"type": "sphere", "radius": 1e9} )"_json;
        Particle a;
        a.pos.setZero();
        a.id = 0;
        ParticleVector pvec({a, a, a});

        Faunus::molecules.at(0).atoms.resize(3);

        spc.addGroup(0, pvec);
        spc.addGroup(0, pvec);
        spc.addGroup(0, pvec);

        spc.groups.at(0).id = 0;
        spc.groups.at(1).id = 1;
        spc.groups.at(2).id = 0;

        for (size_t i = 0; i < spc.particles.size(); i++)
            spc.particles[i].charge = double(i);

        CHECK(spc.particles.size() == 9);
        CHECK(spc.groups.size() == 3);

        CHECK(spc.numMolecules<Space::GroupType::ANY>(0) == 2);
        CHECK(spc.numMolecules<Space::GroupType::ANY>(1) == 1);

        spc.groups[0].deactivate(spc.particles.begin(), spc.particles.begin() + 1);
        spc.groups[1].deactivate(spc.groups[1].begin(), spc.groups[1].end());

        CHECK(spc.groups[0].size() == 2);
        CHECK(spc.groups[1].size() == 0);
        CHECK(spc.groups[2].size() == 3);

        CHECK(spc.numMolecules<Space::GroupType::ACTIVE>(0) == 2);
        CHECK(spc.numMolecules<Space::GroupType::ACTIVE | Space::GroupType::NEUTRAL>(0) == 0);
        CHECK(spc.numMolecules<Space::GroupType::ACTIVE>(1) == 0);
        CHECK(spc.numMolecules<Space::GroupType::INACTIVE>(1) == 1);
        CHECK(spc.numMolecules<Space::GroupType::INACTIVE | Space::GroupType::NEUTRAL>(1) == 0);

        // check the rangev3 implementation in `activeParticles()`:
        auto p2 = spc.activeParticles();
        CHECK(std::distance(p2.begin(), p2.end()) == 5);
        std::vector<int> vals;
        for (const auto& particle : p2) {
            vals.push_back(static_cast<int>(particle.charge));
        }
        CHECK(vals == std::vector<int>({1, 2, 6, 7, 8}));
    }
}

TEST_CASE("[Faunus] Space::updateParticles") {
    using doctest::Approx;
    Space spc;
    Geometry::Cuboid geo({100, 100, 100});
    spc.geometry = Geometry::Chameleon(geo, Geometry::Variant::CUBOID);

    spc.particles.resize(2);

    ParticleVector p(2);
    p[0].pos = {0, 0, 0};
    p[1].pos = {2, 0, 0};

    spc.updateParticles(p.begin(), p.end(), spc.particles.begin());
    CHECK(spc.particles[0].pos.x() == p[0].pos.x());
    CHECK(spc.particles[1].pos.x() == p[1].pos.x());

    std::vector<Point> positions = {{2.1, 0, 0}, {0.9, 0, 0}};
    spc.updateParticles(positions.begin(), positions.end(), spc.particles.begin(),
                        [](const auto& pos, auto& particle) { particle.pos = pos; });
    CHECK(spc.particles[0].pos.x() == 2.1);
    CHECK(spc.particles[1].pos.x() == 0.9);

    SUBCASE("Group update") {
        Space spc;
        SpaceFactory::makeWater(spc, 2, R"( {"type": "cuboid", "length": 20} )"_json);
        CHECK(spc.groups.size() == 2);

        auto copy_function = [](const auto &pos, auto &particle) { particle.pos = pos; };
        std::vector<Point> positions = {{0, 0, 0}, {3, 3, 3}, {6, 6, 6}};

        // first group affected
        spc.groups[0].cm.x() = -1;
        spc.groups[1].cm.x() = -1;
        spc.updateParticles(positions.begin(), positions.end(), spc.groups[1].begin(), copy_function);
        CHECK(spc.groups[0].cm.x() == Approx(-1));
        CHECK(spc.groups[1].cm.x() == Approx(0.5031366235));

        // second group affected
        spc.groups[1].cm.x() = -1;
        spc.updateParticles(positions.begin(), positions.end(), spc.groups[0].begin(), copy_function);
        CHECK(spc.groups[0].cm.x() == Approx(0.5031366235));
        CHECK(spc.groups[1].cm.x() == Approx(-1));

        // both groups affected
        spc.groups[0].cm.x() = -1;
        spc.groups[1].cm.x() = -1;
        spc.updateParticles(positions.begin(), positions.end(), spc.groups[0].begin() + 1, copy_function);
        CHECK(spc.groups[0].cm.x() == Approx(0.1677122078));
        CHECK(spc.groups[1].cm.x() == Approx(5.8322877922));
    }
}

TEST_SUITE_END();

namespace SpaceFactory {

/**
 * @param space Space to insert into (will be overwritten)
 * @param num_particles Number of salt pairs to insert
 * @param geometry Geometry to use
 *
 * Create a system with two atom types, "Na" and "Cl", forming
 * an atomic molecule, "salt". N salt pairs a randomly inserted
 */
void makeNaCl(Space& space, size_t num_particles, const Geometry::Chameleon& geometry) {
    pc::temperature = 298.15_K;
    space.geometry = geometry;

    Faunus::atoms = R"([
             { "Na": { "sigma": 3.8, "eps": 0.1, "q": 1.0 } },
             { "Cl": { "sigma": 4.0, "eps": 0.05, "q": -1.0 } }
             ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([
                { "salt": {"atomic": true, "atoms": ["Na", "Cl"] } }
            ])"_json.get<decltype(molecules)>();

    json j = json::array();
    j.push_back({{"salt", {{"N", num_particles}}}});
    InsertMoleculesInSpace::insertMolecules(j, space);
}

/**
 * @param space Space to insert into (will be overwritten)
 * @param num_particles Number of salt pairs to insert
 * @param geometry Geometry to use
 */
void makeWater(Space& space, size_t num_particles, const Geometry::Chameleon& geometry) {
    pc::temperature = 298.15_K;
    space.geometry = geometry;

    Faunus::atoms = R"([
             { "OW": { "sigma": 3.166, "eps": 0.65, "q": -0.8476, "mw": 15.999 } },
             { "HW": { "sigma": 2.0, "eps": 0.0, "q": 0.4238, "mw": 1.007 } }
             ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([{
             "water": {
                 "structure": [
                     {"OW": [2.3, 6.28, 1.13]},
                     {"HW": [1.37, 6.26, 1.5]},
                     {"HW": [2.31, 5.89, 0.21]}
                 ]
            }}])"_json.get<decltype(molecules)>();

    json j = json::array();
    j.push_back({{"water", {{"N", num_particles}}}});
    InsertMoleculesInSpace::insertMolecules(j, space);
}

} // namespace SpaceFactory

TEST_CASE("SpaceFactory") {
    Space spc;
    SpaceFactory::makeNaCl(spc, 10, R"( {"type": "cuboid", "length": 20} )"_json);
    CHECK(spc.numParticles() == 20);

    SUBCASE("makeWater") {
        Space spc;
        SpaceFactory::makeWater(spc, 2, R"( {"type": "cuboid", "length": 20} )"_json);
        CHECK(spc.numParticles() == 6);
    }
}

/**
 * @param moldata Molecule type to insert. Must be `atomic`.
 * @param spc Space to insert into (at the end)
 * @param num_molecules Number of molecules to insert into the *same* group
 * @param num_inactive_molecules Number of molecules to declare inactive
 * @throw If number if inactive is higher than number of active molecules
 *
 * The atoms in the atomic groups are repeated `num_molecules` times, then inserted
 * into Space.
 */
void InsertMoleculesInSpace::insertAtomicGroups(MoleculeData& moldata, Space& spc, size_t num_molecules,
                                                size_t num_inactive_molecules) {
    assert(moldata.atomic);
    ParticleVector repeated_particles;
    repeated_particles.reserve(num_molecules * moldata.atoms.size()); // prepare memory
    for (size_t i = 0; i < num_molecules; i++) {                      // repeat insertion into the same atomic group
        auto particles = moldata.getRandomConformation(spc.geometry, spc.particles);
        repeated_particles.insert(repeated_particles.end(), particles.begin(), particles.end());
    }

    spc.addGroup(moldata.id(), repeated_particles); // create new group in Space

    if (num_inactive_molecules > num_molecules) {
        throw std::runtime_error("too many inactive molecules requested");
    }
    if (num_inactive_molecules > 0) {
        auto num_active_atoms = (num_molecules - num_inactive_molecules) * moldata.atoms.size();
        spc.groups.back().resize(num_active_atoms); // deactivate atoms in molecule
    }
}

void InsertMoleculesInSpace::insertMolecularGroups(MoleculeData& moldata, Space& spc, size_t num_molecules,
                                                   size_t num_inactive) {
    assert(moldata.atomic == false);
    if (num_inactive > num_molecules) {
        throw std::runtime_error("too many inactive molecules requested");
    }
    for (size_t i = 0; i < num_molecules; i++) { // insert molecules
        spc.addGroup(moldata.id(), moldata.getRandomConformation(spc.geometry, spc.particles));
    }
    // deactivate groups, starting from the back
    std::for_each(spc.groups.rbegin(), spc.groups.rbegin() + num_inactive, [&](auto& group) {
        group.unwrap(spc.geometry.getDistanceFunc()); // make molecules whole (remove PBC) ...
        group.resize(0);                              // ... and then deactivate
    });
}

/**
 * @brief Set positions for num_molecules last groups in Space
 * @param spc Space to insert into
 * @param num_molecules num_molecules Number of groups set affect
 * @param particles Position vector for all particles in the N groups
 * @param offset Translate positions by this offset
 * @throw if num_molecules doesn't match, or if positions are outside simulation cell
 *
 * Sets particle positions in space using a given input
 * particle vector. The vector can span several *identical* groups.
 * The number of groups affected must be given in order to update their mass-centers.
 * Only *positions* are affected.
 */
void InsertMoleculesInSpace::setPositionsForTrailingGroups(Space& spc, size_t num_molecules,
                                                           const Faunus::ParticleVector& particles,
                                                           const Point& offset) {
    assert(spc.groups.size() >= num_molecules);
    if (particles.size() != num_molecules * (spc.groups.rbegin())->traits().atoms.size()) {
        throw std::runtime_error("number of particles doesn't match groups");
    }
    // update positions in space, starting from the back
    std::transform(particles.rbegin(), particles.rend(), spc.particles.rbegin(), spc.particles.rbegin(),
                   [&](const auto& src, auto& dst) {
                       dst.pos = src.pos + offset;            // shift by offset
                       if (spc.geometry.collision(dst.pos)) { // check if position is inside simulation volume
                           throw std::runtime_error("positions outside box");
                       }
                       return dst;
                   });
    // update mass-centers on modified groups; start from the back
    std::for_each(spc.groups.rbegin(), spc.groups.rbegin() + num_molecules, [&](auto& group) {
        group.cm =
            Geometry::massCenter(group.begin(), group.end(), spc.geometry.getBoundaryFunc(), -group.begin()->pos);
    });
}

/**
 * @brief Insert implicit groups into Space
 * @param moldata Molecule to insert. Must be implicit.
 * @param spc Space to insert into.
 * @param num_molecules Nunber of implicit molecules to insert
 */
void InsertMoleculesInSpace::insertImplicitGroups(const MoleculeData& moldata, Space& spc, size_t num_molecules) {
    assert(moldata.isImplicit());
    spc.getImplicitReservoir()[moldata.id()] = num_molecules;
}

/**
 * @param molname Molecule name
 * @param properties json object with insertion properties ('N', 'molarity', 'inactive' etc)
 * @param spc Space to insert into
 */
void InsertMoleculesInSpace::insertItem(const std::string &molname, const json &properties, Space &spc) {
    auto& moldata = findMoleculeByName(molname);
    const auto num_molecules = getNumberOfMolecules(properties, spc.geometry.getVolume(), molname);
    if (num_molecules == 0) {
        if (!moldata.isImplicit()) {
            throw ConfigurationError("one or more {} molecule(s) required; concentration too low?", molname);
        }
    } else {
        const auto num_inactive = getNumberOfInactiveMolecules(properties, num_molecules);
        const auto molarity = (num_molecules - num_inactive) / spc.geometry.getVolume() / 1.0_molar;
        if (moldata.isImplicit()) {
            faunus_logger->info("adding {} implicit {} molecules --> {:.5E} mol/l", num_molecules, molname, molarity);
            insertImplicitGroups(moldata, spc, num_molecules);
        } else {
            faunus_logger->info("adding {} {} molecules --> {:.5E} mol/l ({} inactive)", num_molecules, molname,
                                molarity, num_inactive);
            if (moldata.atomic) {
                insertAtomicGroups(moldata, spc, num_molecules, num_inactive);
            } else {
                insertMolecularGroups(moldata, spc, num_molecules, num_inactive);
                if (auto particles = getExternalPositions(properties, molname); !particles.empty()) {
                    auto offset = properties.value("translate", Point(0, 0, 0));
                    setPositionsForTrailingGroups(spc, num_molecules, particles, offset);
                }
            }
        }
    }
}

/**
 * Look for 'positions' in json object and load structure if found.
 *
 * @param j json object to search for key `positions`
 * @param molname name of molecule to affect; used for logging, only.
 * @returns Particle vector; empty if no external positions are requested
 * @throw If the 'positions' key is there, but could not be loaded
 */
ParticleVector InsertMoleculesInSpace::getExternalPositions(const json &j, const std::string &molname) {
    ParticleVector particles;
    if (auto filename = j.value("positions", ""s); !filename.empty()) {
        particles = loadStructure(filename, false); // throws on error
        faunus_logger->info("{}: loaded position file {}", molname, filename);
    }
    return particles;
}

/**
 * @brief Insert molecules into Space based on JSON input
 * @param j JSON array
 * @param spc Space to insert into
 */
void InsertMoleculesInSpace::insertMolecules(const json &j, Space &spc) {
    if (!j.is_array()) {
        throw ConfigurationError("molecules to insert must be an array");
    }
    spc.clear();
    for (const auto& item : j) { // loop over array of molecules
        const auto& [molecule_name, properties] = jsonSingleItem(item);
        try {
            insertItem(molecule_name, properties, spc);
        } catch (std::exception& e) { throw ConfigurationError("error inserting {}: {}", molecule_name, e.what()); }
    }
}

/**
 * @param j Input json object
 * @param volume Volume of simulation container used to calculate initial concentration
 * @param molecule_name Name of molecule needed for logging
 * @return Number of molecules to insert
 *
 * Looks for json key `N` or `molarity`. For the latter, the nearest corresponding
 * number of particles is calculated based on the given system volume.
 */
size_t InsertMoleculesInSpace::getNumberOfMolecules(const json& j, double volume, const std::string& molecule_name) {
    size_t num_molecules = 0;
    if (j.contains("N")) {
        num_molecules = j.at("N").get<size_t>();
    } else {
        auto density = j.at("molarity").get<double>() * 1.0_molar;
        num_molecules = std::round(density * volume);
        const double error_limit = 0.01; // warn if relative density error is above this
        if (auto error = (density - num_molecules / volume) / density; std::fabs(error) > error_limit) {
            faunus_logger->warn("{}: initial molarity differs by {}% from target value", molecule_name, error * 100);
        }
    }
    return num_molecules;
}

/**
 * @param j Input json object
 * @param number_of_molecules Total number of molecules
 * @return Number of molecules to be inactive (always smaller than `number_of_molecules`)
 * @throw If Inactive molecules is higher than `number_of_molecules`
 *
 * Looks for key "inactive" and if:
 * - boolean true: all molecules are inactive
 * - number: number of molecules to declare inactive
 * - no `inactive` key found, return zero
 */
size_t InsertMoleculesInSpace::getNumberOfInactiveMolecules(const json& j, size_t number_of_molecules) {
    size_t number_of_inactive_molecules = 0; // number of inactive molecules
    if (auto inactive = j.find("inactive"); inactive != j.end()) {
        if (inactive->is_boolean()) {
            if (inactive->get<bool>()) {
                number_of_inactive_molecules = number_of_molecules; // all molecules are inactive
            }
        } else if (inactive->is_number_integer()) {
            number_of_inactive_molecules = inactive->get<size_t>(); // a subset are inactive
            if (number_of_inactive_molecules > number_of_molecules) {
                throw ConfigurationError("too many inactive particles requested");
            }
        }
    }
    return number_of_inactive_molecules;
}
} // namespace Faunus
