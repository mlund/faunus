#include "space.h"
#include "io.h"
#include "aux/iteratorsupport.h"
#include "spdlog/spdlog.h"
#include "aux/eigensupport.h"
#include <memory>
#include <stdexcept>

namespace Faunus {

bool Change::data::operator<(const Faunus::Change::data &other) const { return index < other.index; }

void Change::clear() {
    *this = Change();
    assert(empty());
}
bool Change::empty() const {
    if (dV || all || dN || !groups.empty()) {
        return false;
    } else {
        return true;
    }
}
Change::operator bool() const { return !empty(); }

std::vector<int> Change::touchedParticleIndex(const std::vector<Group<Particle>> &group_vector) const {
    std::vector<int> atom_indexes;                                   // atom index rel. to first particle in system
    for (const auto &changed : groups) {                             // loop over changed groups
        auto begin_first = group_vector.front().begin();             // first particle, first group
        auto begin_current = group_vector.at(changed.index).begin(); // first particle, current group
        auto offset = std::distance(begin_first, begin_current);     // abs. distance from first particle
        atom_indexes.reserve(atom_indexes.size() + changed.atoms.size());
        for (auto index : changed.atoms) {          // atom index relative to group
            atom_indexes.push_back(index + offset); // atom index relative to first
        }
    }
    return atom_indexes;
}

/**
 * @param group_vector Vector of group connected to the change; typically `Space::groups`.
 * @throws If the atoms in the change object is outside range of given group index.
 */
void Change::sanityCheck(const std::vector<Group<Particle>> &group_vector) const {
    for (auto &changed : groups) {
        auto first_particle = group_vector.at(0).begin(); // first particle in first group
        auto &group = group_vector.at(changed.index);     // current group
        for (auto i : changed.atoms) {                    // all atoms must be within `group`
            if (i >= group.capacity()) {
                auto first = std::distance(first_particle, group.begin());
                auto last = std::distance(first_particle, group.trueend()) - 1;
                faunus_logger->error("atom {} outside group capacity: '{}' ({}-{})", i + first,
                                     molecules.at(group.id).name, first, last);
                throw std::runtime_error("insane change object: atom outside group capacity");
            }
        }
    }
}

void to_json(json &j, const Change::data &d) {
    j = {{"all", d.all},           {"internal", d.internal}, {"dNswap", d.dNswap},
         {"dNatomic", d.dNatomic}, {"index", d.index},       {"atoms", d.atoms}};
}

void to_json(json &j, const Change &c) {
    j = {{"dV", c.dV}, {"all", c.all}, {"dN", c.dN}, {"moved2moved", c.moved2moved}, {"groups", c.groups}};
}

TEST_CASE("[Faunus] Change") {
    Change change;
    CHECK(not change);
    change.dV = true;
    CHECK(not change.empty());
    CHECK(change);
    change.clear();
    CHECK(change.empty());
}

void Space::clear() {
    p.clear();
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
 * @param particles Particles to generate group from
 */
void Space::push_back(int molid, const ParticleVector &particles) {
    if (!particles.empty()) {
        auto original_begin = p.begin();                       // used to detect if `p` is relocated
        p.insert(p.end(), particles.begin(), particles.end()); // insert particle into space
        if (p.begin() != original_begin) {                     // update group iterators if `p` is relocated
            std::for_each(groups.begin(), groups.end(),
                          [&](auto &group) { group.relocate(original_begin, p.begin()); });
        }
        Tgroup group(p.end() - particles.size(), p.end()); // create a group
        const auto &moldata = Faunus::molecules.at(molid);
        group.id = molid;
        group.compressible = moldata.compressible;
        group.atomic = moldata.atomic;
        if (group.isMolecular()) {
            if (particles.size() != moldata.atoms.size()) {
                faunus_logger->error("{} requires {} atoms but {} were provided", moldata.name, moldata.atoms.size(),
                                     particles.size());
                throw std::runtime_error("particle size mismatch");
            }
            group.cm = Geometry::massCenter(group.begin(), group.end(), geo.getBoundaryFunc(), -particles.begin()->pos);
        } else {
            if (particles.size() % moldata.atoms.size() != 0) {
                throw std::runtime_error("indivisible by atomic group size: "s + moldata.name);
            }
        }
        groups.push_back(group);
    }
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
        assert(p.size() == other.p.size());
        assert(groups.size() == other.groups.size());
        assert(implicit_reservoir.size() == other.implicit_reservoir.size());
        if (change.dV or change.all) {
            geo = other.geo; // copy simulation geometry
        }
        if (change.all) {                                                   // deep copy *everything*
            implicit_reservoir = other.implicit_reservoir;                  // copy all implicit molecules
            p = other.p;                                                    // copy all positions
            groups = other.groups;                                          // copy all groups
            assert(p.begin() != other.p.begin());                           // check deep copy problem
            assert(groups.front().begin() != other.groups.front().begin()); // check deep copy problem
        } else {
            for (const auto &changed : change.groups) {             // look over changed groups
                auto &group = groups.at(changed.index);             // old group
                auto &other_group = other.groups.at(changed.index); // new group
                assert(group.id == other_group.id);
                if (group.traits().isImplicit()) { // the molecule is implicit
                    implicit_reservoir[group.id] = other.implicit_reservoir.at(group.id);
                } else if (changed.all) {
                    group = other_group;            // copy everything
                } else {                            // copy only a subset
                    group.shallowcopy(other_group); // copy group data but *not* particles
                    for (auto i : changed.atoms) {  // loop over atom index (rel. to group)
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
        group.unwrap(geo.getDistanceFunc()); // ... before scaling the volume
    }
    double Vold = geo.getVolume();             // simulation volume before move
    Point scale = geo.setVolume(Vnew, method); // scale volume of simulation container

    auto scale_position = [&](auto &particle) {
        particle.pos = particle.pos.cwiseProduct(scale);
        geo.boundary(particle.pos);
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
                    Geometry::massCenter(group.begin(), group.end(), geo.getBoundaryFunc(), -original_mass_center);
            } else { // scale mass center; translate positions
                group.cm = group.cm.cwiseProduct(scale);
                geo.boundary(group.cm);
                auto mass_center_displacement = group.cm - original_mass_center;
                std::for_each(group.begin(), group.end(), [&](auto &particle) {
                    particle.pos += mass_center_displacement; // translate internal coordinates
                    geo.boundary(particle.pos);               // apply PBC
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
        Vold = std::pow(Vold, 1. / 3.);  // ?
    }
    for (auto trigger_function : scaleVolumeTriggers) { // external clients may have added function
        trigger_function(*this, Vold, Vnew);            // to be triggered upon each volume change
    }
    return scale;
}

json Space::info() {
    json j = {{"number of particles", p.size()}, {"number of groups", groups.size()}, {"geometry", geo}};
    auto &j_groups = j["groups"];
    for (auto &group : groups) {
        auto &molname = Faunus::molecules.at(group.id).name;
        json tmp, d = group;
        d.erase("cm");
        d.erase("id");
        d.erase("atomic");
        auto ndx = group.to_index(p.begin()); // absolute index
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

Space::Tgvec::iterator Space::randomMolecule(int molid, Random &rand, Space::Selection sel) {
    auto m = findMolecules(molid, sel);
    if (not ranges::cpp20::empty(m))
        return groups.begin() + (&*rand.sample(m.begin(), m.end()) - &*groups.begin());
    return groups.end();
}
const std::map<int, int> &Space::getImplicitReservoir() const { return implicit_reservoir; }

std::map<int, int> &Space::getImplicitReservoir() { return implicit_reservoir; }

std::vector<Space::Tgroup, std::allocator<Space::Tgroup>>::iterator Space::findGroupContaining(const Particle &i) {
    return std::find_if(groups.begin(), groups.end(), [&i](auto &g) { return g.contains(i); });
}

std::vector<Space::Tgroup, std::allocator<Space::Tgroup>>::iterator Space::findGroupContaining(size_t atom_index) {
    assert(atom_index < p.size());
    return std::find_if(groups.begin(), groups.end(),
                        [&](auto &g) { return atom_index < std::distance(p.begin(), g.end()); });
}

size_t Space::numParticles(Space::Selection selection) const {
    if (selection == Selection::ALL) {
        return p.size();
    } else if (selection == Selection::ACTIVE) {
        return std::accumulate(groups.begin(), groups.end(), 0, [](auto sum, auto &g) { return sum + g.size(); });
    } else {
        throw std::runtime_error("invalid selection");
    }
}

/**
 * @throw If group is not part of space
 */
int Space::getGroupIndex(const Space::Tgroup& group) const {
    auto index = std::addressof(group) - std::addressof(groups.front()); // std::ptrdiff_t
    assert(std::abs(index) <= std::numeric_limits<int>::max());
    if (index < 0 or index >= groups.size()) {
        throw std::out_of_range("invalid group index");
    }
    return static_cast<int>(index);
}
int Space::getFirstParticleIndex(const Tgroup& group) const {
    return std::distance<ParticleVector::const_iterator>(p.cbegin(), group.begin());
}

/**
 * Returns the index of the first particle of the group in the range returned by `activeParticles()`
 */
int Space::getFirstActiveParticleIndex(const Tgroup& group) const {
    const auto group_index = getGroupIndex(group);
    int index = 0;
    std::for_each(groups.begin(), groups.begin() + group_index, [&](auto& g) { index += (int)g.size(); });
    return index;
}

TEST_CASE("Space::numParticles") {
    Space spc;
    spc.p.resize(10);
    CHECK(spc.p.size() == spc.numParticles(Space::Selection::ALL));
    CHECK(spc.numParticles(Space::Selection::ACTIVE) == 0);  // zero as there are still no groups
    spc.groups.emplace_back(spc.p.begin(), spc.p.end() - 2); // enclose first 8 particles in group
    CHECK(spc.numParticles(Space::Selection::ACTIVE) == 8);
    spc.groups.emplace_back(spc.p.end() - 2, spc.p.end()); // enclose last 2 particles in group
    CHECK(spc.numParticles(Space::Selection::ACTIVE) == 10);
    spc.groups.front().resize(0); // deactivate first group with 8 particles
    CHECK(spc.numParticles(Space::Selection::ACTIVE) == 2);
}

void to_json(json& j, const Space& spc) {
    j["geometry"] = spc.geo;
    j["groups"] = spc.groups;
    j["particles"] = spc.p;
    j["reactionlist"] = reactions;
    j["implicit_reservoir"] = spc.getImplicitReservoir();
}
void from_json(const json &j, Space &spc) {
    typedef typename Space::Tpvec Tpvec;
    using namespace std::string_literals;

    try {
        if (atoms.empty()) {
            atoms = j.at("atomlist").get<decltype(atoms)>();
        }
        if (molecules.empty()) {
            molecules = j.at("moleculelist").get<decltype(molecules)>();
        }
        if (reactions.empty()) {
            if (j.count("reactionlist") > 0) {
                reactions = j.at("reactionlist").get<decltype(reactions)>();
            }
        }

        spc.clear();
        spc.geo = j.at("geometry");

        if (j.count("groups") == 0) {
            InsertMoleculesInSpace::insertMolecules(j.at("insertmolecules"), spc);
        } else {
            spc.p = j.at("particles").get<Tpvec>();
            if (!spc.p.empty()) {
                auto begin = spc.p.begin();
                Space::Tgroup g(begin, begin); // create new grou[
                for (auto &i : j.at("groups")) {
                    g.begin() = begin;
                    from_json(i, g);
                    spc.groups.push_back(g);
                    begin = g.trueend();
                }
                if (begin != spc.p.end()) {
                    throw ConfigurationError("load error");
                }
            }
        }

        if (auto it = j.find("implicit_reservoir"); it != j.end()) {
            assert(it->is_array());
            for (auto vec : *it) {
                assert(vec.is_array() && vec.size() == 2);
                spc.getImplicitReservoir()[vec[0]] = vec[1];
            }
            faunus_logger->trace("{} implicit molecules loaded", it->size());
        }

        // check correctness of molecular mass centers
        auto active_and_molecular = [](const auto& group) { return (!group.empty() && group.isMolecular()); };
        for (const auto& group : spc.groups | ranges::cpp20::views::filter(active_and_molecular)) {
            const auto should_be_small = spc.geo.sqdist(
                group.cm, Geometry::massCenter(group.begin(), group.end(), spc.geo.getBoundaryFunc(), -group.cm));
            if (should_be_small > 1e-9) {
                throw std::runtime_error(fmt::format(
                    "couldn't calculate mass center for {}; increase periodic box size?", group.traits().name));
            }
        }
    } catch (std::exception& e) {
        std::throw_with_nested(std::runtime_error("error building space"));
    }
}
ActiveParticles::const_iterator::const_iterator(const Space &spc, ActiveParticles::const_iterator::Tparticle_iter it)
    : spc(spc), particle_iter(it) {
    groups_iter = spc.groups.begin();
}

ActiveParticles::const_iterator ActiveParticles::const_iterator::operator++() { // advance particles and groups
    if (++particle_iter == groups_iter->end()) {
        do {
            if (++groups_iter == spc.groups.end()) {
                return *this;
            }
        } while (groups_iter->empty());
        particle_iter = groups_iter->begin();
    }
    return *this;
}
ActiveParticles::const_iterator ActiveParticles::begin() const {
    return const_iterator(spc, spc.p.begin());
}

ActiveParticles::const_iterator ActiveParticles::end() const {
    return spc.groups.empty() ? begin() : const_iterator(spc, spc.groups.back().end());
}

size_t ActiveParticles::size() const {
    return std::accumulate(spc.groups.begin(), spc.groups.end(), 0,
                           [](size_t sum, const auto &g) { return sum + g.size(); });
}

ActiveParticles::ActiveParticles(const Space &spc) : spc(spc) {}

TEST_SUITE_BEGIN("Space");

TEST_CASE("[Faunus] Space") {
    using doctest::Approx;
    Space spc1;
    spc1.geo = R"( {"type": "sphere", "radius": 1e9} )"_json;

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
    spc1.push_back(0, p); // insert molecular group
    CHECK(spc1.p.size() == 2);
    CHECK(spc1.groups.size() == 1);
    CHECK(spc1.groups.back().isMolecular());
    CHECK(spc1.groups.front().id == 0);
    CHECK(spc1.groups.front().cm.x() == doctest::Approx(2.5));

    // check `positions()`
    CHECK(&spc1.positions()[0] == &spc1.p[0].pos);

    SUBCASE("ActiveParticles") {
        // add three groups to space
        Space spc;
        spc.geo = R"( {"type": "sphere", "radius": 1e9} )"_json;
        Particle a;
        a.pos.setZero();
        a.id = 0;
        ParticleVector pvec({a, a, a});

        Faunus::molecules.at(0).atoms.resize(3);

        spc.push_back(0, pvec);
        spc.push_back(0, pvec);
        spc.push_back(0, pvec);

        spc.groups.at(0).id = 0;
        spc.groups.at(1).id = 1;
        spc.groups.at(2).id = 0;

        for (size_t i = 0; i < spc.p.size(); i++)
            spc.p[i].charge = double(i);

        CHECK(spc.p.size() == 9);
        CHECK(spc.groups.size() == 3);

        CHECK(spc.numMolecules<Space::Tgroup::ANY>(0) == 2);
        CHECK(spc.numMolecules<Space::Tgroup::ANY>(1) == 1);

        spc.groups[0].deactivate(spc.p.begin(), spc.p.begin() + 1);
        spc.groups[1].deactivate(spc.groups[1].begin(), spc.groups[1].end());

        CHECK(spc.groups[0].size() == 2);
        CHECK(spc.groups[1].size() == 0);
        CHECK(spc.groups[2].size() == 3);

        CHECK(spc.numMolecules<Space::Tgroup::ACTIVE>(0) == 2);
        CHECK(spc.numMolecules<Space::Tgroup::ACTIVE | Space::Tgroup::NEUTRAL>(0) == 0);
        CHECK(spc.numMolecules<Space::Tgroup::ACTIVE>(1) == 0);
        CHECK(spc.numMolecules<Space::Tgroup::INACTIVE>(1) == 1);
        CHECK(spc.numMolecules<Space::Tgroup::INACTIVE | Space::Tgroup::NEUTRAL>(1) == 0);

        auto p = ActiveParticles(spc);
        size_t size = 0;
        std::vector<int> vals;
        for (const auto &i : p) {
            size++;
            vals.push_back(int(i.charge));
        }

        CHECK(vals == std::vector<int>({1, 2, 6, 7, 8}));
        CHECK(size == p.size());

        // now let's check the rangev3 implementation
        // in `activeParticles()`:
        auto p2 = spc.activeParticles();
        CHECK(std::distance(p2.begin(), p2.end()) == size);
        vals.clear();
        for (const auto &i : p2) {
            vals.push_back(int(i.charge));
        }
        CHECK(vals == std::vector<int>({1, 2, 6, 7, 8}));
    }
}

TEST_CASE("[Faunus] Space::updateParticles") {
    using doctest::Approx;
    Space spc;
    Geometry::Cuboid geo({100, 100, 100});
    spc.geo = Geometry::Chameleon(geo, Geometry::Variant::CUBOID);

    spc.p.resize(2);

    ParticleVector p(2);
    p[0].pos = {0, 0, 0};
    p[1].pos = {2, 0, 0};

    spc.updateParticles(p.begin(), p.end(), spc.p.begin());
    CHECK(spc.p[0].pos.x() == p[0].pos.x());
    CHECK(spc.p[1].pos.x() == p[1].pos.x());

    std::vector<Point> positions = {{2.1, 0, 0}, {0.9, 0, 0}};
    spc.updateParticles(positions.begin(), positions.end(), spc.p.begin(),
                        [](const auto &pos, auto &particle) { particle.pos = pos; });
    CHECK(spc.p[0].pos.x() == 2.1);
    CHECK(spc.p[1].pos.x() == 0.9);

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
void makeNaCl(Space &space, int num_particles, const Geometry::Chameleon &geometry) {
    pc::temperature = 298.15_K;
    space.geo = geometry;

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
void makeWater(Space &space, int num_particles, const Geometry::Chameleon &geometry) {
    pc::temperature = 298.15_K;
    space.geo = geometry;

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
 * @throws If number if inactive is higher than number of active molecules
 *
 * The atoms in the atomic groups are repeated `num_molecules` times, then inserted
 * into Space.
 */
void InsertMoleculesInSpace::insertAtomicGroups(MoleculeData &moldata, Space &spc, int num_molecules,
                                                int num_inactive_molecules) {
    assert(moldata.atomic);
    ParticleVector repeated_particles;
    repeated_particles.reserve(num_molecules * moldata.atoms.size()); // prepare memory
    for (size_t i = 0; i < num_molecules; i++) {                      // repeat insertion into the same atomic group
        auto particles = moldata.getRandomConformation(spc.geo, spc.p);
        repeated_particles.insert(repeated_particles.end(), particles.begin(), particles.end());
    }

    spc.push_back(moldata.id(), repeated_particles); // create new group in Space

    if (num_inactive_molecules > num_molecules) {
        throw std::runtime_error("too many inactive molecules requested");
    } else if (num_inactive_molecules > 0) {
        int num_active_atoms = (num_molecules - num_inactive_molecules) * moldata.atoms.size();
        spc.groups.back().resize(num_active_atoms); // deactivate atoms in molecule
    }
}

void InsertMoleculesInSpace::insertMolecularGroups(MoleculeData &moldata, Space &spc, int num_molecules,
                                                   int num_inactive) {
    assert(moldata.atomic == false);
    for (size_t i = 0; i < num_molecules; i++) { // insert molecules
        spc.push_back(moldata.id(), moldata.getRandomConformation(spc.geo, spc.p));
    }
    if (num_inactive > num_molecules) {
        throw std::runtime_error("too many inactive molecules requested");
    } else { // deactivate groups, starting from the back
        std::for_each(spc.groups.rbegin(), spc.groups.rbegin() + num_inactive, [&](auto &group) {
            group.unwrap(spc.geo.getDistanceFunc()); // make molecules whole (remove PBC) ...
            group.resize(0);                         // ... and then deactivate
        });
    }
}

/**
 * @brief Set positions for num_molecules last groups in Space
 * @param spc Space to insert into
 * @param int num_molecules Number of groups set affect
 * @param particles Position vector for all particles in the N groups
 * @param offset Translate positions by this offset
 * @throws if num_molecules doesn't match, or if positions are outside simulation cell
 *
 * Sets particle positions in space using a given input
 * particle vector. The vector can span several *identical* groups.
 * The number of groups affected must be given in order to update their mass-centers.
 * Only *positions* are affected.
 */
void InsertMoleculesInSpace::setPositionsForTrailingGroups(Space &spc, int num_molecules,
                                                           const Faunus::ParticleVector &particles,
                                                           const Point &offset) {
    assert(spc.groups.size() >= num_molecules);
    if (particles.size() != num_molecules * (spc.groups.rbegin())->traits().atoms.size()) {
        throw std::runtime_error("number of particles doesn't match groups");
    } else {
        // update positions in space, starting from the back
        std::transform(particles.rbegin(), particles.rend(), spc.p.rbegin(), spc.p.rbegin(), [&](auto &src, auto &dst) {
            dst.pos = src.pos + offset;       // shift by offset
            if (spc.geo.collision(dst.pos)) { // check if position is inside simulation volume
                throw std::runtime_error("positions outside box");
            }
            return dst;
        });
        // update mass-centers on modified groups; start from the back
        std::for_each(spc.groups.rbegin(), spc.groups.rbegin() + num_molecules, [&](auto &g) {
            g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.begin()->pos);
        });
    }
}

/**
 * @brief Insert implicit groups into Space
 * @param moldata Molecule to insert. Must be implicit.
 * @param spc Space to insert into.
 * @param num_molecules Nunber of implicit molecules to insert
 */
void InsertMoleculesInSpace::insertImplicitGroups(const MoleculeData &moldata, Space &spc, int num_molecules) {
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
    const auto num_molecules = getNumberOfMolecules(properties, spc.geo.getVolume(), molname);
    if (num_molecules == 0) {
        if (!moldata.isImplicit()) {
            throw ConfigurationError("one or more {} molecule(s) required; concentration too low?", molname);
        }
    } else {
        const auto num_inactive = getNumberOfInactiveMolecules(properties, num_molecules);
        const auto molarity = (num_molecules - num_inactive) / spc.geo.getVolume() / 1.0_molar;
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
 * @throws If the 'positions' key is there, but could not be loaded
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
    spc.clear();
    if (!j.is_array()) {
        throw ConfigurationError("molecules to insert must be an array");
    } else {
        for (const auto &item : j) { // loop over array of molecules
            if (item.is_object() && item.size() == 1) {
                for (auto &[molecule_name, properties] : item.items()) {
                    try {
                        insertItem(molecule_name, properties, spc);
                    } catch (std::exception &e) {
                        throw ConfigurationError("error inserting {}: {}", molecule_name, e.what());
                    }
                }
            } else {
                throw ConfigurationError("syntax error inserting molecules");
            }
        }
    }
}

/**
 * @param j Input json object
 * @param volume Volume of simulation container needed to calculate concentration
 * @param molecule_name Name of molecule needed for logging
 * @return Number of molecules to insert
 *
 * Looks for json key `N` or `molarity`. For the latter, the nearest corresponding
 * number of particles is calculated based on the given system volume.
 */
int InsertMoleculesInSpace::getNumberOfMolecules(const json &j, double volume, const std::string &molecule_name) {
    const double error_limit = 0.01; // warn if relative density error is above this
    int num_molecules = 0;
    if (j.contains("N")) {
        num_molecules = j.at("N").get<int>();
    } else {
        auto density = j.at("molarity").get<double>() * 1.0_molar;
        num_molecules = std::round(density * volume);
        if (double error = (density - num_molecules / volume) / density; error > error_limit) {
            faunus_logger->warn("{}: initial molarity differs by {}% from target value", molecule_name, error * 100);
        }
    }
    return num_molecules;
}

/**
 * @param j Input json object
 * @param number_of_molecules Total number of molecules
 * @return Number of molecules to be inactive (always smaller than `number_of_molecules`)
 * @throws If Inactive molecules is higher than `number_of_molecules`
 *
 * Looks for key "inactive" and if:
 * - boolean true: all molecules are inactive
 * - number: number of molecules to declare inactive
 * - no `inactive` key found, return zero
 */
int InsertMoleculesInSpace::getNumberOfInactiveMolecules(const json &j, int number_of_molecules) {
    int number_of_inactive_molecules = 0; // number of inactive molecules
    if (auto it = j.find("inactive"); it != j.end()) {
        if (it->is_boolean()) {
            if (*it) {
                number_of_inactive_molecules = number_of_molecules; // all molecules are inactive
            }
        } else if (it->is_number_integer()) {
            number_of_inactive_molecules = *it; // a subset are inactive
            if (number_of_inactive_molecules > number_of_molecules) {
                throw ConfigurationError("too many inactive particles requested");
            }
        }
    }
    return number_of_inactive_molecules;
}
} // namespace Faunus
