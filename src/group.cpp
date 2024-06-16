#include "group.h"
#include "aux/eigensupport.h"
#include <range/v3/view/sample.hpp>
#include <range/v3/view/common.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <nlohmann/json.hpp>
#include <doctest/doctest.h>

namespace Faunus {

Group::Group(Group& other)
    : base(other.begin(), other.trueend())
{
    *this = operator=(other);
}

Group::Group(const Group& other)
    : base(other.begin(), other.trueend())
{
    *this = operator=(other);
}

/**
 * @param molid Molecule id the group points to, i.e. a valid index in global `Faunus::molecules`.
 * @param begin Iterator to first particle
 * @param end Iterator to (beyond) end particle
 * @throw if molid is out of range w. respect to `Faunus::molecules`
 */
Group::Group(MoleculeData::index_type molid, Group::iter begin, Group::iter end)
    : base(begin, end)
    , id(molid)
{
    if (id >= Faunus::molecules.size()) {
        throw std::range_error("invalid molecule id");
    }
}

/**
 * Performs a deep copy from another group
 *
 * @throw if the capacities of the two groups differ
 */
Group& Group::operator=(const Group& other)
{
    if (&other != this) {
        shallowCopy(other);
        if (other.begin() != begin()) {
            std::copy(other.begin(), other.trueend(), begin());
        } // copy all particle data
    }
    return *this;
}

/**
 * This copies size, id, mass center etc. from another group
 * but does *not* copy particle data.
 *
 * @throw if the capacities of the two groups differ
 */
Group& Group::shallowCopy(const Group& other)
{
    if (&other != this) {
        if (capacity() != other.capacity()) {
            throw std::runtime_error("Group::shallowCopy: capacity mismatch");
        }
        resize(other.size());
        id = other.id;
        mass_center = other.mass_center;
        conformation_id = other.conformation_id;
    }
    return *this;
}

bool Group::contains(const Particle& particle, bool include_inactive) const
{
    const auto size = (include_inactive ? capacity() : this->size());
    if (size > 0) {
        auto index = std::addressof(particle) - std::addressof(*begin());
        return (index >= 0 and index < size);
    }
    return false;
}

AtomData::index_type Group::getParticleIndex(const Particle& particle, bool include_inactive) const
{
    if (!empty()) {
        const auto index = std::addressof(particle) - std::addressof(*begin()); // std::ptrdiff_t
        const auto group_size = (include_inactive ? capacity() : size());
        if (index >= 0 && index < group_size) {
            return static_cast<AtomData::index_type>(index);
        }
    }
    throw std::out_of_range("invalid particle index or group is empty");
}

double Group::mass() const
{
    return std::accumulate(begin(), end(), 0.0, [](double sum, auto& particle) { return sum + particle.traits().mw; });
}

[[maybe_unused]] void Group::wrap(Geometry::BoundaryFunction boundary)
{
    boundary(mass_center);
    for (auto& particle : *this) {
        boundary(particle.pos);
    }
}

void Group::rotate(const Eigen::Quaterniond& quaternion, Geometry::BoundaryFunction boundary)
{
    Geometry::rotate(begin(), end(), quaternion, boundary, -mass_center);
}

void Group::translate(const Point& displacement, Geometry::BoundaryFunction boundary)
{
    mass_center += displacement;
    boundary(mass_center);
    for (auto& particle : *this) {
        particle.pos += displacement;
        boundary(particle.pos);
    }
}

/**
 * @param boundary_function Function to apply periodic boundaries
 * @param approximate_mass_center Original or appriximate mass center used for PBC removal
 *
 * Only active, molecular groups are affected. Before the mass center is
 * calculated, the molecule is translated towards the center of the simulation
 * box to remove possible periodic boundary conditions; then translated back again.
 * The translation is done by subtracting / adding `approximate_mass_center`.
 * Safely handles empty groups.
 */
void Group::updateMassCenter(Geometry::BoundaryFunction boundary_function, const Point& approximate_mass_center)
{
    if (isMolecular() && !empty()) {
        mass_center = Geometry::massCenter(begin(), end(), boundary_function, -approximate_mass_center);
    }
}

/**
 * @param boundary_function Function to apply periodic boundaries
 *
 * This will approximate the existing mass center by the middle particle which for PBC systems
 * is generally safer then using the old mass center.
 */
void Group::updateMassCenter(Geometry::BoundaryFunction boundary_function)
{
    if (empty()) {
        return;
    }
    const auto& approximate_mass_center = at(size() / 2).pos;
    updateMassCenter(boundary_function, approximate_mass_center);
}

/**
 * @brief Returns i'th element in group
 * @param index index starting at zero
 * @return reference to value at i'th element
 * @throw if out of interval `[0:capacity[`
 */
Particle& Group::at(size_t index)
{
    if (index >= capacity()) {
        throw std::out_of_range("group index out of range");
    }
    return operator[](index);
}

/**
 * @brief Returns i'th element in group
 * @param index index starting at zero
 * @return reference to value at i'th element
 * @throw if out of interval `[0:capacity[`
 */
const Particle& Group::at(size_t index) const
{
    if (index >= capacity()) {
        throw std::out_of_range("group index out of range");
    }
    return operator[](index);
}

/**
 * If the group supports mass centers (i.e. molecular groups), an optional reference
 * will be returned. Example usage:
 *
 *    if (auto mass_center = group.massCenter()) {
 *        (*mass_center).get() = Point(0,1,2);
 *    }
 *
 * @return Optional reference to stored group mass center
 */
std::optional<std::reference_wrapper<Point>> Group::massCenter()
{
    if (isMolecular()) {
        return std::ref(mass_center);
    }
    return std::nullopt;
}

std::optional<std::reference_wrapper<const Point>> Group::massCenter() const
{
    if (isMolecular()) {
        return std::cref(mass_center);
    }
    return std::nullopt;
}

bool Group::isFull() const { return size() == capacity(); }

/**
 * @param mask Bitmask based on enum `Group::Selectors`
 * @return Lambda function that returns true if group matches mask
 *
 * ~~~ cpp
 * auto filter = Group::getSelectionFilter(Select::ACTIVE | Select::NEUTRAL);
 * bool b = filter(mygroup); // true if mygroup is active and uncharged
 * ~~~
 */
// constexpr std::function<bool(const Group &)> getGroupFilter(unsigned int mask) {
//    return [mask = mask](const Group &g) { return g.match<mask>(); };
//}

void to_json(json& j, const Group& group)
{
    j = {{"id", group.id},
         {"cm", group.mass_center},
         {"atomic", group.isAtomic()},
         {"compressible", group.traits().compressible},
         {"size", group.size()}};
    if (group.capacity() > group.size()) {
        j["capacity"] = group.capacity();
    }
    if (group.conformation_id != 0) {
        j["confid"] = group.conformation_id;
    }
}

void from_json(const json& j, Group& group)
{
    group.resize(j.at("size").get<size_t>());
    group.trueend() = group.begin() + j.value("capacity", group.size());
    group.id = j.at("id").get<decltype(group.id)>();
    group.mass_center = j.at("cm").get<Point>();
    group.conformation_id = j.value("confid", 0);
}

using doctest::Approx;

TEST_SUITE_BEGIN("Group");

TEST_CASE("[Faunus] swap_to_back")
{
    using VecInt = std::vector<int>;
    VecInt v = {1, 2, 3, 4};

    swap_to_back(v.begin(), v.end(), v.end());
    CHECK_EQ(v, VecInt({1, 2, 3, 4}));

    std::sort(v.begin(), v.end());
    swap_to_back(v.begin() + 1, v.begin() + 3, v.end());
    CHECK_EQ(v, VecInt({1, 4, 3, 2}));
}

TEST_CASE("[Faunus] ElasticRange")
{
    std::vector<int> v = {10, 20, 30, 40, 50, 60};
    ElasticRange<int> r(v.begin(), v.end());
    CHECK_EQ(r.size(), 6);
    CHECK_EQ(r.empty(), false);
    CHECK_EQ(r.size(), r.capacity());
    *r.begin() += 1;
    CHECK_EQ(v[0], 11);

    r.deactivate(r.begin(), r.end());
    CHECK_EQ(r.size(), 0);
    CHECK_EQ(r.empty(), true);
    CHECK_EQ(r.capacity(), 6);
    CHECK_EQ(r.inactive().size(), 6);
    CHECK_EQ(r.begin(), r.end());

    r.activate(r.inactive().begin(), r.inactive().end());
    CHECK_EQ(r.size(), 6);
    CHECK_EQ(std::is_sorted(r.begin(), r.end()), true); // back to original

    r.deactivate(r.begin() + 1, r.begin() + 3);
    CHECK_EQ(r.size(), 4);
    CHECK_EQ(std::find(r.begin(), r.end(), 20), r.end());
    CHECK_EQ(std::find(r.begin(), r.end(), 30), r.end());
    CHECK_EQ(*r.end(), 20); // deactivated elements can be retrieved from `end()`
    CHECK_EQ(*(r.end() + 1), 30);

    auto ipair = r.to_index(v.begin());
    CHECK_EQ(ipair.first, 0);
    CHECK_EQ(ipair.second, 3);

    r.activate(r.end(), r.end() + 2);
    CHECK_EQ(*(r.end() - 2), 20); // activated elements can be retrieved from `end()-n`
    CHECK_EQ(*(r.end() - 1), 30);
    CHECK_EQ(r.size(), 6);

    // check relocation
    auto v2 = v;
    v2.front() = -7;
    CHECK((*r.begin() != -7));
    r.relocate(v.begin(), v2.begin());
    CHECK_EQ(*r.begin(), -7);
}

TEST_CASE("[Faunus] Group")
{
    Random rand;
    std::vector<Particle> p(3);
    p.reserve(10);
    p[0].id = 0;
    p[1].id = 1;
    p[2].id = 1;
    if (Faunus::molecules.empty()) {
        Faunus::molecules.resize(1);
    }
    Group g(0, p.begin(), p.end());

    SUBCASE("contains()")
    {
        CHECK(g.contains(p[0]));
        CHECK(g.contains(p[1]));
        CHECK(g.contains(p[2]));
        CHECK_EQ(g.size(), 3);
        g.deactivate(g.end() - 1, g.end());
        CHECK_EQ(g.size(), 2);
        CHECK_EQ(g.contains(p[2]), false);
        CHECK_EQ(g.contains(p[2], true), true);
        g.activate(g.end(), g.end() + 1);
        CHECK_EQ(g.size(), 3);
    }

    SUBCASE("getParticleIndex()")
    {
        Group gg(0, p.begin(), p.end());
        CHECK_EQ(gg.getParticleIndex(p[0]), 0);
        CHECK_EQ(gg.getParticleIndex(p[1]), 1);
        CHECK_EQ(gg.getParticleIndex(p[2]), 2);
        CHECK_EQ(gg.size(), 3);
        gg.deactivate(gg.end() - 1, gg.end());
        CHECK_EQ(gg.size(), 2);
        CHECK_THROWS(std::ignore = gg.getParticleIndex(p[2]));
        gg.resize(0);
        CHECK_THROWS(std::ignore = gg.getParticleIndex(p[0]));
    }

    SUBCASE("getGroupFilter(): complete group")
    {
        using T = Group;
        auto filter = getGroupFilter<T::Selectors::ACTIVE>();
        CHECK_EQ(filter(g), true);
        filter = getGroupFilter<T::Selectors::FULL>();
        CHECK_EQ(filter(g), true);
        filter = getGroupFilter<T::Selectors::INACTIVE>();
        CHECK_EQ(filter(g), false);
        filter = getGroupFilter<T::Selectors::ACTIVE | T::Selectors::NEUTRAL>();
        CHECK_EQ(filter(g), true);
        filter = getGroupFilter<T::Selectors::ACTIVE | T::Selectors::MOLECULAR>();
        CHECK_EQ(filter(g), g.isMolecular());
        filter = getGroupFilter<T::Selectors::INACTIVE | T::Selectors::MOLECULAR>();
        CHECK_EQ(filter(g), g.isMolecular());
        filter = getGroupFilter<T::Selectors::ACTIVE | T::Selectors::ATOMIC>();
        CHECK_EQ(filter(g), g.isAtomic());

        g.begin()->charge = 0.1;
        filter = getGroupFilter<T::Selectors::ACTIVE | T::Selectors::NEUTRAL>();
        CHECK_EQ(filter(g), false);
        g.begin()->charge = 0.0;
    }

    // find all elements with id=1
    auto slice1 = g.findAtomID(1);
    CHECK_EQ(std::distance(slice1.begin(), slice1.end()), 2);

    // find *one* random value with id=1
    auto slice2 = slice1 | ranges::views::sample(1, rand.engine) | ranges::views::common;
    CHECK_EQ(std::distance(slice2.begin(), slice2.end()), 1);

    // check rotation
    Eigen::Quaterniond q;
    q = Eigen::AngleAxisd(pc::pi / 2, Point(1, 0, 0));
    p.at(0).pos = p.at(0).getExt().mu = p.at(0).getExt().scdir = {0, 1, 0};

    Geometry::Chameleon geo = R"({"type":"cuboid", "length": [2,2,2]})"_json;
    g.rotate(q, geo.getBoundaryFunc());
    CHECK_EQ(p[0].pos.y(), doctest::Approx(0));
    CHECK_EQ(p[0].pos.z(), doctest::Approx(1));
    CHECK_EQ(p[0].getExt().mu.y(), doctest::Approx(0));
    CHECK_EQ(p[0].getExt().mu.z(), doctest::Approx(1));
    CHECK_EQ(p[0].getExt().scdir.y(), doctest::Approx(0));
    CHECK_EQ(p[0].getExt().scdir.z(), doctest::Approx(1));

    p[0].pos = {1, 2, 3};
    p[1].pos = {4, 5, 6};

    // iterate over positions and modify them
    for (auto& pos : g.positions()) {
        pos *= 2.0;
    }
    CHECK_EQ(p[1].pos.x(), doctest::Approx(8));
    CHECK_EQ(p[1].pos.y(), doctest::Approx(10));
    CHECK_EQ(p[1].pos.z(), doctest::Approx(12));

    SUBCASE("operator[]")
    {
        CHECK_EQ(p.begin(), g.begin());
        CHECK_EQ(p.end(), g.end());

        // a new range by using an index filter
        std::vector<size_t> index = {0, 1};
        auto subset = g[index];
        CHECK_EQ(subset.size(), 2);
        CHECK_EQ(&(*p.begin()), &(*subset.begin()));
        CHECK_EQ(&(*(p.begin() + 1)), &(*(subset.begin() + 1)));
        for (auto& i : subset)
            i.pos *= 2;
        CHECK_EQ(p[1].pos.x(), doctest::Approx(16));
        CHECK_EQ(p[1].pos.y(), doctest::Approx(20));
        CHECK_EQ(p[1].pos.z(), doctest::Approx(24));
    }

    SUBCASE("deep copy and resizing")
    {
        std::vector<Particle> p1(5), p2(5);
        p1.front().id = 1;
        p2.front().id = -1;

        Group g1(0, p1.begin(), p1.end());
        Group g2(0, p2.begin(), p2.end());

        g2.id = 100;
        g2.mass_center = {1, 0, 0};
        g2.conformation_id = 20;
        g1 = g2;

        CHECK_EQ(g1.id, 100);
        CHECK_EQ(g1.mass_center.x(), 1);
        CHECK_EQ(g1.conformation_id, 20);

        CHECK_EQ((*g1.begin()).id, -1);
        CHECK_EQ((*g2.begin()).id, -1);
        CHECK((g1.begin() != g2.begin()));
        CHECK_EQ(g1.size(), g2.size());
        (*g2.begin()).id = 10;
        g2.resize(4);
        g1 = g2;
        CHECK_EQ(g1.size(), 4);
        CHECK_EQ(g1.capacity(), 5);
        CHECK_EQ(p1.front().id, 10);

        g1.id = 0;

        SUBCASE("getGroupFilter(): incomplete group")
        {
            CHECK(!Faunus::molecules.empty());
            using Tgroup = Group;
            auto filter = getGroupFilter<Tgroup::FULL>();
            CHECK_EQ(filter(g1), false);
            filter = getGroupFilter<Tgroup::INACTIVE>();
            CHECK_EQ(filter(g1), false);
            filter = getGroupFilter<Tgroup::ACTIVE>();
            CHECK_EQ(filter(g1), true);
            filter = getGroupFilter<Tgroup::ACTIVE | Tgroup::ATOMIC>();
            CHECK_EQ(filter(g1), g1.isAtomic());
            filter = getGroupFilter<Tgroup::ACTIVE | Tgroup::MOLECULAR>();
            CHECK_EQ(filter(g1), g1.isMolecular());
        }

        std::vector<Group> gvec1, gvec2;
        gvec1.push_back(g1);
        gvec2.push_back(g2);
        p2.front().id = 21;

        CHECK_EQ((*(gvec1.front().begin())).id, 10);
        CHECK_EQ((*(gvec2.front().begin())).id, 21);

        // existing groups point to existing particles when overwritten
        gvec1 = gvec2; // invoke *deep* copy of all contained groups
        CHECK((gvec1[0].begin() != gvec2[0].begin()));
        CHECK_EQ(p1.front().id, 21);

        // new groups point to same particles as original
        auto gvec3 = gvec1;
        CHECK_EQ((*gvec1[0].begin()).id, (*gvec3[0].begin()).id);
    }

    SUBCASE("cerial serialisation")
    {
        std::ostringstream out(std::stringstream::binary);
        { // serialize g2
            std::vector<Particle> p2(5);
            Group g2(0, p2.begin(), p2.end());
            p2.front().id = 8;
            p2.back().pos.x() = -10;
            g2.id = 100;
            g2.mass_center = {1, 0, 0};
            g2.conformation_id = 20;
            g2.resize(4);
            cereal::BinaryOutputArchive archive(out);
            archive(g2);
        }

        {                                     // deserialize into g1
            std::istringstream in(out.str()); // not pretty...
            cereal::BinaryInputArchive archive(in);
            std::vector<Particle> p1(5);
            Group g1(0, p1.begin(), p1.end());
            archive(g1);

            CHECK_EQ(g1.id, 100);
            CHECK_EQ(g1.mass_center.x(), 1);
            CHECK_EQ(g1.conformation_id, 20);
            CHECK_EQ(g1.size(), 4);
            CHECK_EQ(g1.capacity(), 5);
            CHECK_EQ(g1.begin()->id, 8);
            CHECK_EQ(p1.front().id, 8);
            CHECK_EQ(p1.back().pos.x(), -10);
            CHECK_EQ(p1.back().ext, nullptr);
        }
    }
}

TEST_SUITE_END();

} // namespace Faunus