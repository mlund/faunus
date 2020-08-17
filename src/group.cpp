#include "group.h"
#include "aux/eigensupport.h"
#include <Eigen/Geometry>
#include <range/v3/view/sample.hpp>
#include <range/v3/view/bounded.hpp>
#include <cereal/archives/binary.hpp>
#include <nlohmann/json.hpp>

namespace Faunus {

template <class T> Group<T>::Group(Group<T> &o) : base(o.begin(), o.trueend()) { *this = operator=(o); }

template <class T> Group<T>::Group(const Group<T> &o) : base(o.begin(), o.trueend()) { *this = operator=(o); }

template <class T> Group<T>::Group(Group<T>::iter begin, Group<T>::iter end) : base(begin,end) {}

template <class T> Group<T> &Group<T>::operator=(const Group<T> &o) {
    if (&o == this)
        return *this;
    shallowcopy(o);
    if (o.begin()!=begin())
        std::copy(o.begin(), o.trueend(), begin()); // copy all particle data
    return *this;
}

template <class T> Group<T> &Group<T>::shallowcopy(const Group<T> &o) {
    if (&o != this) {
        if (this->capacity() != o.capacity())
            throw std::runtime_error("Group::shallowcopy: capacity mismatch");
        this->resize(o.size());
        id = o.id;
        atomic = o.atomic;
        compressible = o.compressible;
        cm = o.cm;
        confid = o.confid;
    }
    return *this;
}

template <class T> bool Group<T>::contains(const T &a, bool include_inactive) const {
    int size = (include_inactive ? this->capacity() : this->size());
    if (size > 0) {
        int d = &a - &(*(this->begin()));
        if (d >= 0 and d < size) {
            return true;
        }
    }
    return false;
}

template <class T> double Group<T>::mass() const {
    return std::accumulate(begin(), end(), 0.0, [](double sum, auto &i) { return sum + i.traits().mw; });
}

template <class T> std::vector<std::reference_wrapper<Point>> Group<T>::positions() const {
    return Faunus::member_view(begin(), end(), &Particle::pos);
    //return ranges::view::transform(*this, [](auto &i) -> Point& {return i.pos;});
}

template <class T> void Group<T>::wrap(Geometry::BoundaryFunction boundary) {
    boundary(cm);
    for (auto &i : *this)
        boundary(i.pos);
}

template <class T> void Group<T>::rotate(const Eigen::Quaterniond &Q, Geometry::BoundaryFunction boundary) {
    Geometry::rotate(begin(), end(), Q, boundary, -cm);
}

template <class T> void Group<T>::translate(const Point &d, Geometry::BoundaryFunction boundary) {
    cm += d;
    boundary(cm);
    for (auto &i : *this) {
        i.pos += d;
        boundary(i.pos);
    }
}

/**
 * @param boundary_function Function to apply periodic boundaries
 *
 * Only active, molecular groups are affected. Before the mass center can
 * be calculated, the molecule is translated towards the center of the simulation
 * box to remove possible periodic boundary conditions; then translated back again.
 * The translation is done by subtracting / adding the initial mass center position.
 * Control the shift by setting `cm` just before calling this function, or disable
 * it by setting `cm={0,0,0}`.
 */
template <class T> void Group<T>::updateMassCenter(Geometry::BoundaryFunction boundary_function) {
    if (isMolecular() && !empty()) {
        cm = Geometry::massCenter(begin(), end(), boundary_function, -cm);
    }
}

/**
 * @param mask Bitmask based on enum `Group::Selectors`
 * @return Lambda function that returns true if group matches mask
 *
 * ~~~ cpp
 * auto filter = Group::getSelectionFilter(Select::ACTIVE | Select::NEUTRAL);
 * bool b = filter(mygroup); // true if mygroup is active and uncharged
 * ~~~
 */
// constexpr std::function<bool(const Group<Particle> &)> getGroupFilter(unsigned int mask) {
//    return [mask = mask](const Group<Particle> &g) { return g.match<mask>(); };
//}

template class Group<Particle>;

void to_json(json &j, const Group<Particle> &g) {
    j = {
        {"id", g.id}, {"cm", g.cm}, {"atomic", g.atomic}, {"compressible", g.compressible},  {"size", g.size()}
    };
    if (g.capacity()>g.size())
        j["capacity"] = g.capacity();
    if (g.confid!=0)
        j["confid"] = g.confid;
}

void from_json(const json &j, Group<Particle> &g) {
    g.resize( j.at("size").get<int>() );
    g.trueend() = g.begin() + j.value("capacity", g.size());
    g.id = j.at("id").get<unsigned int>();
    g.cm = j.at("cm").get<Point>();
    g.atomic = j.at("atomic").template get<bool>();
    g.compressible = j.value("compressible", false);
    g.confid = j.value("confid", 0);
}

using doctest::Approx;

TEST_SUITE_BEGIN("Group");

TEST_CASE("[Faunus] swap_to_back") {
    typedef std::vector<int> _T;
    _T v = {1, 2, 3, 4};

    swap_to_back(v.begin(), v.end(), v.end());
    CHECK(v == _T({1, 2, 3, 4}));

    std::sort(v.begin(), v.end());
    swap_to_back(v.begin() + 1, v.begin() + 3, v.end());
    CHECK(v == _T({1, 4, 3, 2}));
}

TEST_CASE("[Faunus] ElasticRange") {
    std::vector<int> v = {10, 20, 30, 40, 50, 60};
    ElasticRange<int> r(v.begin(), v.end());
    CHECK(r.size() == 6);
    CHECK(r.empty() == false);
    CHECK(r.size() == r.capacity());
    *r.begin() += 1;
    CHECK(v[0] == 11);

    r.deactivate(r.begin(), r.end());
    CHECK(r.size() == 0);
    CHECK(r.empty() == true);
    CHECK(r.capacity() == 6);
    CHECK(r.inactive().size() == 6);
    CHECK(r.begin() == r.end());

    r.activate(r.inactive().begin(), r.inactive().end());
    CHECK(r.size() == 6);
    CHECK(std::is_sorted(r.begin(), r.end()) == true); // back to original

    r.deactivate(r.begin() + 1, r.begin() + 3);
    CHECK(r.size() == 4);
    CHECK(std::find(r.begin(), r.end(), 20) == r.end());
    CHECK(std::find(r.begin(), r.end(), 30) == r.end());
    CHECK(*r.end() == 20); // deactivated elements can be retrieved from `end()`
    CHECK(*(r.end() + 1) == 30);

    auto ipair = r.to_index(v.begin());
    CHECK(ipair.first == 0);
    CHECK(ipair.second == 3);

    r.activate(r.end(), r.end() + 2);
    CHECK(*(r.end() - 2) == 20); // activated elements can be retrieved from `end()-n`
    CHECK(*(r.end() - 1) == 30);
    CHECK(r.size() == 6);

    // check relocation
    auto v2 = v;
    v2.front() = -7;
    CHECK(*r.begin() != -7);
    r.relocate(v.begin(), v2.begin());
    CHECK(*r.begin() == -7);
}

TEST_CASE("[Faunus] Group") {
    Random rand;
    std::vector<Particle> p(3);
    p.reserve(10);
    p[0].id = 0;
    p[1].id = 1;
    p[2].id = 1;
    Group<Particle> g(p.begin(), p.end());

    SUBCASE("contains()") {
        CHECK(g.contains(p[0]));
        CHECK(g.contains(p[1]));
        CHECK(g.contains(p[2]));
        CHECK(g.size() == 3);
        g.deactivate(g.end() - 1, g.end());
        CHECK(g.size() == 2);
        CHECK(g.contains(p[2]) == false);
        CHECK(g.contains(p[2], true) == true);
        g.activate(g.end(), g.end() + 1);
        CHECK(g.size() == 3);
    }

    SUBCASE("getGroupFilter(): complete group") {
        typedef Group<Particle> T;
        auto filter = getGroupFilter<T::Selectors::ACTIVE>();
        CHECK(filter(g) == true);
        filter = getGroupFilter<T::Selectors::FULL>();
        CHECK(filter(g) == true);
        filter = getGroupFilter<T::Selectors::INACTIVE>();
        CHECK(filter(g) == false);
        filter = getGroupFilter<T::Selectors::ACTIVE | T::Selectors::NEUTRAL>();
        CHECK(filter(g) == true);
        filter = getGroupFilter<T::Selectors::ACTIVE | T::Selectors::MOLECULAR>();
        CHECK(filter(g) == true);
        filter = getGroupFilter<T::Selectors::INACTIVE | T::Selectors::MOLECULAR>();
        CHECK(filter(g) == false);
        filter = getGroupFilter<T::Selectors::ACTIVE | T::Selectors::ATOMIC>();
        CHECK(filter(g) == false);

        g.begin()->charge = 0.1;
        filter = getGroupFilter<T::Selectors::ACTIVE | T::Selectors::NEUTRAL>();
        CHECK(filter(g) == false);
        g.begin()->charge = 0.0;
    }

    // find all elements with id=1
    auto slice1 = g.find_id(1);
    CHECK(std::distance(slice1.begin(), slice1.end()) == 2);

    // find *one* random value with id=1
    auto slice2 = slice1 | ranges::views::sample(1, rand.engine) | ranges::views::common;
    CHECK(std::distance(slice2.begin(), slice2.end()) == 1);

    // check rotation
    Eigen::Quaterniond q;
    q = Eigen::AngleAxisd(pc::pi / 2, Point(1, 0, 0));
    p.at(0).pos = p.at(0).getExt().mu = p.at(0).getExt().scdir = {0, 1, 0};

    Geometry::Chameleon geo = R"({"type":"cuboid", "length": [2,2,2]})"_json;
    g.rotate(q, geo.getBoundaryFunc());
    CHECK(p[0].pos.y() == doctest::Approx(0));
    CHECK(p[0].pos.z() == doctest::Approx(1));
    CHECK(p[0].getExt().mu.y() == doctest::Approx(0));
    CHECK(p[0].getExt().mu.z() == doctest::Approx(1));
    CHECK(p[0].getExt().scdir.y() == doctest::Approx(0));
    CHECK(p[0].getExt().scdir.z() == doctest::Approx(1));

    p[0].pos = {1, 2, 3};
    p[1].pos = {4, 5, 6};

    // iterate over positions and modify them
    for (Point &i : g.positions())
        i = 2 * i;
    CHECK(p[1].pos.x() == doctest::Approx(8));
    CHECK(p[1].pos.y() == doctest::Approx(10));
    CHECK(p[1].pos.z() == doctest::Approx(12));

    SUBCASE("operator[]") {
        CHECK(p.begin() == g.begin());
        CHECK(p.end() == g.end());

        // a new range by using an index filter
        std::vector<size_t> index = {0, 1};
        auto subset = g[index];
        CHECK(subset.size() == 2);
        CHECK(&(*p.begin()) == &(*subset.begin()));
        CHECK(&(*(p.begin() + 1)) == &(*(subset.begin() + 1)));
        for (auto &i : subset)
            i.pos *= 2;
        CHECK(p[1].pos.x() == doctest::Approx(16));
        CHECK(p[1].pos.y() == doctest::Approx(20));
        CHECK(p[1].pos.z() == doctest::Approx(24));
    }

    SUBCASE("deep copy and resizing") {
        std::vector<Particle> p1(5), p2(5);
        p1.front().id = 1;
        p2.front().id = -1;

        Group<Particle> g1(p1.begin(), p1.end());
        Group<Particle> g2(p2.begin(), p2.end());

        g2.id = 100;
        g2.atomic = true;
        g2.cm = {1, 0, 0};
        g2.confid = 20;
        g1 = g2;

        CHECK(g1.id == 100);
        CHECK(g1.atomic == true);
        CHECK(g1.cm.x() == 1);
        CHECK(g1.confid == 20);

        CHECK((*g1.begin()).id == -1);
        CHECK((*g2.begin()).id == -1);
        CHECK(g1.begin() != g2.begin());
        CHECK(g1.size() == g2.size());
        (*g2.begin()).id = 10;
        g2.resize(4);
        g1 = g2;
        CHECK(g1.size() == 4);
        CHECK(g1.capacity() == 5);
        CHECK(p1.front().id == 10);

        SUBCASE("getGroupFilter(): incomplete group") {
            typedef Group<Particle> Tgroup;
            auto filter = getGroupFilter<Tgroup::FULL>();
            CHECK(filter(g1) == false);
            filter = getGroupFilter<Tgroup::INACTIVE>();
            CHECK(filter(g1) == false);
            filter = getGroupFilter<Tgroup::ACTIVE>();
            CHECK(filter(g1) == true);
            filter = getGroupFilter<Tgroup::ACTIVE | Tgroup::ATOMIC>();
            CHECK(filter(g1) == true);
            filter = getGroupFilter<Tgroup::ACTIVE | Tgroup::MOLECULAR>();
            CHECK(filter(g1) == false);
        }

        std::vector<Group<Particle>> gvec1, gvec2;
        gvec1.push_back(g1);
        gvec2.push_back(g2);
        p2.front().id = 21;

        CHECK((*(gvec1.front().begin())).id == 10);
        CHECK((*(gvec2.front().begin())).id == 21);

        // existing groups point to existing particles when overwritten
        gvec1 = gvec2; // invoke *deep* copy of all contained groups
        CHECK(gvec1[0].begin() != gvec2[0].begin());
        CHECK(p1.front().id == 21);

        // new groups point to same particles as original
        auto gvec3 = gvec1;
        CHECK((*gvec1[0].begin()).id == (*gvec3[0].begin()).id);
    }

    SUBCASE("cerial serialisation") {
        std::ostringstream out(std::stringstream::binary);
        { // serialize g2
            std::vector<Particle> p2(5);
            Group<Particle> g2(p2.begin(), p2.end());
            p2.front().id = 8;
            p2.back().pos.x() = -10;
            g2.id = 100;
            g2.atomic = true;
            g2.compressible = true;
            g2.cm = {1, 0, 0};
            g2.confid = 20;
            g2.resize(4);
            cereal::BinaryOutputArchive archive(out);
            archive(g2);
        }

        {                                     // deserialize into g1
            std::istringstream in(out.str()); // not pretty...
            cereal::BinaryInputArchive archive(in);
            std::vector<Particle> p1(5);
            Group<Particle> g1(p1.begin(), p1.end());
            archive(g1);

            CHECK(g1.id == 100);
            CHECK(g1.atomic == true);
            CHECK(g1.compressible == true);
            CHECK(g1.cm.x() == 1);
            CHECK(g1.confid == 20);
            CHECK(g1.size() == 4);
            CHECK(g1.capacity() == 5);
            CHECK(g1.begin()->id == 8);
            CHECK(p1.front().id == 8);
            CHECK(p1.back().pos.x() == -10);
            CHECK(p1.back().ext == nullptr);
        }
    }
}

TEST_SUITE_END();

} // namespace Faunus
