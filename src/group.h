#pragma once

#include "core.h"
#include "particle.h"
#include "molecule.h"
#include "geometry.h"
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>

#ifdef DOCTEST_LIBRARY_INCLUDED
#include "units.h"
#include <Eigen/Geometry>
#include <range/v3/view/sample.hpp>
#include <range/v3/view/bounded.hpp>
#endif

namespace Faunus {

    template<typename T>
        void inline swap_to_back(T first, T last, T end) {
            while (end-- > last)
                std::iter_swap(first++,end);
        } //!< Move range [first:last] to [end] by swapping elements

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] swap_to_back") {
        typedef std::vector<int> _T;
        _T v = {1,2,3,4};

        swap_to_back( v.begin(), v.end(), v.end() );
        CHECK( v==_T({1,2,3,4}) );

        std::sort(v.begin(), v.end());
        swap_to_back( v.begin()+1, v.begin()+3, v.end() );
        CHECK( v==_T({1,4,3,2}) );
    }
#endif

    template<class T>
        struct IterRange : std::pair<T,T> {
            using std::pair<T,T>::pair;
            T& begin() { return this->first; }
            T& end() { return this->second; }
            const T& begin() const { return this->first; }
            const T& end() const { return this->second; }
            size_t size() const { return std::distance(this->first, this->second); }
            void resize(size_t n) { end() += n-size(); assert(size()==n); }
            bool empty() const { return this->first==this->second; }
            void clear() { this->second = this->first; assert(empty()); }
            std::pair<int,int> to_index(T reference) {
                return {std::distance(reference, begin()), std::distance(reference, end()-1)};
            } //!< Returns index pair relative to reference
        }; //!< Turns a pair of iterators into a range

    /**
     * @brief Turns a pair of iterators into an elastic range
     *
     * The elastic range is a range where elements can be deactivated
     * and later activated without inserting/erasing:
     *
     * - Just deactivated elements are moved to `end()` and can be retrieved from there.
     * - Just activated elements are placed at `end()-n`.
     * - The true size is given by `capacity()`
     */
    template<class T>
        class ElasticRange : public IterRange<typename std::vector<T>::iterator> {
            public:
                typedef typename std::vector<T>::iterator Titer;
            private:
                Titer _trueend;
            public:
                typedef IterRange<typename std::vector<T>::iterator> base;
                using base::begin;
                using base::end;
                using base::size;

                ElasticRange(Titer begin, Titer end);
                size_t capacity() const;
                auto inactive() const; //!< Range of inactive elements
                void deactivate(Titer first,
                                Titer last); //!< Deactivate particles by moving to end, reducing the effective size
                void activate(Titer first, Titer last); //!< Activate previously deactivated elements
                Titer &trueend();
                const Titer &trueend() const;
                void relocate(
                    Titer oldorigin,
                    Titer neworigin); //!< Shift all iterators to new underlying container; useful when resizing vectors
        };

        template <class T>
        ElasticRange<T>::ElasticRange(ElasticRange::Titer begin, ElasticRange::Titer end)
            : base({begin, end}), _trueend(end) {}

        template <class T> size_t ElasticRange<T>::capacity() const { return std::distance(begin(), _trueend); }

        template <class T> auto ElasticRange<T>::inactive() const { return base({end(), _trueend}); }

        template <class T> void ElasticRange<T>::deactivate(ElasticRange::Titer first, ElasticRange::Titer last) {
            size_t n = std::distance(first, last);
            assert(n >= 0);
            assert(first >= begin() && last <= end());
            std::rotate(begin(), last, end());
            end() -= n;
            assert(size() + inactive().size() == capacity());
        }

        template <class T> void ElasticRange<T>::activate(ElasticRange::Titer first, ElasticRange::Titer last) {
            size_t n = std::distance(first, last);
            std::rotate(end(), first, _trueend);
            end() += n;
            assert(size() + inactive().size() == capacity());
        }

        template <class T> typename ElasticRange<T>::Titer &ElasticRange<T>::trueend() { return _trueend; }

        template <class T> const typename ElasticRange<T>::Titer &ElasticRange<T>::trueend() const { return _trueend; }

        template <class T>
        void ElasticRange<T>::relocate(ElasticRange::Titer oldorigin, ElasticRange::Titer neworigin) {
            begin() = neworigin + std::distance(oldorigin, begin());
            end() = neworigin + std::distance(oldorigin, end());
            trueend() = neworigin + std::distance(oldorigin, trueend());
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] ElasticRange") {
        std::vector<int> v = {10,20,30,40,50,60};
        ElasticRange<int> r(v.begin(), v.end());
        CHECK( r.size() == 6 );
        CHECK( r.empty() == false );
        CHECK( r.size() == r.capacity() );
        *r.begin() += 1;
        CHECK( v[0]==11 );

        r.deactivate( r.begin(), r.end() );
        CHECK( r.size() == 0 );
        CHECK( r.empty() == true );
        CHECK( r.capacity() == 6 );
        CHECK( r.inactive().size() == 6 );
        CHECK( r.begin() == r.end() );

        r.activate( r.inactive().begin(), r.inactive().end() );
        CHECK( r.size() == 6);
        CHECK( std::is_sorted(r.begin(), r.end() ) == true ); // back to original

        r.deactivate( r.begin()+1, r.begin()+3 );
        CHECK( r.size() == 4 );
        CHECK( std::find(r.begin(), r.end(), 20)==r.end() );
        CHECK( std::find(r.begin(), r.end(), 30)==r.end() );
        CHECK( *r.end()==20); // deactivated elements can be retrieved from `end()`
        CHECK( *(r.end()+1)==30);

        auto ipair = r.to_index( v.begin() );
        CHECK( ipair.first == 0);
        CHECK( ipair.second == 3);

        r.activate( r.end(), r.end()+2 );
        CHECK( *(r.end()-2)==20); // activated elements can be retrieved from `end()-n`
        CHECK( *(r.end()-1)==30);
        CHECK( r.size() == 6);

        // check relocation
        auto v2 = v;
        v2.front()=-7;
        CHECK( *r.begin()!=-7 );
        r.relocate(v.begin(), v2.begin());
        CHECK( *r.begin()==-7 );
    }
#endif

    template<class T /** Particle type */>
        struct Group : public ElasticRange<T> {
            typedef ElasticRange<T> base;
            typedef typename base::Titer iter;
            typedef typename std::vector<T> Tpvec;
            using base::begin;
            using base::end;
            using base::size;
            int id=-1;           //!< Molecule id
            int confid=0;        //!< Conformation index / id
            Point cm={0,0,0};    //!< Mass center
            bool atomic=false;   //!< Is it an atomic group?

            const auto &traits() const { return molecules.at(id); } //!< Convenient access to molecule properties

            template<class Trange>
                Group(Trange &rng) : base(rng.begin(), rng.end()) {
                    assert(1==2);
                    // WARNING: Only iterators are copied!
                } //!< Constructor from range. WARNING: Only iterators are copied

            Group(Group &o) : base(o.begin(), o.trueend()) { *this = operator=(o); }
            Group(const Group &o) : base(o.begin(), o.trueend()) { *this = operator=(o); }
            Group(iter begin, iter end) : base(begin,end) {} //!< Constructor

            Group& operator=(const Group &o) {
                if (&o == this)
                    return *this;
                shallowcopy(o);
                if (o.begin()!=begin())
                    std::copy(o.begin(), o.trueend(), begin()); // copy all particle data
                return *this;
            } //!< Deep copy contents from another Group

            Group& shallowcopy(const Group &o) {
                if (&o != this) {
                    if (this->capacity() != o.capacity())
                        throw std::runtime_error("Group::shallowcopy: capacity mismatch");
                    this->resize(o.size());
                    id = o.id;
                    atomic = o.atomic;
                    cm = o.cm;
                    confid = o.confid;
                }
                return *this;
            } //!< copy group data from `other` but *not* particle data

            bool contains(const T &a, bool include_inactive=false) const {
                if (not this->empty()) {
                    int size = (include_inactive ? this->capacity() : this->size());
                    int d = &a - &(*(this->begin()));
                    if (d>=0 and d<size)
                        return true;
                }
                return false;
            } //!< Determines if particle belongs to group (complexity: constant)

            template<class UnaryPredicate>
                auto filter( UnaryPredicate f ) const {
                    return ranges::view::filter(*this,f);
                } //!< Filtered range according to unary predicate, `f`

            auto find_id(int id) const {
                return *this | ranges::view::filter( [id](T &i){ return (i.id==id); } );
            } //!< Range of all elements with matching particle id

            /*
             * @brief Reference to subset of given index, where 0 is the start of the group
             * @note do not parse index as `const&` which would create a dangling reference
             */
            template<class Tint=size_t>
                auto operator[](std::vector<Tint> &index) {
#ifndef NDEBUG
                    // check that range is within group
                    if (not index.empty()) 
                        assert( *std::max_element(index.begin(), index.end()) < size() );
#endif
                    return index | ranges::view::transform([this](Tint i) -> T& { return *(this->begin()+i); });
                }

            double mass() const {
                double m=0;
                for (auto &i : *this)
                    m += atoms[i.id].mw;
                return m;
            } //!< Sum of all active masses

            auto positions() const {
                return ranges::view::transform(*this, [](auto &i) -> Point& {return i.pos;});
            } //!< Iterable range with positions of active particles

            // warning! this should be tested -- do not use yet.
            template<typename TdistanceFunc>
                void unwrap(const TdistanceFunc &vdist) {
                    for (auto &i : *this)
                        i.pos = cm + vdist( i.pos, cm );
                } //!< Remove periodic boundaries with respect to mass center (Order N complexity).

            void wrap(Geometry::BoundaryFunction boundary) {
                boundary(cm);
                for (auto &i : *this)
                    boundary(i.pos);
            } //!< Apply periodic boundaries (Order N complexity).

            void translate(const Point &d, Geometry::BoundaryFunction boundary=[](Point&){}) {
                cm += d;
                boundary(cm);
                for (auto &i : *this) {
                    i.pos += d;
                    boundary(i.pos);
                }
            } //!< Translate particle positions and mass center

            void rotate(const Eigen::Quaterniond &Q, Geometry::BoundaryFunction boundary) {
                Geometry::rotate(begin(), end(), Q, boundary, -cm);
            } //!< Rotate all particles in group incl. internal coordinates (dipole moment etc.)

        }; //!< Groups of particles

    template<class T /** Particle type */>
        void to_json(json &j, const Group<T> &g) {
            j = {
                {"id", g.id}, {"cm", g.cm}, {"atomic", g.atomic}, {"size", g.size()}
            };
            if (g.capacity()>g.size())
                j["capacity"] = g.capacity();
            if (g.confid!=0)
                j["confid"] = g.confid;
        }

    template<class T /** Particle type */>
        void from_json(const json &j, Group<T> &g) {
            g.resize( j.at("size").get<int>() );
            g.trueend() = g.begin() + j.value("capacity", g.size());
            g.id = j.at("id").get<unsigned int>();
            g.cm = j.at("cm").get<Point>();
            g.atomic = j.at("atomic").template get<bool>();
            g.confid = j.value("confid", 0);
        }

    template<class Trange>
        auto positions(Trange &r) {
            return ranges::view::transform(r, [](auto &i) -> Point& {return i.pos;});
        } //!< Iterable range with positions (works for groups and particle vectors)

extern template struct Group<Particle>;

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Group") {
        Random rand;
        std::vector<Particle> p(3);
        p.reserve(10);
        p[0].id=0;
        p[1].id=1;
        p[2].id=1;
        Group<Particle> g(p.begin(), p.end());

        SUBCASE("contains()") {
            CHECK( g.contains(p[0]) );
            CHECK( g.contains(p[1]) );
            CHECK( g.contains(p[2]) );
            CHECK( g.contains(p[3]) == false );
        }
 
        // find all elements with id=1
        auto slice1 = g.find_id(1);
        CHECK( std::distance(slice1.begin(), slice1.end()) == 2 );

        // find *one* random value with id=1
        auto slice2 = slice1 | ranges::view::sample(1, rand.engine) | ranges::view::bounded;
        CHECK( std::distance(slice2.begin(), slice2.end()) == 1 );

        // check rotation
        Eigen::Quaterniond q;
        q = Eigen::AngleAxisd(pc::pi/2, Point(1,0,0));
        p.at(0).pos = p.at(0).getExt().mu = p.at(0).getExt().scdir = {0, 1, 0};

        Geometry::Chameleon geo = R"({"type":"cuboid", "length": [2,2,2]})"_json;
        g.rotate(q, geo.getBoundaryFunc());
        CHECK( p[0].pos.y() == doctest::Approx(0) );
        CHECK( p[0].pos.z() == doctest::Approx(1) );
        CHECK(p[0].getExt().mu.y() == doctest::Approx(0));
        CHECK(p[0].getExt().mu.z() == doctest::Approx(1));
        CHECK(p[0].getExt().scdir.y() == doctest::Approx(0));
        CHECK(p[0].getExt().scdir.z() == doctest::Approx(1));

        p[0].pos = {1,2,3};
        p[1].pos = {4,5,6};

        // iterate over positions and modify them
        for (Point &i : g.positions() )
            i = 2*i;
        CHECK( p[1].pos.x() == doctest::Approx(8) );
        CHECK( p[1].pos.y() == doctest::Approx(10) );
        CHECK( p[1].pos.z() == doctest::Approx(12) );

        SUBCASE("operator[]") {
            CHECK( p.begin() == g.begin() );
            CHECK( p.end() == g.end() );

            // a new range by using an index filter
            std::vector<size_t> index = {0,1};
            auto subset = g[index];
            CHECK( subset.size()==2 );
            CHECK( &(*p.begin()) == &(*subset.begin()) );
            CHECK( &(*(p.begin()+1)) == &(*(subset.begin()+1)) );
            for (auto& i : subset)
                i.pos *= 2;
            CHECK( p[1].pos.x() == doctest::Approx(16) );
            CHECK( p[1].pos.y() == doctest::Approx(20) );
            CHECK( p[1].pos.z() == doctest::Approx(24) );
        }

        // check deep copy and resizing
        std::vector<int> p1(5), p2(5);
        p1.front() = 1;
        p2.front() = -1;

        Group<int> g1(p1.begin(), p1.end());
        Group<int> g2(p2.begin(), p2.end());

        g2.id=100;
        g2.atomic=true;
        g2.cm={1,0,0};
        g2.confid=20;
        g1=g2;

        CHECK(g1.id==100);
        CHECK(g1.atomic==true);
        CHECK(g1.cm.x()==1);
        CHECK(g1.confid==20);

        CHECK( *g1.begin()==-1);
        CHECK( *g2.begin()==-1);
        CHECK( g1.begin() != g2.begin() );
        CHECK( g1.size() == g2.size());
        *g2.begin()=10;
        g2.resize(4);
        g1 = g2;
        CHECK( g1.size() == 4);
        CHECK( g1.capacity() == 5);
        CHECK( p1.front()==10);

        std::vector<Group<int>> gvec1, gvec2;
        gvec1.push_back(g1);
        gvec2.push_back(g2);
        p2.front()=21;

        CHECK( *(gvec1.front().begin())==10 );
        CHECK( *(gvec2.front().begin())==21 );

        // existing groups point to existing particles when overwritten
        gvec1 = gvec2; // invoke *deep* copy of all contained groups
        CHECK( gvec1[0].begin() != gvec2[0].begin());
        CHECK( p1.front() == 21);

        // new groups point to same particles as original
        auto gvec3 = gvec1;
        CHECK( gvec1[0].begin() == gvec3[0].begin());
    }
#endif

}//end of faunus namespace
