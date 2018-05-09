#pragma once

#include "core.h"
#include "molecule.h"
#include "geometry.h"

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
                return {std::distance(reference, begin()), std::distance(reference, end())};
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

                ElasticRange(Titer begin, Titer end) : base({begin,end}), _trueend(end) {}

                size_t capacity() const { return std::distance(begin(), _trueend); }

                auto inactive() const {
                    return base({ end(), _trueend});
                } //!< Range of inactive elements

                void deactivate(Titer first, Titer last) {
                    size_t n = std::distance(first,last);
                    assert(n>=0);
                    assert(first>=begin() && last<=end() );
                    std::rotate( begin(), last, end() );
                    end() -= n;
                    assert(size() + inactive().size() == capacity());
                } //!< Deactivate particles by moving to end, reducing the effective size

                void activate(Titer first, Titer last) {
                    size_t n = std::distance(first,last);
                    std::rotate( end(), first, _trueend );
                    end() += n;
                    assert(size() + inactive().size() == capacity());
                } //!< Activate previously deactivated elements

                Titer& trueend() { return _trueend; }
                const Titer& trueend() const { return _trueend; }

                void relocate(Titer oldorigin, Titer neworigin) {
                    begin() = neworigin + std::distance(oldorigin, begin());
                    end() = neworigin + std::distance(oldorigin, end());
                    trueend() = neworigin + std::distance(oldorigin, trueend());
                } //!< Shift all iterators to new underlying container; useful when resizing vectors
        };

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
        CHECK( ipair.second == 4);

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
            using base::begin;
            using base::end;
            int id=-1;
            bool atomic=false;   //!< Is it an atomic group?
            Point cm={0,0,0};    //!< Mass center

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
                std::copy(o.begin(), o.end(), begin()); // copy all particle data
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
                }
                return *this;
            } //!< copy group data from `other` but *not* particle data

            bool contains(const T &a) const {
                int d = &a - &(*(this->begin()));
                if (d>=0)
                    if (d<this->size())
                        return true;
                return false;
            } //!< Determines if particle belongs to (active part) of group

            template<class UnaryPredicate>
                auto filter( UnaryPredicate f ) const {
                    return ranges::view::filter(*this,f);
                } //!< Filtered range according to unary predicate, `f`

            auto find_id(int id) const {
                return *this | ranges::view::filter( [id](T &i){ return (i.id==id); } );
            } //!< Range of all elements with matching particle id

            template<class Trange=std::vector<int>>
                auto find_index( const Trange &index ) {
                    return ranges::view::indirect(
                            ranges::view::transform(index, [this](int i){return this->begin()+i;}) );
                } //!< Group subset matching given `index` (Complexity: linear with index size)

            auto positions() const {
                return ranges::view::transform(*this, [](auto &i) -> Point& {return i.pos;});
            } //!< Iterable range with positions

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
                {"id", g.id}, {"cm", g.cm}, {"atomic", g.atomic}, {"size", g.size()},
                {"capacity", g.capacity()}
            };
        }

    template<class T /** Particle type */>
        void from_json(const json &j, Group<T> &g) {
            g.trueend() = g.begin() + j.at("capacity"); 
            g.resize( j.at("size") );
            g.id = j.at("id").get<unsigned int>();
            g.cm = j.at("cm").get<Point>();
            g.atomic = j.at("atomic").template get<bool>();
        }

    template<class Trange>
        auto positions(Trange &r) {
            return ranges::view::transform(r, [](auto &i) -> Point& {return i.pos;});
        } //!< Iterable range with positions (works for groups and particle vectors)

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Group") {
        Random rand;
        typedef ParticleAllProperties particle;
        std::vector<particle> p(3);
        p[0].id=0;
        p[1].id=1;
        p[2].id=1;
        Group<particle> g(p.begin(), p.end());

        // find all elements with id=1
        auto slice1 = g.find_id(1);
        CHECK( std::distance(slice1.begin(), slice1.end()) == 2 );

        // find *one* random value with id=1
        auto slice2 = slice1 | ranges::view::sample(1, rand.engine) | ranges::view::bounded;
        CHECK( std::distance(slice2.begin(), slice2.end()) == 1 );

        // check rotation
        Eigen::Quaterniond q;
        q = Eigen::AngleAxisd(pc::pi/2, Point(1,0,0));
        p.at(0).pos = p.at(0).mu = p.at(0).scdir = {0,1,0};

        Geometry::Cuboid geo = R"({"length": [2,2,2]})"_json;
        g.rotate(q, geo.boundaryFunc);
        CHECK( p[0].pos.y() == doctest::Approx(0) );
        CHECK( p[0].pos.z() == doctest::Approx(1) );
        CHECK( p[0].mu.y() == doctest::Approx(0) );
        CHECK( p[0].mu.z() == doctest::Approx(1) );
        CHECK( p[0].scdir.y() == doctest::Approx(0) );
        CHECK( p[0].scdir.z() == doctest::Approx(1) );

        p[0].pos = {1,2,3};
        p[1].pos = {4,5,6};

        // iterate over positions and modify them
        for (auto &i : g.positions() )
            i *= 2;
        CHECK( p[1].pos.x() == doctest::Approx(8) );
        CHECK( p[1].pos.y() == doctest::Approx(10) );
        CHECK( p[1].pos.z() == doctest::Approx(12) );

#pragma OPTIMIZE OFF 
        SUBCASE("test find_index") {
            CHECK( p.begin() == g.begin() );
            CHECK( p.end() == g.end() );

            // a new range by using an index filter
            auto subset = g.find_index( {0,1} );
            CHECK( std::distance(subset.begin(), subset.end()) == 2 );
            // On some compilers, this may FAIL if not in Debug mode
            //for (auto& i : subset) // On some compilers, this may FAIL if not in Debug mode
            //    i.pos *= 2;
            //CHECK( p[1].pos.x() == doctest::Approx(16) );
            //CHECK( p[1].pos.y() == doctest::Approx(20) );
            //CHECK( p[1].pos.z() == doctest::Approx(24) );
        }
#pragma OPTIMIZE ON

        // check deep copy and resizing
        std::vector<int> p1(5), p2(5);
        p1.front() = 1;
        p2.front() = -1;

        Group<int> g1(p1.begin(), p1.end());
        Group<int> g2(p2.begin(), p2.end());

        g2.id=100;
        g2.atomic=true;
        g2.cm={1,0,0};
        g1=g2;

        CHECK(g1.id==100);
        CHECK(g1.atomic==true);
        CHECK(g1.cm.x()==1);

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
