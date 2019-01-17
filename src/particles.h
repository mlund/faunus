#pragma once

#include <vector>
#include <cassert>
#include <range/v3/all.hpp>

namespace Faunus {

    template<class T>
        class VectorDormantBasic {
            protected:
                typedef std::vector<T> Tvec;
                typedef typename Tvec::iterator Titerator;
                std::vector<int> active, inactive;
                Tvec vec;
            private:
                virtual void update() {} //!< invoked after each call to `move` or `push_back`

                void move( const Titerator &begin, const Titerator &end, std::vector<int> &src, std::vector<int> &dst ) {
                    assert( std::distance(vec.begin(), begin)>=0 && "Invalid iterators -- perhaps memory has been reallocated?");

                    // generate index interval in vec
                    std::vector<int> v( std::distance( begin,end ) );
                    std::iota(v.begin(), v.end(), std::distance( vec.begin(), begin ) );

                    // find and remove from source list
                    auto it = std::search( src.begin(), src.end(), v.begin(), v.end() ); 
                    assert(it!=src.end());
                    src.erase( it, it+v.size() );

                    // insert into sorted destination list
                    it = std::upper_bound( dst.begin(), dst.end(), v.front() );
                    dst.insert( it, v.begin(), v.end() );

                    assert( vec.size() == active.size() + inactive.size() );
                    assert( std::is_sorted(dst.begin(), dst.end() ) );
                    assert( std::is_sorted(src.begin(), src.end() ) );

                    update();
                } // move index from one vector to another while keeping them sorted

            public:
                Tvec& getVector() { return vec; } //!< Get full vector w. active and inactive particles
                const Tvec& getVector() const { return vec; } //!< Get full vector w. active and inactive particles

                void reserve(size_t size) {
                    vec.reserve( size );
                    active.reserve( size );
                    inactive.reserve( size );
                } //!< Allocate memory for elements

                int push_back(const Tvec &in) {
                    reserve( vec.size() + in.size() );

                    size_t index = vec.size();
                    for (size_t i=0; i<in.size(); i++)
                        active.push_back(index++);

                    vec.insert( vec.end(), in.begin(), in.end() ); 

                    assert( vec.size() == active.size() + inactive.size() );
                    assert( std::is_sorted( active.begin(), active.end() ) );

                    update();

                    return vec.size() - in.size();
                } //!< Insert elements at back

                void deactivate( const Titerator &begin, const Titerator &end ) {
                    move( begin, end, active, inactive );
                } //!< Deactivate elements

                void activate( const Titerator &begin, const Titerator &end ) {
                    move( begin, end, inactive, active );
                } //!< Activate elements

                size_t size() const { return active.size(); } //!< Size of active range
        };

    /**
     * @brief STL-like vector container with range support and element (de)activation
     *
     * This manages a private `std::vector` and has the following properties:
     *
     * - Elements can be added only to the end
     * - Elements cannot be removed, but _deactivated_
     * - An internal `std::vector<int>` keeps track of active/inactive elements
     * - The class acts as an iterable range pointing only to _active_ elements
     * - Range functionality is provided by the range-V3 project
     *
     * Example:
     *
     *     VectorDormant<int> v;
     *     v.push_back( {10,20,30,40} );
     *     v.size(); // --> 4
     *     auto it = v.getVector().begin();
     *     v.deactivate(it+1, it+2);
     *     v.size(); // --> 2
     *     for (auto i : v):
     *        cout << i << " "; // --> 10 40
     */
    template<class T>
        class VectorDormant : public ranges::view_facade<VectorDormant<T>>, public VectorDormantBasic<T>
    {
        private:
            friend ranges::range_access;
            using typename VectorDormantBasic<T>::Tvec;
            typename std::vector<int>::const_iterator _it, _end;
            typename Tvec::iterator _data;

            T& read() const { return *(_data+(*_it)); }
            bool equal(ranges::default_sentinel) const { return _it == _end; }
            void next() { ++_it; }

            void update() override {
                _it = this->active.begin();
                _end = this->active.end();
                _data = this->vec.begin();
            }

        public:
            size_t size() const { return VectorDormantBasic<T>::size(); }
    };

    struct Group : public Range::range {
        Point cm;     //!< Mass center
        int id;       //!< Group ID
        bool atomic;  //!< True if container species are atomic
        Group() {}
        Group(int beg, int end) : Range::range(beg, end+1, 1) {}

        int front() const {
            assert( !empty() );
            return *begin();
        }
        int back() const {
            assert( !empty() );
            return *(end()-1);
        }

        int random(Random &r) const {
            assert(!empty());
            return *r.element(begin(), end()); 
        } //!< Random index
    };


}//namespace
