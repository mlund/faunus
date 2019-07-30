#pragma once
#include "core.h"
#include "particle.h"
#include "molecule.h"
#include "geometry.h"
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <nlohmann/json_fwd.hpp>

namespace Faunus {

    template<typename T>
        void inline swap_to_back(T first, T last, T end) {
            while (end-- > last)
                std::iter_swap(first++,end);
        } //!< Move range [first:last] to [end] by swapping elements

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

            Group(Group &o);
            Group(const Group &o);
            Group(iter begin, iter end); //!< Constructor

            Group& operator=(const Group &o); //!< Deep copy contents from another Group

            Group& shallowcopy(const Group &o); //!< copy group data from `other` but *not* particle data

            bool contains(const T &a, bool include_inactive=false) const; //!< Determines if particle belongs to group (complexity: constant)

            auto find_id(int id) const {
                //return Faunus::filter(begin(), end(), [id](T &i){return (i.id==id);} );
                return *this | ranges::view::filter( [id](T &i){ return (i.id==id); } );
            } //!< Range of all (active) elements with matching particle id

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

            double mass() const; //!< Sum of all active masses

            std::vector<std::reference_wrapper<Point>> positions() const; //!< Iterable range with positions of active particles

            // warning! this should be tested -- do not use yet.
            template<typename TdistanceFunc>
                void unwrap(const TdistanceFunc &vdist) {
                    for (auto &i : *this)
                        i.pos = cm + vdist( i.pos, cm );
                } //!< Remove periodic boundaries with respect to mass center (Order N complexity).

            void wrap(Geometry::BoundaryFunction); //!< Apply periodic boundaries (Order N complexity).

            void translate(const Point&, Geometry::BoundaryFunction=[](Point&){}); //!< Translate particle positions and mass center

            void rotate(const Eigen::Quaterniond&, Geometry::BoundaryFunction); //!< Rotate all particles in group incl. internal coordinates (dipole moment etc.)

        }; //!< End of Group struct

        // Group<Particle> is instantiated elsewhere (group.cpp)
    extern template struct Group<Particle>;

void to_json(json&, const Group<Particle>&);
void from_json(const json&, Group<Particle>&);

}//end of faunus namespace
