#pragma once
#include "core.h"
#include "particle.h"
#include "molecule.h"
#include "geometry.h"
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <nlohmann/json.hpp>
#include <optional>

namespace Faunus {

template <typename T> void inline swap_to_back(T first, T last, T end) {
    while (end-- > last) {
        std::iter_swap(first++, end);
    }
} //!< Move range [first:last] to [end] by swapping elements

template <class T> struct IterRange : std::pair<T, T> {
    using std::pair<T, T>::pair;
    T& begin() { return this->first; }
    T& end() { return this->second; }
    const T& begin() const { return this->first; }
    const T& end() const { return this->second; }
    size_t size() const { return std::distance(this->first, this->second); }
    void resize(size_t n) {
        end() += n - size();
        assert(size() == n);
    }
    bool empty() const { return this->first == this->second; }
    void clear() {
        this->second = this->first;
        assert(empty());
    }
    std::pair<int, int> to_index(T reference) const {
        return {std::distance(reference, begin()), std::distance(reference, end() - 1)};
    } //!< Returns particle index pair relative to given reference
};    //!< Turns a pair of iterators into a range

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
template <class T> class ElasticRange : public IterRange<typename std::vector<T>::iterator> {
  public:
    using Titer = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

  private:
    Titer _trueend;

  public:
    using base = IterRange<typename std::vector<T>::iterator>;
    using base::begin;
    using base::end;
    using base::size;

    ElasticRange(Titer begin, Titer end);
    size_t capacity() const;
    auto inactive() const; //!< Range of inactive elements
    void deactivate(Titer first,
                    Titer last);            //!< Deactivate particles by moving to end, reducing the effective size
    void activate(Titer first, Titer last); //!< Activate previously deactivated elements
    Titer& trueend();
    const Titer& trueend() const;
    void relocate(const_iterator oldorigin,
                  Titer neworigin); //!< Shift all iterators to new underlying container; useful when resizing vectors
    bool isFull() const { return end() == trueend(); }
    auto numInactive() const; //!< Number of inactive elements
};

template <class T>
ElasticRange<T>::ElasticRange(ElasticRange::Titer begin, ElasticRange::Titer end) : base({begin, end}), _trueend(end) {}

template <class T> size_t ElasticRange<T>::capacity() const { return std::distance(begin(), _trueend); }
template <class T> auto ElasticRange<T>::inactive() const { return base({end(), _trueend}); }

template <class T> void ElasticRange<T>::deactivate(ElasticRange::Titer first, ElasticRange::Titer last) {
    size_t n = std::distance(first, last);
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

template <class T> typename ElasticRange<T>::Titer& ElasticRange<T>::trueend() { return _trueend; }
template <class T> const typename ElasticRange<T>::Titer& ElasticRange<T>::trueend() const { return _trueend; }

template <class T>
void ElasticRange<T>::relocate(ElasticRange::const_iterator oldorigin, ElasticRange::Titer neworigin) {
    begin() = neworigin + std::distance(oldorigin, const_iterator(begin()));
    end() = neworigin + std::distance(oldorigin, const_iterator(end()));
    trueend() = neworigin + std::distance(oldorigin, const_iterator(trueend()));
}

template <class T> auto ElasticRange<T>::numInactive() const { return inactive().size(); }

class Group : public ElasticRange<Particle> {
  public:
    using base = ElasticRange<Particle>;
    using iter = typename base::Titer;
    int id = -1;                         //!< Molecule id
    int conformation_id = 0;             //!< Conformation index / id
    Point mass_center = {0.0, 0.0, 0.0}; //!< Mass center

    inline bool isAtomic() const { return traits().atomic; }     //!< Is it an atomic group?
    inline bool isMolecular() const { return !traits().atomic; } //!< is it a molecular group?

    std::optional<std::reference_wrapper<Point>> massCenter();             //!< Optional reference to mass center
    std::optional<std::reference_wrapper<const Point>> massCenter() const; //!< Optional ref. to mass center

    //! Selections to filter groups using `getSelectionFilter()`
    enum Selectors : unsigned int {
        ANY = (1U << 1U),       //!< Match any group (disregards all other flags)
        ACTIVE = (1U << 2U),    //!< Only active groups (non-zero size)
        INACTIVE = (1U << 3U),  //!< Only inactive groups (zero size)
        NEUTRAL = (1U << 4U),   //!< Only groups with zero net charge
        ATOMIC = (1U << 5U),    //!< Only atomic groups
        MOLECULAR = (1U << 6U), //!< Only molecular groups (atomic=false)
        FULL = (1U << 7U)       //!< Only groups where size equals capacity
    };

    /**
     * @brief Determines if given `Selectors` bitmask matches group
     * @tparam mask Bitmask build from enum `Group::Selectors`
     * @return true if ALL enabled bits in the mask are satisfied
     *
     * Note that for `INACTIVE | NEUTRAL`, the criterion is applied
     * to all active and inactive particles, i.e. until `trueend()`.
     */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++11-narrowing"
    template <unsigned int mask> bool match() const {
        static_assert(mask >= ANY && mask <= FULL);
        if constexpr (mask & ANY) {
            static_assert(mask == ANY, "don't mix ANY with other flags");
            return true;
        }
        if constexpr (mask & ACTIVE) {
            static_assert(!(mask & INACTIVE), "don't mix ACTIVE and INACTIVE");
            if (size() == 0) {
                return false;
            }
        } else if constexpr (mask & INACTIVE) {
            if (!empty()) {
                return false;
            }
        }
        if constexpr (mask & FULL) {
            if (end() != trueend()) {
                return false;
            }
        }
        if constexpr (mask & ATOMIC) {
            static_assert(!(mask & MOLECULAR), "don't mix ATOMIC and MOLECULAR");
            if (isMolecular()) {
                return false;
            }
        } else if constexpr (mask & MOLECULAR) {
            if (isAtomic()) {
                return false;
            }
        }
        if constexpr (mask & NEUTRAL) {
            auto _end = (mask & INACTIVE) ? trueend() : end();
            auto _charge = std::accumulate(begin(), _end, 0.0,
                                           [](auto sum, const auto& particle) { return sum + particle.charge; });
            if (std::fabs(_charge) > pc::epsilon_dbl) {
                return false;
            }
        }
        return true;
    }
#pragma clang diagnostic pop

    inline const MoleculeData& traits() const {
        assert(id >= 0 && id < Faunus::molecules.size());
        return Faunus::molecules[id];
    } //!< Convenient access to molecule properties

    Group(Group& other);
    Group(const Group& other);
    Group(MoleculeData::index_type molid, iter begin, iter end); //!< Constructor
    Group& operator=(const Group& other);                        //!< Deep copy contents from another Group
    Group& shallowCopy(const Group& other);                      //!< copy from `other` but *not* particle data
    bool contains(const Particle& particle, bool include_inactive = false) const; //!< Does particle belong?
    double mass() const;                                                          //!< Sum of all active masses
    auto positions(); //!< Range of positions of active particles

    AtomData::index_type getParticleIndex(
        const Particle& particle,
        bool include_inactive = false) const; //!< Finds index of particle within group. Throws if not part of group

    auto findAtomID(AtomData::index_type atomid) const {
        return *this | ranges::cpp20::views::filter([atomid](auto& particle) { return (particle.id == atomid); });
    } //!< Range of all (active) elements with matching particle id

    /**
     * @brief Returns i'th element in group
     * @param i index starting at zero
     * @return reference to value at i'th element
     * @note No range-checking and i must be in interval `[0:size[`
     */
    inline auto& operator[](size_t index) { return *(begin() + index); }
    inline const auto& operator[](size_t index) const { return *(begin() + index); }

    Particle& at(size_t index);
    const Particle& at(size_t index) const;

    /**
     * @brief Reference to subset of given indices, where 0 is the start of the group
     * @param indices Range of indices relative to group
     * @warning Do not change `indices` to `const&` which would create a dangling reference
     */
    template <class Tint = size_t> auto operator[](std::vector<Tint>& indices) {
        static_assert(std::is_integral<Tint>::value, "integer indices expected");
#ifndef NDEBUG
        if (not indices.empty()) {
            assert(*std::max_element(indices.begin(), indices.end()) < size());
        }
#endif
        return indices | ranges::cpp20::views::transform([this](auto i) -> Particle& { return *(begin() + i); });
    }

    /**
     * @brief Remove PBC for molecular groups w. respect to mass center
     * @tparam TdistanceFunc
     * @param vector_distance Distance calculation function
     * @remarks Atomic groups are not touched
     * @warning Is this fit for use?
     */
    template <typename TdistanceFunc> void unwrap(const TdistanceFunc& vector_distance) {
        if (isMolecular()) {
            for (auto& particle : *this) {
                particle.pos = mass_center + vector_distance(particle.pos, mass_center);
            }
        }
    }

    void wrap(Geometry::BoundaryFunction boundary); //!< Apply periodic boundaries (Order N complexity).

    void translate(
        const Point& displacement,
        Geometry::BoundaryFunction boundary = [](Point&) {}); //!< Translate particle positions and mass center

    void
    rotate(const Eigen::Quaterniond& quaternion,
           Geometry::BoundaryFunction boundary); //!< Rotate particles incl. internal coordinates (dipole moment etc.)

    void updateMassCenter(Geometry::BoundaryFunction boundary,
                          const Point& approximate_mass_center); //!< Calculates mass center

    void updateMassCenter(Geometry::BoundaryFunction boundary); //!< Calculates mass center

}; //!< End of Group class

void to_json(json& j, const Group& group);
void from_json(const json& j, Group& group);

//! Get lambda function matching given enum Select mask
template <unsigned int mask> std::function<bool(const Group&)> getGroupFilter() {
    return [](const Group& group) { return group.match<mask>(); };
}

/*
 * The following two functions are used to perform a complete
 * serialisation of a group to/from a Cereal archive. When loading,
 * the group *must* match the capacity of the saved data or an exception
 * is thrown.
 */

template <class Archive> void save(Archive& archive, const Group& group, std::uint32_t const version) {
    switch (version) {
    case 0:
        archive(group.id, group.conformation_id, group.mass_center, group.size(), group.capacity());
        std::for_each(group.begin(), group.trueend(), [&archive](auto& particle) { archive(particle); });
        break;
    default:
        throw std::runtime_error("unknown serialisation version");
    };
} //!< Cereal serialisation

template <class Archive> void load(Archive& archive, Group& group, std::uint32_t const version) {
    size_t size = 0;
    size_t capacity = 0;
    switch (version) {
    case 0:
        archive(group.id, group.conformation_id, group.mass_center, size, capacity);
        if (capacity != group.capacity()) {
            throw std::runtime_error("capacity mismatch of archived group");
        }
        group.resize(size);
        std::for_each(group.begin(), group.trueend(), [&archive](auto& particle) { archive(particle); });
        break;
    default:
        throw std::runtime_error("unknown serialisation version");
    }
} //!< Cereal serialisation

} // namespace Faunus
