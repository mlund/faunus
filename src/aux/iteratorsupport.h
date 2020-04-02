#pragma once
#include <iterator>
#include <range/v3/distance.hpp>

namespace Faunus {

template<class T1, class T2>
int distance(T1 first, T2 last) {
    return std::distance( &(*first), &(*last) );
} //!< Distance between two arbitrary contiguous iterators

template <class Trange> size_t range_size(Trange &rng) {
    return ranges::cpp20::distance(rng.begin(), rng.end());
} //!< Size of a Ranges object using the `cpp20::distance()`

/**
 * Constant view to pairs in two containers
 * Calling `std::distance` is of O(N) complexity while `size` has constant complexity
 */
template <class Iter1, class Iter2> class cartesian_product {
  private:
    struct iterator {
        // these five are useful for stl
        using iterator_category = std::forward_iterator_tag; // we can only move forward
        using value_type = std::tuple<const typename Iter1::value_type &, const typename Iter2::value_type &>;
        using difference_type = std::ptrdiff_t;
        using reference = value_type &;
        using pointer = value_type *;

        Iter1 pos1, last1;
        Iter2 pos2, first2, last2;
        iterator(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2)
            : pos1(first1), last1(last1), pos2(first2), first2(first2), last2(last2) {}

        inline value_type operator*() { return {*pos1, *pos2}; }
        inline bool operator!=(const iterator &other) const { return (pos1 != other.pos1) or (pos2 != other.pos2); }
        inline iterator &operator++() {
            if (++pos2 == last2) {
                pos2 = first2;
                pos1++;
                if (pos1 == last1)
                    pos2 = last2;
            }
            return *this;
        }
    };
    iterator _begin, _end;

  public:
    cartesian_product(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2)
        : _begin(first1, last1, first2, last2), _end(last1, last1, last2, last2) {
        if (size() == 0)
            _begin = _end;
    }
    auto begin() const { return _begin; }
    auto end() const { return _end; }
    size_t size() const {
        return std::distance(_begin.pos1, _begin.last1) * std::distance(_begin.first2, _begin.last2);
    }
};

}//end of namespace
