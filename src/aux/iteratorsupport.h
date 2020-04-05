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
 * @brief Iterator view to unique, self-avoiding pairs in a container
 *
 * The operation corresponds to iterating through a triangular matrix sans the trace.
 * Deferencing the iterator yields a tuple with (const) references to the
 * two data elements.
 *
 * Example:
 *
 * ~~~ cpp
 * std::vector<std::string> v = {"one","two", "three"};
 * for (auto [i,j] : internal_pairs(v)):
 *    cout << i << "-" << j << " ";
 * // --> "one-two" "one-three" "two-three"
 * ~~~
 *
 * @todo: implement non-const version
 */
template <class T, bool Const = std::is_const<T>::value> class internal_pairs {
  private:
    T &vec;
    struct iterator {
        // required for e.g. std::distance
        using iter = typename std::conditional<Const, typename T::const_iterator, typename T::iterator>::type;
        using pointer = void;
        using reference = typename std::conditional<Const, typename T::const_reference, typename T::reference>::type;
        using value_type = std::tuple<reference, reference>;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag; // we can only move forward

        iter _end, i, j;
        iterator(iter end, iter i, iter j) : _end(end), i(i), j(j){};
        inline value_type operator*() { return {*i, *j}; }
        inline bool operator!=(const iterator &other) const { return (i != other.i) or (j != other.j); }
        inline iterator &operator++() {
            if (++j == _end)
                j = std::next(++i);
            return *this;
        }
    };

  public:
    internal_pairs(T &vec) : vec(vec) {}

    template <bool _Const = Const> std::enable_if_t<_Const, iterator> begin() const {
        return iterator(vec.end(), vec.begin(), std::next(vec.begin()));
    } // first pair

    template <bool _Const = Const> std::enable_if_t<_Const, iterator> end() const {
        return iterator(vec.end(), std::prev(vec.end()), vec.end());
    } // one iteration after last pair

    template <bool _Const = Const> std::enable_if_t<not _Const, iterator> begin() {
        return iterator(vec.end(), vec.begin(), std::next(vec.begin()));
    } // first pair

    template <bool _Const = Const> std::enable_if_t<not _Const, iterator> end() {
        return iterator(vec.end(), std::prev(vec.end()), vec.end());
    } // one iteration after last pair

    size_t size() const { return vec.size() * (vec.size() - 1) / 2; }
};

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
