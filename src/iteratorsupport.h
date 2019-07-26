#pragma once
#include <iterator>
#include <range/v3/distance.hpp>

namespace Faunus {

template<class T1, class T2>
int distance(T1 first, T2 last) {
    return std::distance( &(*first), &(*last) );
} //!< Distance between two arbitrary contiguous iterators

template<class T>
int size(T &rng) {
    return ranges::distance(rng.begin(), rng.end());
} //!< Size of arbitrary range

}//end of namespace
