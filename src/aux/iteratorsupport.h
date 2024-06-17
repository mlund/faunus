#pragma once
#include <iterator>
#include <range/v3/iterator/operations.hpp>

namespace Faunus {

template <class T1, class T2> auto distance(T1 first, T2 last)
{
    return std::distance(&(*first), &(*last));
} //!< Distance between two arbitrary contiguous iterators

template <class Trange> size_t range_size(Trange& rng)
{
    return ranges::cpp20::distance(rng.begin(), rng.end());
} //!< Size of a Ranges object using the `cpp20::distance()`

} // namespace Faunus
