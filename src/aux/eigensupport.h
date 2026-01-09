#pragma once

#include <iterator>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

// Eigen<->JSON (de)serialization
namespace Eigen {

template <typename T> void to_json(nlohmann::json& j, const T& p)
{
    auto d = p.data();
    for (int i = 0; i < (int)p.size(); ++i)
        j.push_back(d[i]);
}

template <class T> void from_json(const nlohmann::json& j, Eigen::Matrix<T, 3, 1>& p)
{
    if (j.size() == 3) {
        int i = 0;
        for (auto d : j.get<std::vector<T>>())
            p[i++] = d;
        return;
    }
    throw std::runtime_error("JSON->Eigen conversion error");
}
} // namespace Eigen

namespace Faunus {

/**
 * @brief Eigen::Map facade to data members in STL container
 *
 * No data is copied and modifications of the Eigen object
 * modifies the original container and vice versa.
 *
 * Example:
 *
 *    std::vector<Tparticle> v(10);
 *    auto m1 = asEigenVector(v.begin, v.end(), &Tparticle::pos);    --> 10x3 maxtrix view
 *    auto m2 = asEigenMatrix(v.begin, v.end(), &Tparticle::charge); --> 10x1 vector view
 *
 * @warning Be careful that objects are properly aligned and divisible with `sizeof<double>`
 */
template <typename dbl = double, class iter, class memberptr>
auto asEigenMatrix(iter begin, iter end, memberptr m)
{
    using T = typename std::iterator_traits<iter>::value_type;
    static_assert(sizeof(T) % sizeof(dbl) == 0, "value_type size must multiples of double");
    constexpr size_t s = sizeof(T) / sizeof(dbl);
    constexpr size_t cols = sizeof((static_cast<T*>(0))->*m) / sizeof(dbl);
    using Tmatrix = Eigen::Matrix<dbl, Eigen::Dynamic, cols>;
    return Eigen::Map<Tmatrix, 0, Eigen::Stride<1, s>>((dbl*)&(*begin.*m), end - begin, cols)
        .array();
}

template <typename dbl = double, class iter, class memberptr>
auto asEigenVector(iter begin, iter end, memberptr m)
{
    using T = typename std::iterator_traits<iter>::value_type;
    static_assert(std::is_same<dbl&, decltype((static_cast<T*>(0))->*m)>::value,
                  "member must be a scalar");
    return asEigenMatrix<dbl>(begin, end, m).col(0);
}

} // namespace Faunus
