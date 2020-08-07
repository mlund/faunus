#pragma once

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include "aux/eigen_cerealisation.h"

// Eigen<->JSON (de)serialization
namespace Eigen {

template<typename T>
void to_json(nlohmann::json& j, const T &p) {
    auto d = p.data();
    for (int i=0; i<(int)p.size(); ++i)
        j.push_back( d[i] );
}

template<class T>
void from_json(const nlohmann::json& j, Eigen::Matrix<T,3,1> &p) {
    if ( j.size()==3 ) {
        int i=0;
        for (auto d : j.get<std::vector<T>>())
            p[i++] = d;
        return;
    }
    throw std::runtime_error("JSON->Eigen conversion error");
}
} // Eigen namespace
