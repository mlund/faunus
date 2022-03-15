#pragma once

#include <string>
#include <nlohmann/json.hpp>

namespace exprtk { // exprtk.hpp
template <typename T> class parser;
template <typename T> class expression;
template <typename T> class symbol_table;
}

/**
 * Since parser<T> is non-copyable we instantiate it
 * with a shared pointer, allowing `ExprFunction`
 * to be directly assigned to `std::function`.
 */
template <std::floating_point T = double> class ExprFunction {
    std::shared_ptr<exprtk::parser<T>> parser;
    std::shared_ptr<exprtk::expression<T>> expression;
    std::shared_ptr<exprtk::symbol_table<T>> symbols;
    typedef std::vector<std::pair<std::string, T*>> Tvarvec;
    typedef std::vector<std::pair<std::string, T>> Tconstvec;

  public:
    void set(const std::string &exprstr, const Tvarvec &vars = {}, const Tconstvec &consts = {});
    void set(const nlohmann::json &, const Tvarvec &vars = {});
    T operator()() const;
    T derivative(T& variable) const;
};

extern template class ExprFunction<double>;

