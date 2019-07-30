#pragma once

#include <string>
#include <nlohmann/json_fwd.hpp>

namespace exprtk { // exprtk.hpp
template <typename T> class parser;
template <typename T> class expression;
template <typename T> class symbol_table;
}

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <nlohmann/json.hpp>
#else
#include <nlohmann/json_fwd.hpp>
#endif

/**
 * Since parser<T> is non-copyable we instantiate it
 * with a shared pointer, allowing `ExprFunction`
 * to be directly assigned to `std::function`.
 */
template<typename T=double>
class ExprFunction {
    std::shared_ptr<exprtk::parser<T>> parser;
    std::shared_ptr<exprtk::expression<T>> expression;
    std::shared_ptr<exprtk::symbol_table<T>> symbols;
    typedef std::vector<std::pair<std::string, T*>> Tvarvec;
    typedef std::vector<std::pair<std::string, T>> Tconstvec;

  public:
    void set(const std::string &exprstr, const Tvarvec &vars = {}, const Tconstvec &consts = {});
    void set(const nlohmann::json &, const Tvarvec &vars = {});
    T operator()() const;
};

extern template class ExprFunction<double>;

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] ExprFunction") {
    double x = 0, y = 0;
    ExprFunction<double> expr;
    nlohmann::json j = R"({ "function": "x*x+kappa", "constants": {"kappa": 0.4, "f": 2} })"_json;
    expr.set(j, {{"x", &x}, {"y", &y}});
    std::function<double()> f = expr;
    x = 4;
    CHECK( f() == doctest::Approx(4*4+0.4) );
}
#endif
