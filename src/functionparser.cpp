#include <doctest/doctest.h>
#include "functionparser.h"
#include <exprtk.hpp> // https://github.com/ArashPartow/exprtk
#include <nlohmann/json.hpp>

template <std::floating_point T>
void ExprFunction<T>::set(const std::string& exprstr, const Tvarvec& vars, const Tconstvec& consts)
{
    if (not parser) {
        parser = std::make_shared<exprtk::parser<T>>();
        symbols = std::make_shared<exprtk::symbol_table<T>>();
        expression = std::make_shared<exprtk::expression<T>>();
    }
    symbols->clear();
    for (auto& v : vars)
        symbols->add_variable(v.first, *v.second);
    for (auto& v : consts)
        symbols->add_constant(v.first, v.second);
    symbols->add_constants();
    expression->register_symbol_table(*symbols);
    if (!parser->compile(exprstr, *expression))
        throw std::runtime_error("error passing function/expression");
}

template <std::floating_point T>
void ExprFunction<T>::set(const nlohmann::json& j, const Tvarvec& vars)
{
    Tconstvec consts;
    auto it = j.find("constants");
    if (it != j.end())
        for (auto i = it->begin(); i != it->end(); ++i)
            consts.push_back({i.key(), i.value()});
    set(j.at("function"), vars, consts);
}

template <std::floating_point T> T ExprFunction<T>::operator()() const
{
    return expression->value();
}

template <std::floating_point T> T ExprFunction<T>::derivative(T& variable) const
{
    return exprtk::derivative(*expression, variable);
}

template class ExprFunction<double>;

TEST_CASE("[Faunus] ExprFunction")
{
    double x = 0, y = 0;
    ExprFunction<double> expr;
    nlohmann::json j = R"({ "function": "x*x+kappa", "constants": {"kappa": 0.4, "f": 2} })"_json;
    expr.set(j, {{"x", &x}, {"y", &y}});
    std::function<double()> f = expr;
    x = 4;
    CHECK(f() == doctest::Approx(4 * 4 + 0.4));
    CHECK(expr.derivative(x) == doctest::Approx(2.0 * 4.0));
}
