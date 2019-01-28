#pragma once
#include <string>
#include <functional>
#include <iostream>
#include <exprtk.hpp> // https://github.com/ArashPartow/exprtk
#include <nlohmann/json.hpp>

/**
 * Since parser<T> is non-copyable we instantiate it
 * with a shared pointer, allowing `ExprFunction`
 * to be directly assigned to `std::function`.
 */
template <typename T=double> class ExprFunction {
    private:
        std::shared_ptr<exprtk::parser<T>> parser;
        exprtk::expression<T> expression;
        exprtk::symbol_table<T> symbols;
        typedef std::vector<std::pair<std::string,T*>> Tvarvec;
        typedef std::vector<std::pair<std::string,T>> Tconstvec;
    public:
        void set(const std::string &exprstr, const Tvarvec &vars={}, const Tconstvec &consts={}) {
            if (not parser)
                parser = std::make_shared<exprtk::parser<T>>();
            symbols.clear();
            for (auto &v : vars)
                symbols.add_variable(v.first,*v.second);
            for (auto &v : consts)
                symbols.add_constant(v.first,v.second);
            symbols.add_constants();
            expression.register_symbol_table(symbols);
            if (not parser->compile( exprstr, expression ))
                throw std::runtime_error("error passing function/expression");
        }

        void set(const nlohmann::json &j, const Tvarvec &vars={}) {
            Tconstvec consts;
            auto it = j.find("constants");
            if (it!=j.end())
                for (auto i=it->begin(); i!=it->end(); ++i)
                    consts.push_back( {i.key(), i.value()} );
            set(j.at("function"), vars, consts);
        }

        T operator()() const { return expression.value(); }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] ExprFunction") {
    double x=0, y=0;
    ExprFunction<double> expr;
    nlohmann::json j = R"( { "function":"x*x+kappa", "constants":{ "kappa":0.4, "f":2 } })"_json;
    expr.set(j, {{"x",&x}, {"y",&y}} );
    std::function<double()> f = expr;
    x=4;
    CHECK(f() == doctest::Approx(4*4+0.4));
}
#endif
