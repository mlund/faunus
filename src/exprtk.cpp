#include <string>
#include <functional>
#include <iostream>
#include "exprtk.hpp"
#include "json.hpp"

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
            if (!parser)
                parser = std::make_shared<exprtk::parser<T>>();
            symbols.clear();
            for (auto &v : vars) symbols.add_variable(v.first,*v.second);
            for (auto &v : consts) symbols.add_constant(v.first,v.second);
            symbols.add_constants();
            expression.register_symbol_table(symbols);
            parser->compile(exprstr,expression);
        }
        T operator()() { return expression.value(); }
#ifdef NLOHMANN_JSON_HPP
        ExprFunction(const nlohmann::json &j, const Tvarvec &vars={}) {
            Tconstvec consts;
            auto it = j.find("constants");
            if (it!=j.end())
                for (auto i=it->begin(); i!=it->end(); ++i)
                    consts.push_back( {i.key(), i.value()} );
            set(j.at("expr"), vars, consts);
        }
#endif
};

int main() {
    double x=0, y=0;
    nlohmann::json j = R"( { "expr":"x*x+kappa", "constants":{ "kappa":0.4, "f":2 } })"_json;
    std::function<double()> f = ExprFunction<double>( j, {{"x",&x}, {"y",&y}} );
    for (x=0; x<5; x+=1)
        std::cout << x << " " << f() << std::endl;
}
