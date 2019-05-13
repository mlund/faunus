#include "functionparser.h"

template<typename T>
void ExprFunction<T>::set(const std::string &exprstr, const Tvarvec &vars, const Tconstvec &consts) {
    if (not parser)
        parser = std::make_shared<exprtk::parser<T>>();
    symbols.clear();
    for (auto &v : vars)
        symbols.add_variable(v.first, *v.second);
    for (auto &v : consts)
        symbols.add_constant(v.first, v.second);
    symbols.add_constants();
    expression.register_symbol_table(symbols);
    if (! parser->compile(exprstr, expression))
        throw std::runtime_error("error passing function/expression");
}

template<typename T>
void ExprFunction<T>::set(const nlohmann::json &j, const Tvarvec &vars) {
    Tconstvec consts;
    auto it = j.find("constants");
    if (it != j.end())
        for (auto i = it->begin(); i != it->end(); ++i)
            consts.push_back({i.key(), i.value()});
    set(j.at("function"), vars, consts);
}

template<typename T>
T ExprFunction<T>::operator()() const {
    return expression.value();
}

template class ExprFunction<double>;
