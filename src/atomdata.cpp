#include <spdlog/spdlog.h>
#include <cmath>
#include "atomdata.h"
#include "units.h"
#include "aux/eigensupport.h"

namespace Faunus {

bool InteractionData::has(const Tkey name) const {
    auto it = data.find(name);
    return it != data.end() && !std::isnan(it->second);
}

double InteractionData::get(const Tkey name) const {
    try {
        return data.at(name);
    } catch (const std::out_of_range &e) {
        // cannot throw until non-used atomic pairs are eliminated from the potential matrices
        faunus_logger->error("Unknown/unset atom property {} required.", name);
        return std::numeric_limits<double>::signaling_NaN();
    }
}

double& InteractionData::get(const Tkey name) {
    if(data.find(name) == data.end()) {
        set(name, std::numeric_limits<double>::signaling_NaN());
    }
    return data.at(name);
}

void InteractionData::set(const Tkey name, const double value) {
    auto it = data.find(name);
    if(it != data.end()) {
        it->second = value;
    } else {
        data.insert({name, value});
    }
}

void from_json(const json &j, InteractionData &a) {
    for (const auto &j_it : j.items()) {
        if (j_it.value().is_number()) {
            a.set(j_it.key(), j_it.value());
        }
    }
}

void from_single_use_json(SingleUseJSON &j, InteractionData &a) {
    auto j_copy = j;
    for (const auto &j_it : j_copy.items()) {
        if (j_it.value().is_number()) {
            a.set(j_it.key(), j_it.value());
            j.erase(j_it.key());
        }
    }
}

void to_json(json &j, const InteractionData &a) {
    for (const auto &kv : a.data) {
        j[kv.first] = kv.second;
    }
}

AtomData::Tid &AtomData::id() { return _id; }

const AtomData::Tid &AtomData::id() const { return _id; }

void to_json(json &j, const AtomData &a) {
    auto &_j = j[a.name];
    _j = {{"activity", a.activity / 1.0_molar},
          {"pactivity", -std::log10(a.activity / 1.0_molar)},
          {"alphax", a.alphax},
          {"mw", a.mw},
          {"q", a.charge},
          {"dp", a.dp / 1.0_angstrom},
          {"dprot", a.dprot / 1.0_rad},
          {"tension", a.tension * 1.0_angstrom * 1.0_angstrom / 1.0_kJmol},
          {"tfe", a.tfe * 1.0_angstrom * 1.0_angstrom * 1.0_molar / 1.0_kJmol},
          {"mu", a.mu},
          {"mulen", a.mulen},
          {"scdir", a.scdir},
          {"sclen", a.sclen},
          {"id", a.id()}};
    to_json(_j, a.interaction); // append other interactions
    if (a.hydrophobic)
        _j["hydrophobic"] = a.hydrophobic;
    if (a.implicit)
        _j["implicit"] = a.implicit;
}

void from_json(const json &j, AtomData &a) {
    if (j.is_object() == false || j.size() != 1)
        throw std::runtime_error("Invalid JSON data for AtomData");
    for (auto atom_iter : j.items()) {
        a.name = atom_iter.key();
        SingleUseJSON val = atom_iter.value();
        a.alphax = val.value("alphax", a.alphax);
        a.charge = val.value("q", a.charge);
        a.dp = val.value("dp", a.dp) * 1.0_angstrom;
        a.dprot = val.value("dprot", a.dprot) * 1.0_rad;
        a.id() = val.value("id", a.id());
        a.mu = val.value("mu", a.mu);
        a.mulen = val.value("mulen", a.mulen);
        if (a.mu.norm() > 1e-6) {          // if mu is given...
            if (std::fabs(a.mulen) < 1e-6) // if mulen is not given ...
                a.mulen = a.mu.norm();     // ... then set mulen
            a.mu = a.mu / a.mu.norm();     // normalize mu
        }
        a.scdir = val.value("scdir", a.scdir);
        a.sclen = val.value("sclen", a.sclen);
        a.mw = val.value("mw", a.mw);
        a.tension = val.value("tension", a.tension) * 1.0_kJmol / (1.0_angstrom * 1.0_angstrom);
        a.tfe = val.value("tfe", a.tfe) * 1.0_kJmol / (1.0_angstrom * 1.0_angstrom * 1.0_molar);
        a.hydrophobic = val.value("hydrophobic", false);
        a.implicit = val.value("implicit", false);
        if (val.count("activity") == 1)
            a.activity = val.at("activity").get<double>() * 1.0_molar;
        if (val.count("pactivity") == 1) {
            if (val.count("activity") == 1) {
                throw std::runtime_error("Specify either activity or pactivity for atom '"s + a.name + "'!");
            }
            a.activity = std::pow(10, -val.at("pactivity").get<double>()) * 1.0_molar;
        }

        if (val.count("r") == 1) {
            faunus_logger->warn("Atom property `r` is obsolete; use `sigma = 2*r` instead on atom {}", a.name);
        }
        a.interaction.set("sigma", val.value("sigma", 0.0) * 1.0_angstrom);
        if (fabs(a.interaction.get("sigma")) < 1e-20)
            a.interaction.set("sigma", 2.0 * val.value("r", 0.0) * 1.0_angstrom);
        // an ugly temporal hack needed until the refactorization is finished
        // as sigma is unfortunately accessed in loops
        a.sigma = a.interaction.get("sigma");

        from_single_use_json(val, a.interaction);
        if (!val.empty()) {
            throw std::runtime_error("unused key(s) for atom '"s + a.name + "':\n" + val.items().begin().key() +
                                     usageTip["atomlist"]);
        }

        if (std::isnan(a.interaction.get("sigma")))
            throw ConfigurationError("no sigma parameter defined");
    }
}

void from_json(const json &j, std::vector<AtomData> &v) {
    auto j_list_iter = j.find("atomlist");
    auto j_list = (j_list_iter == j.end()) ? j : *j_list_iter;

    v.reserve(v.size() + j_list.size());
    for (auto &j_atom : j_list) {
        if (j_atom.is_string()) // treat ax external file to load
            from_json(openjson(j_atom.get<std::string>()), v);
        else if (j_atom.is_object()) {
            AtomData a = j_atom;
            auto atom_iter = findName(v, a.name);
            if (atom_iter == v.end()) {
                v.push_back(a); // add new atom
            } else {
                faunus_logger->warn("Redefining atomic properties of atom {}.", atom_iter->name);
                *atom_iter = a;
            }
        }
    }
    for (size_t i = 0; i < v.size(); i++) {
        if (std::numeric_limits<AtomData::Tid>::max() < i) {
            throw std::overflow_error("Number of atoms to high.");
        }
        v[i].id() = i; // id must match position in vector
    }
}

std::vector<AtomData> atoms;

} // namespace Faunus
