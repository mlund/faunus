#include <doctest/doctest.h>
#include <stdexcept>
#include "atomdata.h"
#include "units.h"
#include "aux/eigensupport.h"
#include <spdlog/spdlog.h>

namespace Faunus {

bool InteractionData::contains(const key_type& name) const {
    if (auto it = data.find(name); it != data.end()) {
        return !std::isnan(it->second);
    }
    return false;
}

double InteractionData::at(const key_type& name) const {
    try {
        return data.at(name);
    } catch (const std::out_of_range& e) {
        // cannot throw until non-used atomic pairs are eliminated from the potential matrices
        faunus_logger->error("Unknown/unset atom property {} required.", name);
        return std::numeric_limits<double>::signaling_NaN();
    }
}

double& InteractionData::at(const key_type& name) {
    if (data.find(name) == data.end()) {
        insert_or_assign(name, std::numeric_limits<double>::signaling_NaN());
    }
    return data.at(name);
}

void InteractionData::insert_or_assign(const key_type& name, const double value) {
    auto it = data.find(name);
    if (it != data.end()) {
        it->second = value;
    } else {
        data.insert({name, value});
    }
}

void from_json(const json& j, InteractionData& a) {
    for (const auto& j_it : j.items()) {
        if (j_it.value().is_number()) {
            a.insert_or_assign(j_it.key(), j_it.value());
        }
    }
}

void from_single_use_json(SingleUseJSON& j, InteractionData& a) {
    auto j_copy = j;
    for (const auto& j_it : j_copy.items()) {
        if (j_it.value().is_number()) {
            a.insert_or_assign(j_it.key(), j_it.value());
            j.erase(j_it.key());
        }
    }
}

void to_json(json& j, const InteractionData& a) {
    for (const auto& kv : a.data) {
        j[kv.first] = kv.second;
    }
}

AtomData::index_type& AtomData::id() { return _id; }

const AtomData::index_type& AtomData::id() const { return _id; }

void to_json(json& j, const AtomData& a) {
    auto& _j = j[a.name];
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
          // sphero cylinders
          {"psc_length", a.psc_length / 1.0_angstrom},
          {"patch_type", a.patch_type},
          {"patch_angle", a.patch_angle / 1.0_deg},
          {"patch_angle_switch", a.patch_angle_switch / 1.0_deg},
          {"patch_attraction_range", a.patch_attraction_range / 1.0_angstrom},
          {"patch_cutoff", a.patch_cutoff / 1.0_angstrom},
          {"patch_chiral_angle", a.patch_chiral_angle / 1.0_deg},
          {"id", a.id()}};
    to_json(_j, a.interaction); // append other interactions
    if (a.hydrophobic)
        _j["hydrophobic"] = a.hydrophobic;
    if (a.implicit)
        _j["implicit"] = a.implicit;
}

void from_json(const json& j, AtomData& a) {
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

        // spherocylindrical stuff
        a.psc_length = val.value("psc_length", 0.0) * 1.0_angstrom;
        a.patch_type = val.value("patch_type", 0);
        a.patch_angle = val.value("patch_angle", 0.0) * 1.0_deg;
        a.patch_angle_switch = val.value("patch_angle_switch", 0.0) * 1.0_deg;
        a.patch_attraction_range = val.value("patch_attraction_range", 0.0) * 1.0_angstrom;
        a.patch_cutoff = val.value("patch_cutoff", 0.0) * 1.0_angstrom;
        a.patch_chiral_angle = val.value("patch_chiral_angle", 0.0) * 1.0_deg;

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
        a.interaction.insert_or_assign("sigma", val.value("sigma", 0.0) * 1.0_angstrom);
        if (fabs(a.interaction.at("sigma")) < 1e-20)
            a.interaction.insert_or_assign("sigma", 2.0 * val.value("r", 0.0) * 1.0_angstrom);
        // an ugly temporal hack needed until the refactorization is finished
        // as sigma is unfortunately accessed in loops
        a.sigma = a.interaction.at("sigma");

        from_single_use_json(val, a.interaction);
        if (!val.empty()) {
            usageTip.pick("atomlist");
            throw ConfigurationError("unused key(s) for atom '{}':\n{}", a.name, val.items().begin().key());
        }

        if (std::isnan(a.interaction.at("sigma")))
            throw ConfigurationError("no sigma parameter defined");
    }
}

void from_json(const json& j, std::vector<Faunus::AtomData>& atom_vector) {
    // List of atoms can be provided as an array or wrapped in an object containing an array – {"atomlist": […], …}.
    // Here we attempt to unwrap the object envelope.
    auto j_atomlist = j.find("atomlist");
    const auto& new_atoms = (j_atomlist == j.end()) ? j : *j_atomlist;
    if (!new_atoms.is_array()) {
        throw ConfigurationError("`atomlist` must be an array");
    }

    atom_vector.reserve(atom_vector.size() + new_atoms.size()); // reserve memory
    try {
        for (const auto& element : new_atoms) { // loop over elements in json array
            if (element.is_object()) {
                const auto atomdata = Faunus::AtomData(element);
                if (auto it = Faunus::findName(atom_vector, atomdata.name); it != atom_vector.end()) {
                    faunus_logger->warn("overwriting existing properties of {}", it->name);
                    *it = atomdata;
                } else { // add new atom
                    atom_vector.push_back(atomdata);
                }
            } else if (element.is_string()) { // treat element as filename
                const auto filename = element.get<std::string>();
                faunus_logger->info("reading atoms from external file '{}'", filename);
                from_json(Faunus::loadJSON(filename), atom_vector);
            } else {
                throw ConfigurationError("atom entry must be string or object").attachJson(element);
            }
        }
        assert(atom_vector.size() < std::numeric_limits<Faunus::AtomData::index_type>::max());
        // the id exactly matches it's position (index) in the atom vector
        for (size_t i = 0; i < atom_vector.size(); ++i) {
            atom_vector[i].id() = i;
        }
    } catch (std::exception& e) {
        std::throw_with_nested(ConfigurationError("error in atoms definition").attachJson(new_atoms));
    }
}

std::vector<Faunus::AtomData> atoms;

TEST_SUITE_BEGIN("AtomData");

TEST_CASE("[Faunus] AtomData") {
    using doctest::Approx;

    json j = R"({ "atomlist" : [
             { "A": { "sigma": 2.5, "pactivity":2, "eps_custom": 0.1 } },
             { "B": { "r":1.1, "activity":0.2, "eps":0.05, "dp":9.8, "dprot":3.14, "mw":1.1, "tfe":0.98, "tension":0.023 } }
             ]})"_json;

    pc::temperature = 298.15_K;
    atoms = j["atomlist"].get<decltype(atoms)>();
    auto& v = atoms; // alias to global atom list

    CHECK_EQ(v.size(), 2);
    CHECK_EQ(v.front().id(), 0);
    CHECK_EQ(v.front().name, "A");                                 // alphabetic order in std::map
    CHECK(v.front().interaction.at("sigma") == Approx(2.5));      // raw number, no units
    CHECK(v.front().interaction.at("eps_custom") == Approx(0.1)); // raw number, no units

    CHECK_EQ(std::isnan(v.front().interaction.at("eps_unknown")), true);
    // CHECK_THROWS_AS_MESSAGE(v.front().interaction.get("eps_unknown"), std::runtime_error, "unknown atom property");
    CHECK(v.front().sigma == Approx(2.5e-10_m));
    CHECK(v.front().activity == Approx(0.01_molar));
    CHECK(v.back().tfe == Approx(0.98_kJmol / (1.0_angstrom * 1.0_angstrom * 1.0_molar)));

    AtomData a = json(v.back()); // AtomData -> JSON -> AtomData

    CHECK_EQ(a.name, "B");
    CHECK_EQ(a.id(), 1);
    CHECK(a.activity == Approx(0.2_molar));
    CHECK(a.interaction.at("sigma") == Approx(2.2)); // raw number, no units
    CHECK(a.interaction.at("eps") == Approx(0.05));  // raw number, no units
    CHECK(a.dp == Approx(9.8));
    CHECK(a.dprot == Approx(3.14));
    CHECK(a.mw == Approx(1.1));
    CHECK(a.tfe == Approx(0.98_kJmol / 1.0_angstrom / 1.0_angstrom / 1.0_molar));
    CHECK(a.tension == Approx(0.023_kJmol / 1.0_angstrom / 1.0_angstrom));

    auto it = findName(v, "B");
    CHECK_EQ(it->id(), 1);
    it = findName(v, "unknown atom");
    CHECK_EQ(it, v.end());
}
TEST_SUITE_END();

UnknownAtomError::UnknownAtomError(std::string_view atom_name)
    : GenericError("unknown atom: '{}'", atom_name) {}

AtomData& findAtomByName(std::string_view name) {
    const auto result = findName(Faunus::atoms, name);
    if (result == Faunus::atoms.end()) {
        throw UnknownAtomError(name);
    }
    return *result;
}

} // namespace Faunus
