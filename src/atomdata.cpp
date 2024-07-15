#pragma clang diagnostic push
#pragma ide diagnostic ignored "altera-unroll-loops"
#include <doctest/doctest.h>
#include <stdexcept>
#include "atomdata.h"
#include "units.h"
#include "aux/eigensupport.h"
#include <spdlog/spdlog.h>

namespace Faunus {

bool InteractionData::contains(const key_type& name) const
{
    if (auto it = data.find(name); it != data.end()) {
        return !std::isnan(it->second);
    }
    return false;
}

double InteractionData::at(const key_type& name) const
{
    try {
        return data.at(name);
    }
    catch (const std::out_of_range& e) {
        // cannot throw until non-used atomic pairs are eliminated from the potential matrices
        faunus_logger->error("Unknown/unset atom property {} required.", name);
        return std::numeric_limits<double>::signaling_NaN();
    }
}

double& InteractionData::at(const key_type& name)
{
    if (data.find(name) == data.end()) {
        insert_or_assign(name, std::numeric_limits<double>::signaling_NaN());
    }
    return data.at(name);
}

void InteractionData::insert_or_assign(const key_type& name, const double value)
{
    auto it = data.find(name);
    if (it != data.end()) {
        it->second = value;
    }
    else {
        data.insert({name, value});
    }
}

void from_json(const json& j, InteractionData& a)
{
    for (const auto& [key, value] : j.items()) {
        if (value.is_number()) {
            a.insert_or_assign(key, value);
        }
    }
}

void from_single_use_json(SingleUseJSON& j, InteractionData& a)
{
    auto j_copy = j;
    for (const auto& [key, value] : j_copy.items()) {
        if (value.is_number()) {
            a.insert_or_assign(key, value);
            j.erase(key);
        }
    }
}

void to_json(json& j, const InteractionData& a)
{
    for (const auto& [key, value] : a.data) {
        j[key] = value;
    }
}

AtomData::index_type& AtomData::id()
{
    return _id;
}

const AtomData::index_type& AtomData::id() const
{
    return _id;
}

void to_json(json& j, const AtomData& a)
{
    auto& _j = j[a.name];
    _j = {{"activity", a.activity / 1.0_molar},
          {"pactivity", -std::log10(a.activity / 1.0_molar)},
          {"alphax", a.alphax},
          {"mw", a.mw},
          {"q", a.charge},
          {"tension", a.tension * 1.0_angstrom * 1.0_angstrom / 1.0_kJmol},
          {"tfe", a.tfe * 1.0_angstrom * 1.0_angstrom * 1.0_molar / 1.0_kJmol},
          {"mu", a.mu},
          {"mulen", a.mulen},
          {"psc", a.sphero_cylinder},
          {"id", a.id()}};
    if (a.dp.has_value()) {
        _j["dp"] = a.dp.value() / 1.0_angstrom;
    }
    if (a.dprot.has_value()) {
        _j["dprot"] = a.dprot.value() / 1.0_rad;
    }
    to_json(_j, a.interaction); // append other interactions
    if (a.hydrophobic) {
        _j["hydrophobic"] = a.hydrophobic;
    }
    if (a.implicit) {
        _j["implicit"] = a.implicit;
    }
}

/// Handles optional translational and rotational displacement
void set_dp_and_dprot(const json& j, AtomData& atomdata)
{
    // todo: use `std::optional::and_then()` when C++23 is available
    if (const auto dp = get_optional<double>(j, "dp")) {
        atomdata.dp = dp.value() * 1.0_angstrom;
    }
    if (const auto dprot = get_optional<double>(j, "dprot")) {
        atomdata.dprot = dprot.value() * 1.0_rad;
        if (std::fabs(atomdata.dprot.value()) > 2.0 * pc::pi) {
            faunus_logger->warn("rotational displacement should be between [0:2π]");
        }
    }
}

void from_json(const json& j, AtomData& a)
{
    if (!j.is_object() || j.size() != 1) {
        throw std::runtime_error("Invalid JSON data for AtomData");
    }
    for (const auto& atom_iter : j.items()) {
        a.name = atom_iter.key();
        auto val = SingleUseJSON(atom_iter.value());
        a.alphax = val.value("alphax", a.alphax);
        a.charge = val.value("q", a.charge);
        a.id() = val.value("id", a.id());
        a.mu = val.value("mu", a.mu);
        a.mulen = val.value("mulen", a.mulen);
        if (a.mu.norm() > 1e-6) {            // if mu is given...
            if (std::fabs(a.mulen) < 1e-6) { // if mulen is not given ...
                a.mulen = a.mu.norm();       // ... then set mulen
            }
            a.mu = a.mu / a.mu.norm(); // normalize mu
        }
        a.sphero_cylinder = val.value("psc", SpheroCylinderData());
        a.mw = val.value("mw", a.mw);
        a.tension = val.value("tension", a.tension) * 1.0_kJmol / (1.0_angstrom * 1.0_angstrom);
        a.tfe = val.value("tfe", a.tfe) * 1.0_kJmol / (1.0_angstrom * 1.0_angstrom * 1.0_molar);
        a.hydrophobic = val.value("hydrophobic", false);
        a.implicit = val.value("implicit", false);
        set_dp_and_dprot(val, a);
        if (val.contains("activity")) {
            a.activity = val.at("activity").get<double>() * 1.0_molar;
        }
        if (val.contains("pactivity")) {
            if (val.contains("activity")) {
                throw std::runtime_error("Specify either activity or pactivity for atom '"s +
                                         a.name + "'!");
            }
            a.activity = std::pow(10, -val.at("pactivity").get<double>()) * 1.0_molar;
        }

        if (val.contains("r")) {
            faunus_logger->warn(
                "Atom property `r` is obsolete; use `sigma = 2*r` instead on atom {}", a.name);
        }
        a.interaction.insert_or_assign("sigma", val.value("sigma", 0.0) * 1.0_angstrom);
        if (fabs(a.interaction.at("sigma")) < 1e-20) {
            a.interaction.insert_or_assign("sigma", 2.0 * val.value("r", 0.0) * 1.0_angstrom);
        }
        // an ugly temporal hack needed until the refactorization is finished
        // as sigma is unfortunately accessed in loops
        a.sigma = a.interaction.at("sigma");

        from_single_use_json(val, a.interaction);
        if (!val.empty()) {
            usageTip.pick("atomlist");
            throw ConfigurationError("unused key(s) for atom '{}':\n{}", a.name,
                                     val.items().begin().key());
        }

        if (std::isnan(a.interaction.at("sigma"))) {
            throw ConfigurationError("no sigma parameter defined");
        }
    }
}

void from_json(const json& j, std::vector<Faunus::AtomData>& atom_vector)
{
    // List of atoms can be provided as an array or wrapped in an object containing an array –
    // {"atomlist": […], …}. Here we attempt to unwrap the object envelope.
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
                if (auto it = Faunus::findName(atom_vector, atomdata.name);
                    it != atom_vector.end()) {
                    faunus_logger->warn("overwriting existing properties of {}", it->name);
                    *it = atomdata;
                }
                else { // add new atom
                    atom_vector.push_back(atomdata);
                }
            }
            else if (element.is_string()) { // treat element as filename
                const auto filename = element.get<std::string>();
                faunus_logger->info("reading atoms from external file '{}'", filename);
                from_json(Faunus::loadJSON(filename), atom_vector);
            }
            else {
                throw ConfigurationError("atom entry must be string or object").attachJson(element);
            }
        }
        assert(atom_vector.size() < std::numeric_limits<Faunus::AtomData::index_type>::max());
        // the id exactly matches it's position (index) in the atom vector
        for (size_t i = 0; i < atom_vector.size(); ++i) {
            atom_vector[i].id() = i;
        }
    }
    catch (std::exception& e) {
        std::throw_with_nested(
            ConfigurationError("error in atoms definition").attachJson(new_atoms));
    }
}

std::vector<Faunus::AtomData> atoms;

TEST_SUITE_BEGIN("AtomData");

TEST_CASE("[Faunus] AtomData")
{
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
    CHECK_EQ(v.front().interaction.at("sigma"), Approx(2.5));      // raw number, no units
    CHECK_EQ(v.front().interaction.at("eps_custom"), Approx(0.1)); // raw number, no units

    CHECK_EQ(std::isnan(v.front().interaction.at("eps_unknown")), true);
    // CHECK_THROWS_AS_MESSAGE(v.front().interaction.get("eps_unknown"), std::runtime_error,
    // "unknown atom property");
    CHECK_EQ(v.front().sigma, Approx(2.5e-10_m));
    CHECK_EQ(v.front().activity, Approx(0.01_molar));
    CHECK_EQ(v.back().tfe, Approx(0.98_kJmol / (1.0_angstrom * 1.0_angstrom * 1.0_molar)));

    AtomData a = json(v.back()); // AtomData -> JSON -> AtomData

    CHECK_EQ(a.name, "B");
    CHECK_EQ(a.id(), 1);
    CHECK_EQ(a.activity, Approx(0.2_molar));
    CHECK_EQ(a.interaction.at("sigma"), Approx(2.2)); // raw number, no units
    CHECK_EQ(a.interaction.at("eps"), Approx(0.05));  // raw number, no units
    CHECK_EQ(a.dp, Approx(9.8));
    CHECK_EQ(a.dprot, Approx(3.14));
    CHECK_EQ(a.mw, Approx(1.1));
    CHECK_EQ(a.tfe, Approx(0.98_kJmol / 1.0_angstrom / 1.0_angstrom / 1.0_molar));
    CHECK_EQ(a.tension, Approx(0.023_kJmol / 1.0_angstrom / 1.0_angstrom));

    auto it = findName(v, "B");
    CHECK_EQ(it->id(), 1);
    it = findName(v, "unknown atom");
    CHECK_EQ(it, v.end());
}

TEST_SUITE_END();

UnknownAtomError::UnknownAtomError(std::string_view atom_name)
    : GenericError("unknown atom: '{}'", atom_name)
{
}

AtomData& findAtomByName(std::string_view name)
{
    const auto iter = findName(Faunus::atoms, name);
    if (iter == Faunus::atoms.end()) {
        throw UnknownAtomError(name);
    }
    return *iter;
}

void from_json(const json& j, SpheroCylinderData& psc)
{
    psc.length = j.value("length", 0.0) * 1.0_angstrom;
    psc.type = j.value("type", SpheroCylinderData::PatchType::None);
    psc.patch_angle = j.value("patch_angle", 0.0) * 1.0_deg;
    psc.patch_angle_switch = j.value("patch_angle_switch", 0.0) * 1.0_deg;
    psc.chiral_angle = j.value("chiral_angle", 0.0) * 1.0_deg;
}

void to_json(json& j, const SpheroCylinderData& psc)
{
    j = {{"length", psc.length / 1.0_angstrom},
         {"patch_angle", psc.patch_angle / 1.0_deg},
         {"patch_angle_switch", psc.patch_angle_switch / 1.0_deg},
         {"type", psc.type},
         {"chiral_angle", psc.chiral_angle / 1.0_deg}};
}

} // namespace Faunus

#pragma clang diagnostic pop