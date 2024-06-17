#include "potentials.h"
#include "multipole.h"
#include "units.h"
#include "auxiliary.h"
#include "aux/arange.h"
#include "smart_montecarlo.h"
#include <coulombgalore.h>
#include <spdlog/spdlog.h>

#include <utility>

namespace Faunus::pairpotential {

// =============== PairMixer ===============

TCombinatorFunc PairMixer::getCombinator(CombinationRuleType combination_rule, CoefficientType coefficient) {
    TCombinatorFunc combinator;
    switch (combination_rule) {
    case CombinationRuleType::UNDEFINED:
        combinator = &combUndefined;
        break;
    case CombinationRuleType::ARITHMETIC:
        combinator = &combArithmetic;
        break;
    case CombinationRuleType::GEOMETRIC:
        combinator = &combGeometric;
        break;
    case CombinationRuleType::LORENTZ_BERTHELOT:
        switch (coefficient) {
        case CoefficientType::SIGMA:
            combinator = &combArithmetic;
            break;
        case CoefficientType::EPSILON:
            combinator = &combGeometric;
            break;
        default:
            throw std::logic_error("unsupported mixer initialization");
        }
        break;
    default:
        throw std::logic_error("unsupported mixer initialization");
    }
    return combinator;
}

TPairMatrixPtr PairMixer::createPairMatrix(const std::vector<AtomData>& atoms) {
    size_t n = atoms.size(); // number of atom types
    TPairMatrixPtr matrix = std::make_shared<TPairMatrix>(n, n);
    for (const auto& i : atoms) {
        for (const auto& j : atoms) {
            if (i.implicit || j.implicit) {
                // implicit atoms are ignored as the missing properties, e.g., sigma and epsilon, might raise errors
                (*matrix)(i.id(), j.id()) = combUndefined();
            } else if (i.id() == j.id()) {
                // if the combinator is "undefined" the homogeneous interaction is still well-defined
                (*matrix)(i.id(), j.id()) = modifier(extractor(i.interaction));
            } else {
                (*matrix)(i.id(), j.id()) = modifier(combinator(extractor(i.interaction), extractor(j.interaction)));
            }
        }
    }
    return matrix;
}

TPairMatrixPtr PairMixer::createPairMatrix(const std::vector<AtomData> &atoms,
                                           const std::vector<CustomInteractionData> &interactions) {
    TPairMatrixPtr matrix = PairMixer::createPairMatrix(atoms);
    auto dimension = std::min(matrix->rows(), matrix->cols());
    for (const auto &i : interactions) {
        if (i.atom_id[0] < dimension && i.atom_id[1] < dimension) {
            // interaction is always symmetric
            (*matrix)(i.atom_id[0], i.atom_id[1]) = (*matrix)(i.atom_id[1], i.atom_id[0]) =
                modifier(extractor(i.interaction));
        } else {
            throw std::range_error("atomtype index out of range");
        }
    }
    return matrix;
}
PairMixer::PairMixer(TExtractorFunc extractor, TCombinatorFunc combinator, TModifierFunc modifier)
    : extractor(std::move(extractor)), combinator(std::move(combinator)), modifier(std::move(modifier)){}

TEST_CASE("[Faunus] PairMixer") {
    using namespace std::string_literals;
    using doctest::Approx;
    SUBCASE("Enumerated potential") {
        REQUIRE((PairMixer::combArithmetic(2.0, 8.0) == Approx(5.0)));
        REQUIRE((PairMixer::combGeometric(2.0, 8.0) == Approx(4.0)));
        CHECK_EQ(PairMixer::getCombinator(CombinationRuleType::LORENTZ_BERTHELOT, PairMixer::CoefficientType::SIGMA)(
                  2.0, 8.0), PairMixer::combArithmetic(2.0, 8.0));
        CHECK_EQ(PairMixer::getCombinator(CombinationRuleType::LORENTZ_BERTHELOT, PairMixer::CoefficientType::EPSILON)(
                  2.0, 8.0), PairMixer::combGeometric(2.0, 8.0));
        CHECK_THROWS_AS(PairMixer::getCombinator(CombinationRuleType::LORENTZ_BERTHELOT), std::logic_error);

        SUBCASE("") {
            atoms =
                R"([{"A": {"sigma":2.0}}, {"B": {"sigma":8.0}}, {"C": {"sigma":18.0}}])"_json.get<decltype(atoms)>();
            REQUIRE((atoms.front().interaction.at("sigma") == Approx(2.0)));
            std::vector<CustomInteractionData> pairs = R"([{"A C": {"sigma": 9.5}}, {"C B": {"sigma": 12.5}}])"_json;
            TExtractorFunc sigma = [](InteractionData a) -> double { return a.at("sigma"); };

            SUBCASE("") {
                PairMixer mixer(sigma, &PairMixer::combArithmetic);
                SUBCASE("Atom pairs") {
                    auto matrix = mixer.createPairMatrix(atoms);
                    CHECK(matrix->isApprox(matrix->transpose())); // symmetric
                    CHECK_EQ((*matrix)(0, 0), Approx(2.0));
                    CHECK_EQ((*matrix)(0, 1), Approx(5.0));
                }
                SUBCASE("Custom pairs") {
                    auto matrix = mixer.createPairMatrix(atoms, pairs);
                    CHECK(matrix->isApprox(matrix->transpose())); // symmetric
                    CHECK_EQ((*matrix)(0, 0), Approx(2.0));
                    CHECK_EQ((*matrix)(0, 1), Approx(5.0));
                    CHECK_EQ((*matrix)(2, 0), Approx(9.5));
                    CHECK_EQ((*matrix)(2, 1), Approx(12.5));
                }
            }
            SUBCASE("Modifier") {
                PairMixer mixer(sigma, &PairMixer::combArithmetic, [](double x) { return 10 * x; });
                auto matrix = mixer.createPairMatrix(atoms, pairs);
                CHECK(matrix->isApprox(matrix->transpose())); // symmetric
                CHECK_EQ((*matrix)(0, 0), Approx(20.0));
                CHECK_EQ((*matrix)(0, 1), Approx(50.0));
                CHECK_EQ((*matrix)(2, 0), Approx(95.0));
                CHECK_EQ((*matrix)(2, 1), Approx(125.0));
            }
            SUBCASE("Alternative JSON") {
                CHECK_NOTHROW(R"({"A C": {"sigma": 9.5}})"_json.get<std::vector<CustomInteractionData>>());
                std::vector<CustomInteractionData> alt_pairs =
                    R"({"A C": {"sigma": 9.5}, "C B": {"sigma": 12.5}})"_json;
                CHECK_EQ(alt_pairs.size(), pairs.size());
            }
        }
    }
}

// =============== CustomInteractionData ===============

void from_json(const json &j, CustomInteractionData &c) {
    if (!j.is_object() || j.size() != 1) {
        throw ConfigurationError("invalid JSON for custom interaction parameters");
    }
    auto j_item = j.items().begin();
    try {
        const auto atom_names = splitConvert<std::string>(j_item.key());
        const auto atom_ids = names2ids(atoms, atom_names);
        if (atom_ids.size() != c.atom_id.size()) {
            faunus_logger->error("Custom interaction parameters require exactly {} space-separated atoms: {}.",
                                 c.atom_id.size(), joinToString(atom_names));
            throw ConfigurationError("wrong number of atoms in custom interaction parameters");
        }
        std::copy(atom_ids.begin(), atom_ids.end(), c.atom_id.begin());
    } catch (std::out_of_range &e) {
        faunus_logger->error("Unknown atom pair [{}] in custom interaction parameters.", j_item.key());
        throw ConfigurationError("unknown atom pair in custom interaction parameters");
    }
    c.interaction = j_item.value();
}

void to_json(json &j, const CustomInteractionData &interaction) {
    std::vector<std::string> atom_names;
    for(auto atom_id : interaction.atom_id) {
        atom_names.push_back(atoms[atom_id].name);
    }
    j = {{joinToString(atom_names), interaction.interaction}};
}

void from_json(const json& j, std::vector<CustomInteractionData>& interactions) {
    auto append = [&](const auto& key, const auto& value) {
        interactions.push_back(json{{key, value}});
        faunus_logger->debug("Custom interaction for particle pair {}: {}", key, value.dump(-1));
    };
    if (j.is_array()) {
        for (const auto& j_pair : j) {
            if (j_pair.size() != 1) {
                throw ConfigurationError("Custom interaction input error");
            }
            const auto& [key, value] = j_pair.items().begin();
            append(key, value);
        }
    } else if (j.is_object()) {
        for (const auto& [key, value] : j.items()) {
            append(key, value);
        }
    } else {
        throw ConfigurationError("invalid JSON for custom interaction parameters");
    }
}

// =============== PairPotentialBase ===============

PairPotential::PairPotential(std::string name, std::string cite, bool isotropic)
    : name(std::move(name))
    , cite(std::move(cite))
    , isotropic(isotropic) {}

/**
 * @brief Calculates force on particle a due to another particle, b
 * @param particle_a Particle a
 * @param particle_b Particle b
 * @param squared_distance Squared norm |𝐚-𝐛|²
 * @param b_towards_a Distance vector 𝐛 -> 𝐚 = 𝐚 - 𝐛
 * @return Force on particle a due to particle b
 */
Point PairPotential::force([[maybe_unused]] const Particle& a, [[maybe_unused]] const Particle& b,
                           [[maybe_unused]] double squared_distance, [[maybe_unused]] const Point& b_towards_a) const {
    throw(std::logic_error("Force computation not implemented for this setup!"));
}

void to_json(json& j, const PairPotential& base) { base.name.empty() ? base.to_json(j) : base.to_json(j[base.name]); }

void from_json(const json& j, PairPotential& base) {
    try {
        if (not base.name.empty()) {
            if (j.count(base.name) == 1) {
                base.from_json(j.at(base.name));
                return;
            }
        }
        base.from_json(j);
    } catch (std::exception &e) {
        usageTip.pick(base.name);
        throw ConfigurationError(e);
    }
}

// =============== MixerPairPotentialBase ===============

void MixerPairPotentialBase::init() {
    json j_combination_rule = combination_rule;
    faunus_logger->debug("Combination rule {} in effect for the {} potential.", j_combination_rule.get<std::string>(),
                         name);
    initPairMatrices();
}

void MixerPairPotentialBase::from_json(const json &j) {
    try {
        if (j.contains("mixing")) {
            const json& mixing = j.at("mixing");
            combination_rule = mixing.get<CombinationRuleType>();
            if (combination_rule == CombinationRuleType::UNDEFINED && mixing != "undefined") {
                // an ugly hack because the first pair in the json ↔ enum mapping is silently selected by default
                throw PairPotentialException("unknown combination rule " + mixing.get<std::string>());
            }
        }
        if (j.contains("custom")) {
            // *custom_pairs = j["custom"]; // does not work, perhaps as from_json is also a method (a namespace conflict)
            pairpotential::from_json(j["custom"], *custom_pairs);
        }
    } catch (const ConfigurationError &e) {
        faunus_logger->error(std::string(e.what()) + " in potential " + name);
        throw std::runtime_error("error reading potential " + name + " from json");
    }
    extractorsFromJson(j);
    init();
}

void MixerPairPotentialBase::to_json(json &j) const {
    j["mixing"] = combination_rule;
    j.update(json_extra_params);
    if (!custom_pairs->empty()) {
        j["custom"] = *custom_pairs;
    }
}
void MixerPairPotentialBase::extractorsFromJson(const json&) {}
MixerPairPotentialBase::MixerPairPotentialBase(const std::string& name, const std::string& cite,
                                               CombinationRuleType combination_rule, bool isotropic)
    : PairPotential(name, cite, isotropic)
    , combination_rule(combination_rule) {}

// =============== RepulsionR3 ===============

void RepulsionR3::from_json(const json &j) {
    f = j.value("prefactor", 1.0);
    e = j.value("lj-prefactor", 1.0);
    s = j.value("sigma", 1.0);
}

void RepulsionR3::to_json(json &j) const { j = {{"prefactor", f}, {"lj-prefactor", e}, {"sigma", s}}; }

// =============== CosAttract ===============

void CosAttract::to_json(json &j) const {
    j = {{"eps", eps / 1.0_kJmol}, {"rc", rc / 1.0_angstrom}, {"wc", wc / 1.0_angstrom}};
}

void CosAttract::from_json(const json &j) {
    eps = j.at("eps").get<double>() * 1.0_kJmol;
    rc = j.at("rc").get<double>() * 1.0_angstrom;
    wc = j.at("wc").get<double>() * 1.0_angstrom;
    rc2 = rc * rc;
    c = pc::pi / 2 / wc;
    rcwc2 = pow((rc + wc), 2);
}
double CosAttract::cutOffSquared() const {
    return rcwc2;
}
CosAttract::CosAttract(const std::string& name)
    : PairPotential(name) {}

TEST_CASE("[Faunus] CosAttract") {
    Particle a, b;
    Point r1 = {0.5 - 0.001, 0.0, 0.0};       // r < r_c (-epsilon)
    Point r2 = {2.1 + 0.5 + 0.001, 0.0, 0.0}; // r > r_c + w_c (zero)
    Point r3 = {1.0, 0.0, 0.0};               // r_c > r < r_c + w_c (switching region)
    a.id = 0;
    b.id = 1;

    SUBCASE("basic") {
        auto j = R"({ "atomlist" : [
                 { "A": { "eps": 1.0, "rc": 0.5, "wc": 2.1 } },
                 { "B": { "eps": 1.0, "rc": 0.5, "wc": 2.1 } }]})"_json;
        Faunus::atoms = j["atomlist"].get<decltype(atoms)>();
        CosAttract pairpot;
        pairpotential::from_json(R"({ "cos2": {"eps": 1.0, "rc": 0.5, "wc": 2.1}})"_json, pairpot);

        CHECK_EQ(pairpot(a, b, r1.squaredNorm(), r1), doctest::Approx(-0.4033930777));
        CHECK_EQ(pairpot(a, b, r2.squaredNorm(), r2), doctest::Approx(0));
        CHECK_EQ(pairpot(a, b, r3.squaredNorm(), r3), doctest::Approx(-0.3495505642));

        CHECK_EQ(pairpot.force(a, b, r1.squaredNorm(), r1).x(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r1.squaredNorm(), r1).y(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r1.squaredNorm(), r1).y(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r2.squaredNorm(), r2).x(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r2.squaredNorm(), r2).y(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r2.squaredNorm(), r2).z(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r3.squaredNorm(), r3).x(), doctest::Approx(-0.2052334967));
        CHECK_EQ(pairpot.force(a, b, r3.squaredNorm(), r3).y(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r3.squaredNorm(), r3).z(), doctest::Approx(0));
    }

    SUBCASE("mixed - symmetric") {
        auto j = R"({ "atomlist" : [
                 { "A": { "eps": 1.0, "rc": 0.5, "wc": 2.1 } },
                 { "B": { "eps": 1.0, "rc": 0.5, "wc": 2.1 } }]})"_json;
        Faunus::atoms = j["atomlist"].get<decltype(atoms)>();
        CosAttractMixed pairpot;
        pairpotential::from_json(R"({ "cos2mix": {"mixing": "LB"}})"_json, pairpot);
        CHECK_EQ(pairpot(a, b, r1.squaredNorm(), r1), doctest::Approx(-0.4033930777));
        CHECK_EQ(pairpot(a, b, r2.squaredNorm(), r2), doctest::Approx(0));
        CHECK_EQ(pairpot(a, b, r3.squaredNorm(), r3), doctest::Approx(-0.3495505642));
        CHECK_EQ(pairpot.cutOffSquared(a.id, b.id), doctest::Approx(std::pow(0.5 + 2.1, 2)));

        CHECK_EQ(pairpot.force(a, b, r1.squaredNorm(), r1).x(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r1.squaredNorm(), r1).y(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r1.squaredNorm(), r1).y(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r2.squaredNorm(), r2).x(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r2.squaredNorm(), r2).y(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r2.squaredNorm(), r2).z(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r3.squaredNorm(), r3).x(), doctest::Approx(-0.2052334967));
        CHECK_EQ(pairpot.force(a, b, r3.squaredNorm(), r3).y(), doctest::Approx(0));
        CHECK_EQ(pairpot.force(a, b, r3.squaredNorm(), r3).z(), doctest::Approx(0));
    }

    SUBCASE("mixed - asymmetric") {
        auto j = R"({ "atomlist" : [
                 { "A": { "eps": 1.0, "rc": 0.5, "wc": 2.1 } },
                 { "B": { "eps": 0.5, "rc": 0.6, "wc": 1.9 } }]})"_json;
        Faunus::atoms = j["atomlist"].get<decltype(atoms)>();
        CosAttractMixed pairpot;
        pairpotential::from_json(R"({ "cos2mix": {"mixing": "LB"}})"_json, pairpot);
        CHECK_EQ(pairpot(a, b, r1.squaredNorm(), r1), doctest::Approx(-0.2852419807));
        CHECK_EQ(pairpot(a, b, r2.squaredNorm(), r2), doctest::Approx(0));
        CHECK_EQ(pairpot(a, b, r3.squaredNorm(), r3), doctest::Approx(-0.2510708423));
        CHECK_EQ(pairpot(a, b, std::pow(0.55 + 2.0 + 0.001, 2), Point::Zero()), doctest::Approx(0));
        CHECK_EQ(pairpot.cutOffSquared(a.id, b.id), doctest::Approx(6.5025));
    }

}

// =============== Coulomb (old) ===============

void Coulomb::to_json(json &j) const {
    j["epsr"] = pc::relativeDielectricFromBjerrumLength(bjerrum_length);
    j["lB"] = bjerrum_length;
}

void Coulomb::from_json(const json &j) {
    if (j.size() == 1 && j.is_object()) {
        bjerrum_length = pc::bjerrumLength(j.at("epsr").get<double>());
    } else {
        throw ConfigurationError("Plain Coulomb potential expects 'epsr' key (only)");
    }
}
Coulomb::Coulomb(const std::string& name)
    : PairPotential(name) {}

// =============== DipoleDipole (old) ===============

void DipoleDipole::to_json(json &j) const {
    j["epsr"] = pc::relativeDielectricFromBjerrumLength(bjerrum_length);
    j["lB"] = bjerrum_length;
}

void DipoleDipole::from_json(const json& j) { bjerrum_length = pc::bjerrumLength(j.at("epsr")); }

DipoleDipole::DipoleDipole(const std::string& name, const std::string& cite)
    : PairPotential(name, cite, false) {}

// =============== FENE ===============

void FENE::from_json(const json &j) {
    k = j.at("stiffness");
    r02 = std::pow(double(j.at("maxsep")), 2);
    r02inv = 1 / r02;
}

void FENE::to_json(json &j) const { j = {{"stiffness", k}, {"maxsep", std::sqrt(r02)}}; }
FENE::FENE(const std::string& name)
    : PairPotential(name) {}

// =============== SASApotential ===============

/**
 *
 * @param R Radius of first sphere
 * @param r Readius of second sphere
 * @param center_center_distance_squared Distance between spheres centers
 * @return Total surface area. If `shift=true` then the areas of the two individual
 *         spheres are subtracted, meaning that the returned surface area is always
 *         smaller than or equal to zero. I.e. the surface area is relative to the
 *         maximum surface area (fully separated).
 *
 * Used formulas:
 * - https://mathworld.wolfram.com/Sphere-SphereIntersection.html
 * - https://mathworld.wolfram.com/SphericalCap.html
 */
double SASApotential::area(double R, double r, double center_center_distance_squared) const {
    R += proberadius;
    r += proberadius;
    const auto spheres_area = 4.0 * pc::pi * (R * R + r * r); // full surface area of both spheres
    const auto offset = shift ? spheres_area : 0.0;
    if (center_center_distance_squared > (R + r) * (R + r)) {
        return spheres_area - offset; // spheres do not overlap
    }
    if (r > R) { // always make R the bigger sphere
        std::swap(r, R);
    }
    const auto d = sqrt(center_center_distance_squared);
    if (d + r <= R) {
        return 4.0 * pc::pi * R * R - offset; // full volume of biggest sphere
    }
    auto h1 = (r - R + d) * (r + R - d) / (2.0 * d);         // height of first spherical cap
    auto h2 = (R - r + d) * (R + r - d) / (2.0 * d);         // height of second spherical cap
    const auto lens_area = 2.0 * pc::pi * (R * h1 + r * h2); // area of lens defined by caps
    return spheres_area - lens_area - offset;
}

void SASApotential::from_json(const json &j) {
    shift = j.value("shift", true);
    conc = j.at("molarity").get<double>() * 1.0_molar;
    proberadius = j.value("radius", 1.4) * 1.0_angstrom;
}

void SASApotential::to_json(json &j) const {
    j["molarity"] = conc / 1.0_molar;
    j["radius"] = proberadius / 1.0_angstrom;
    j["shift"] = shift;
}
SASApotential::SASApotential(const std::string& name, const std::string& cite)
    : PairPotential(name, cite) {}

TEST_CASE("[Faunus] SASApotential") {
    using doctest::Approx;
    json j = R"({ "atomlist" : [
                 { "A": { "r": 1.5, "tension": 0.023} },
                 { "B": { "r": 2.1, "tfe": 0.98 } }]})"_json;

    atoms = j["atomlist"].get<decltype(atoms)>();
    Particle a, b;
    a.id = 0;
    b.id = 1;
    SASApotential pot;
    json in = R"({ "sasa": {"molarity": 1.0, "radius": 0.0, "shift":false}})"_json;
    pairpotential::from_json(in["sasa"], pot);
    double conc = 1.0 * 1.0_molar;
    double tension = atoms.at(a.id).tension / 2;
    double tfe = atoms[b.id].tfe / 2;
    double f = tension + conc * tfe;
    CHECK((tension > 0.0));
    CHECK((conc > 0.0));
    CHECK((tfe > 0.0));
    CHECK((f > 0.0));
    CHECK_EQ(in, json(pot));
    CHECK_EQ(pot(a, b, 0, {0, 0, 0}), Approx(f * 4 * pc::pi * 2.1 * 2.1));                      // complete overlap
    CHECK_EQ(pot(a, b, 10 * 10, {10, 0, 0}), Approx(f * 4 * pc::pi * (2.1 * 2.1 + 1.5 * 1.5))); // far apart
    CHECK_EQ(pot(a, b, 2.5 * 2.5, {2.5, 0, 0}), Approx(f * 71.74894965974514));                 // partial overlap
}

// =============== CustomPairPotential ===============

void CustomPairPotential::from_json(const json& j) {
    squared_cutoff_distance = std::pow(j.value("cutoff", pc::infty), 2);
    original_input = j;
    auto& constants = original_input["constants"];
    if (constants == nullptr) {
        constants = json::object();
    }
    constants["e0"] = pc::vacuum_permittivity;
    constants["kB"] = pc::boltzmann_constant;
    constants["kT"] = pc::kT();
    constants["Nav"] = pc::avogadro;
    constants["Rc"] = std::sqrt(squared_cutoff_distance);
    constants["T"] = pc::temperature;
    expression.set(original_input, {{"r", &symbols->distance},
                                    {"charge1", &symbols->charge1},
                                    {"charge2", &symbols->charge2},
                                    {"s1", &symbols->sigma1},
                                    {"s2", &symbols->sigma2}});
}

void CustomPairPotential::to_json(json& j) const {
    j = original_input;
    if (std::isfinite(squared_cutoff_distance)) {
        j["cutoff"] = std::sqrt(squared_cutoff_distance);
    }
}
CustomPairPotential::CustomPairPotential(const std::string& name)
    : PairPotential(name)
    , symbols(std::make_shared<Symbols>()) {}

TEST_CASE("[Faunus] CustomPairPotential") {
    using doctest::Approx;
    json j = R"({ "atomlist" : [
                 {"A": { "q":1.0,  "r":3, "eps":0.1 }},
                 {"B": { "q":-1.0, "r":4, "eps":0.05 }} ]})"_json;
    Faunus::atoms = j["atomlist"].get<decltype(Faunus::atoms)>();
    Particle a, b;
    a = Faunus::atoms[0];
    b = Faunus::atoms[1];

    SUBCASE("energy") {
        CustomPairPotential pot;
        pairpotential::from_json(R"({"constants": { "kappa": 30, "lB": 7},
                              "function": "lB * charge1 * charge2 / (s1+s2) * exp(-kappa/r) * kT + pi"})"_json,
                                 pot);
        CHECK_EQ(pot(a, b, 2 * 2, {0, 0, 2}), Approx(-7.0 / (3.0 + 4.0) * std::exp(-30.0 / 2.0) * pc::kT() + pc::pi));
    }
    SUBCASE("force") {
        CustomPairPotential pot;
        pairpotential::from_json(
            R"({"constants": { "lB": 7.0056973292 }, "function": "lB * charge1 * charge2 / r"})"_json, pot);
        NewCoulombGalore coulomb;
        pairpotential::from_json(R"({ "coulomb": {"epsr": 80.0, "type": "plain"} } )"_json, coulomb);
        Point r = {coulomb.bjerrum_length, 0.2, -0.1};
        auto r2 = r.squaredNorm();
        auto force_ref = coulomb.force(a, b, r2, r);
        auto force = pot.force(a, b, r2, r);
        CHECK_EQ(coulomb.bjerrum_length, Approx(7.0056973292));
        CHECK_EQ(force.norm(), Approx(0.1425956964));
        CHECK_EQ(force.x(), Approx(force_ref.x()));
        CHECK_EQ(force.y(), Approx(force_ref.y()));
        CHECK_EQ(force.z(), Approx(force_ref.z()));
    }
}

// =============== Dummy ===============

// =============== LennardJones ===============

void LennardJones::initPairMatrices() {
    const TCombinatorFunc comb_sigma = PairMixer::getCombinator(combination_rule, PairMixer::CoefficientType::SIGMA);
    const TCombinatorFunc comb_epsilon =
        PairMixer::getCombinator(combination_rule, PairMixer::CoefficientType::EPSILON);

    sigma_squared = PairMixer(extract_sigma, comb_sigma, &PairMixer::modSquared).createPairMatrix(atoms, *custom_pairs);
    epsilon_quadruple = PairMixer(extract_epsilon, comb_epsilon, [](double x) -> double {
                            return 4 * x;
                        }).createPairMatrix(atoms, *custom_pairs);

    faunus_logger->trace("Pair matrices for {} sigma ({}×{}) and epsilon ({}×{}) created using {} custom pairs.", name,
                         sigma_squared->rows(), sigma_squared->cols(), epsilon_quadruple->rows(),
                         epsilon_quadruple->cols(), custom_pairs->size());
}

void LennardJones::extractorsFromJson(const json &j) {
    auto sigma_name = j.value("sigma", "sigma");
    json_extra_params["sigma"] = sigma_name;
    extract_sigma = [sigma_name](const InteractionData &a) -> double { return a.at(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const InteractionData &a) -> double { return a.at(epsilon_name) * 1.0_kJmol; };
}
LennardJones::LennardJones(const std::string& name, const std::string& cite, CombinationRuleType combination_rule)
    : MixerPairPotentialBase(name, cite, combination_rule) {}

TEST_CASE("[Faunus] LennardJones") {
    atoms = R"([{"A": {"sigma":2.0, "eps":0.9}},
                {"B": {"sigma":8.0, "eps":0.1}},
                {"C": {"sigma":5.0, "eps":1.1}}])"_json.get<decltype(atoms)>();
    Particle a = atoms[0];
    Particle b = atoms[1];

    const auto d = 0.9_nm;
    auto lj_func = [d](double sigma, double eps) -> double {
        return 4 * eps * (std::pow(sigma / d, 12) - std::pow(sigma / d, 6));
    };

    SUBCASE("JSON initilization") {
        CHECK_THROWS_AS(makePairPotential<LennardJones>(R"({"mixing": "unknown"})"_json), std::runtime_error);
        CHECK_NOTHROW(makePairPotential<LennardJones>(R"({})"_json));
        // alternative notation for custom as an array: custom: []
        CHECK_NOTHROW(
            makePairPotential<LennardJones>(R"({"mixing": "LB", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json));
    }
    SUBCASE("JSON output custom") {
        json j_in = R"({"mixing": "LB", "sigma": "sigma", "eps": "eps",
            "custom": [{"A B": {"eps": 0.5, "sigma": 8}}, {"A C": {"eps": 1.0, "sigma": 5}}]})"_json;
        auto lj = makePairPotential<LennardJones>(j_in);
        json j_out = lj;
        // The custom pair data in the output contain a lot of ballast among the original data.
        // If the original data match the output data, patching of the output shall not change it.
        json j_custom = j_out["lennardjones"]["custom"];
        j_custom[0].merge_patch(j_in["custom"][0]);
        j_custom[1].merge_patch(j_in["custom"][1]);
        CHECK_EQ(j_out["lennardjones"]["custom"], j_custom);
    }
    SUBCASE("Lorentz-Berthelot mixing") {
        using doctest::Approx;
        auto lj = makePairPotential<LennardJones>(R"({"mixing": "LB"})"_json);
        CHECK_EQ(lj(a, a, d * d, {0, 0, d}), Approx(lj_func(0.2_nm, 0.9_kJmol)));
        CHECK_EQ(lj(a, b, d * d, {0, 0, d}), Approx(lj_func(0.5_nm, 0.3_kJmol)));
    }
    SUBCASE("Geometric mixing") {
        using doctest::Approx;
        auto lj = makePairPotential<LennardJones>(R"({"mixing": "geometric"})"_json);
        CHECK_EQ(lj(a, a, d * d, {0, 0, d}), Approx(lj_func(0.2_nm, 0.9_kJmol)));
        CHECK_EQ(lj(a, b, d * d, {0, 0, d}), Approx(lj_func(0.4_nm, 0.3_kJmol)));
    }
    SUBCASE("Custom pairs") {
        using doctest::Approx;
        auto lj =
            makePairPotential<LennardJones>(R"({"mixing": "LB", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json);
        CHECK_EQ(lj(a, b, d * d, {0, 0, d}), Approx(lj_func(0.8_nm, 0.5_kJmol)));
        CHECK_EQ(lj(a, a, d * d, {0, 0, d}), Approx(lj_func(0.2_nm, 0.9_kJmol)));
    }

    SUBCASE("Force") {
        using doctest::Approx;
        auto lj =
            makePairPotential<LennardJones>(R"({"mixing": "LB", "custom": [{"A B": {"eps": 2.0, "sigma": 8}}]})"_json);
        a.pos = {0, 0, 0};
        b.pos = {9, 0, 0};
        Point b_towards_a = a.pos - b.pos;
        Point force = lj.force(a, b, b_towards_a.squaredNorm(), b_towards_a);
        CHECK_EQ(force.x(), Approx(0.0142838474)); // force on particle a
    }
}

// =============== HardSphere ===============

void HardSphere::initPairMatrices() {
    sigma_squared = PairMixer(extract_sigma, PairMixer::getCombinator(combination_rule), &PairMixer::modSquared)
                        .createPairMatrix(atoms, *custom_pairs);
    faunus_logger->trace("Pair matrix for {} sigma ({}×{}) created using {} custom pairs.", name, sigma_squared->rows(),
                         sigma_squared->cols(), custom_pairs->size());
}

void HardSphere::extractorsFromJson(const json &j) {
    auto sigma_name = j.value("sigma", "sigma");
    json_extra_params["sigma"] = sigma_name;
    extract_sigma = [sigma_name](const InteractionData &a) -> double { return a.at(sigma_name) * 1.0_angstrom; };
}
HardSphere::HardSphere(const std::string& name)
    : MixerPairPotentialBase(name, std::string(), CombinationRuleType::ARITHMETIC){}

TEST_CASE("[Faunus] HardSphere") {
    atoms = R"([{"A": {"sigma": 2}}, {"B": {"sigma": 8}}])"_json.get<decltype(atoms)>();
    Particle a = atoms[0], b = atoms[1];

    SUBCASE("JSON initialization") {
        CHECK_NOTHROW(makePairPotential<HardSphere>(R"({})"_json));
        CHECK_THROWS_AS(makePairPotential<HardSphere>(R"({"mixing": "unknown"})"_json), std::runtime_error);
    }
    SUBCASE("Undefined mixing") {
        auto hs = makePairPotential<HardSphere>(R"({"mixing": "undefined"})"_json);
        CHECK_EQ(hs(a, a, 1.99 * 1.99, {0, 0, 1.99}), pc::infty);
        CHECK_EQ(hs(a, a, 2.01 * 2.01, {0, 0, 2.01}), 0.0);
        // CHECK(std::isnan(hs(a, b, {0, 0, 4.99}))); // fails
        // CHECK(std::isnan(hs(a, b, {0, 0, 5.01}))); // fails
    }
    SUBCASE("Arithmetic mixing") {
        auto hs = makePairPotential<HardSphere>(R"({"mixing": "arithmetic"})"_json);
        CHECK_EQ(hs(a, a, 2.01_angstrom * 2.01_angstrom, {0, 0, 2.01_angstrom}), 0);
        CHECK_EQ(hs(a, a, 1.99_angstrom * 1.99_angstrom, {0, 0, 1.99_angstrom}), pc::infty);
        CHECK_EQ(hs(a, b, 5.01_angstrom * 5.01_angstrom, {0, 0, 5.01_angstrom}), 0);
        CHECK_EQ(hs(a, b, 4.99_angstrom * 4.99_angstrom, {0, 0, 4.99_angstrom}), pc::infty);
    }
    SUBCASE("Custom pairs with implicit mixing") {
        auto hs = makePairPotential<HardSphere>(R"({"custom": [{"A B": {"sigma": 6}}]})"_json);
        CHECK_EQ(hs(a, a, 2.01_angstrom * 2.01_angstrom, {0, 0, 2.01_angstrom}), 0);
        CHECK_EQ(hs(a, a, 1.99_angstrom * 1.99_angstrom, {0, 0, 1.99_angstrom}), pc::infty);
        CHECK_EQ(hs(a, b, 6.01_angstrom * 6.01_angstrom, {0, 0, 6.01_angstrom}), 0);
        CHECK_EQ(hs(a, b, 5.99_angstrom * 5.99_angstrom, {0, 0, 5.99_angstrom}), pc::infty);
    }
}

// =============== Hertz ===============

void Hertz::initPairMatrices() {
    const TCombinatorFunc comb_diameter = PairMixer::getCombinator(combination_rule, PairMixer::CoefficientType::SIGMA);
    const TCombinatorFunc comb_epsilon =
        PairMixer::getCombinator(combination_rule, PairMixer::CoefficientType::EPSILON);

    sigma_squared =
        PairMixer(extract_sigma, comb_diameter, &PairMixer::modSquared).createPairMatrix(atoms, *custom_pairs);
    epsilon = PairMixer(extract_epsilon, comb_epsilon).createPairMatrix(atoms, *custom_pairs);

    faunus_logger->trace("Pair matrix for {} radius ({}×{}) and epsilon ({}×{}) created using {} custom pairs.", name,
                         sigma_squared->rows(), sigma_squared->cols(), epsilon->rows(), epsilon->cols(),
                         custom_pairs->size());
}

void Hertz::extractorsFromJson(const json &j) {
    auto sigma_name = j.value("sigma", "sigma");
    json_extra_params["sigma"] = sigma_name;
    extract_sigma = [sigma_name](const InteractionData &a) -> double { return a.at(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const InteractionData &a) -> double { return a.at(epsilon_name) * 1.0_kJmol; };
}
Hertz::Hertz(const std::string& name)
    : MixerPairPotentialBase(name) {}

TEST_CASE("[Faunus] Hertz") {
    json j = R"({ "atomlist" : [
                 { "A": { "eps": 1.0, "sigma": 1.3} },
                 { "B": { "eps": 2.0, "sigma": 1.0 } }]})"_json;

    atoms = j["atomlist"].get<decltype(atoms)>();
    Particle a = atoms[0], b = atoms[1];

    SUBCASE("JSON serialization") {
        json hertz_json = R"({ "hertz": {"mixing": "lorentz_berthelot", "eps": "eps", "sigma": "sigma"}})"_json;
        auto hertz = makePairPotential<Hertz>(hertz_json);
        CHECK_EQ(hertz_json, json(hertz));
    }
    SUBCASE("Lorentz-Berthelot mixing") {
        using doctest::Approx;
        pc::temperature = 298.15_K;
        auto hertz = makePairPotential<Hertz>(R"({"mixing": "lorentz_berthelot"})"_json);
        CHECK_EQ(hertz(a, b, 0.7 * 0.7, {0.7, 0, 0}), Approx(0.0546424449)); // within cut-off
        CHECK_EQ(hertz(a, b, 1.15 * 1.15, {1.15, 0, 0}), Approx(0.0));       // at cut-off
        CHECK_EQ(hertz(a, b, 2.0 * 2.0, {2.0, 0, 0}), Approx(0.0));          // outside of cut-off
    }
}

// =============== SquareWell ===============

void SquareWell::initPairMatrices() {
    const TCombinatorFunc comb_diameter = PairMixer::getCombinator(combination_rule, PairMixer::CoefficientType::SIGMA);
    const TCombinatorFunc comb_depth = PairMixer::getCombinator(combination_rule, PairMixer::CoefficientType::EPSILON);

    sigma_squared =
        PairMixer(extract_sigma, comb_diameter, &PairMixer::modSquared).createPairMatrix(atoms, *custom_pairs);
    epsilon = PairMixer(extract_epsilon, comb_depth).createPairMatrix(atoms, *custom_pairs);

    faunus_logger->trace("Pair matrix for {} diameter ({}×{}) and depth ({}×{}) created using {} custom pairs.", name,
                         sigma_squared->rows(), sigma_squared->cols(), epsilon->rows(), epsilon->cols(),
                         custom_pairs->size());
}

void SquareWell::extractorsFromJson(const json &j) {
    auto sigma_name = j.value("sigma", "sigma");
    json_extra_params["sigma"] = sigma_name;
    extract_sigma = [sigma_name](const InteractionData &a) -> double { return a.at(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const InteractionData &a) -> double { return a.at(epsilon_name) * 1.0_kJmol; };
}
SquareWell::SquareWell(const std::string& name)
    : MixerPairPotentialBase(name) {}

TEST_CASE("[Faunus] SquareWell") {
    atoms = R"([{"A": { "r":5,  "sigma":4, "eps":0.2 }},
                {"B": { "r":10, "sigma":2, "eps":0.1 }} ])"_json.get<decltype(atoms)>();
    Particle a = atoms[0];
    Particle b = atoms[1];

    SUBCASE("JSON initilization") {
        CHECK_THROWS_AS(makePairPotential<SquareWell>(R"({"mixing": "unknown"})"_json), std::runtime_error);
        CHECK_NOTHROW(makePairPotential<SquareWell>(R"({})"_json));
    }
    SUBCASE("Undefined mixing") {
        using doctest::Approx;
        auto sw = makePairPotential<SquareWell>(R"({"mixing": "undefined"})"_json);
        CHECK_EQ(sw(a, a, 3.99 * 3.99, {0, 0, 0}), Approx(-0.2_kJmol));
        CHECK_EQ(sw(a, a, 4.01 * 4.01, {0, 0, 0}), Approx(0.0));
        // CHECK(std::isnan(sw(a, b, {0, 0, 5.99}))); // fails
        // CHECK(std::isnan(sw(a, b, {0, 0, 6.01}))); // fails
    }
    SUBCASE("Lorentz-Berthelot mixing") {
        using doctest::Approx;
        auto sw = makePairPotential<SquareWell>(R"({"mixing": "LB"})"_json);
        CHECK_EQ(sw(a, b, 2.99 * 2.99, {0, 0, 0}), Approx(-std::sqrt(0.2_kJmol * 0.1_kJmol)));
        CHECK_EQ(sw(a, b, 3.01 * 3.01, {0, 0, 0}), Approx(0));
    }
}

// =============== Polarizability ===============

void Polarizability::from_json(const json& j) {
    epsr = j.at("epsr").get<double>();
    const auto bjerrum_length = pc::bjerrumLength(epsr);
    for (const auto& i : Faunus::atoms) {
        for (const auto& j : Faunus::atoms) {
            m_neutral->set(i.id(), j.id(), -3 * i.alphax * pow(0.5 * i.sigma, 3) * j.alphax * pow(0.5 * j.sigma, 3));
            m_charged->set(i.id(), j.id(),
                           -bjerrum_length * 0.5 *
                               (pow(i.charge, 2) * j.alphax * pow(0.5 * j.sigma, 3) +
                                pow(j.charge, 2) * i.alphax * pow(0.5 * i.sigma, 3)));
        }
    }
}
Polarizability::Polarizability(const std::string& name)
    : Coulomb(name) {
    m_neutral = std::make_shared<PairMatrix<double>>();
    m_charged = std::make_shared<PairMatrix<double>>();
}
void Polarizability::to_json(json& j) const { j = {{"epsr", epsr}}; }

// =============== FunctorPotential ===============

void FunctorPotential::registerSelfEnergy(PairPotential* pot) {
    if (pot->selfEnergy) {
        if (not selfEnergy) { // no self energy is defined
            selfEnergy = pot->selfEnergy;
        } else // accumulate self energies
        {
            selfEnergy = [pot = pot, &selfEnergy = selfEnergy](const Particle& p) {
                return pot->selfEnergy(p) + selfEnergy(p);
            };
        }
        faunus_logger->debug("Added selfEnergy function from {} to {}", pot->name, name);
    } else {
        faunus_logger->trace("Failed to register non-defined selfEnergy() for {}", pot->name);
    }
}

FunctorPotential::EnergyFunctor FunctorPotential::combinePairPotentials(json& potential_array) {
    if (!potential_array.is_array()) {
        throw std::runtime_error("potential array required");
    }
    EnergyFunctor func = [](auto&, auto&, auto, auto&) { return 0.0; };
    for (auto& single_record : potential_array) { // loop over all defined potentials in array
        if (!single_record.is_object() || single_record.size() != 1) {
            continue;
        }
        for (auto& [name, j_config] : single_record.items()) {
            EnergyFunctor new_func = nullptr;
            try {
                if (name == "custom") {
                    new_func = makePairPotential<CustomPairPotential>(j_config);
                }
                // add Coulomb potential and self-energy
                // terms if not already added
                else if (name == "coulomb") { // temporary key
                    auto coulomb = makePairPotential<NewCoulombGalore>(j_config);
                    pairpotential::to_json(j_config, coulomb); // store output
                    if (!have_monopole_self_energy) {
                        registerSelfEnergy(&coulomb);
                        have_monopole_self_energy = true;
                    }
                    new_func = coulomb;
                } else if (name == "cos2") {
                    new_func = makePairPotential<pairpotential::CosAttract>(single_record);
                } else if (name == "polar") {
                    new_func = makePairPotential<pairpotential::Polarizability>(single_record);
                } else if (name == "hardsphere") {
                    new_func = makePairPotential<pairpotential::HardSphere>(single_record);
                } else if (name == "lennardjones") {
                    new_func = makePairPotential<LennardJones>(single_record);
                } else if (name == "repulsionr3") {
                    new_func = makePairPotential<RepulsionR3>(single_record);
                } else if (name == "sasa") {
                    new_func = makePairPotential<SASApotential>(single_record);
                } else if (name == "wca") {
                    new_func = makePairPotential<WeeksChandlerAndersen>(single_record);
                } else if (name == "pm") {
                    new_func = makePairPotential<PrimitiveModel>(j_config);
                } else if (name == "pmwca") {
                    new_func = makePairPotential<PrimitiveModelWCA>(j_config);
                } else if (name == "hertz") {
                    new_func = makePairPotential<Hertz>(single_record);
                } else if (name == "squarewell") {
                    new_func = makePairPotential<SquareWell>(single_record);
                } else if (name == "dipoledipole") {
                    faunus_logger->error("'{}' is deprecated, use 'multipole' instead", name);
                } else if (name == "stockmayer") {
                    faunus_logger->error("'{}' is deprecated, use 'lennardjones'+'multipole' instead", name);
                } else if (name == "multipole") {
                    auto multipole = makePairPotential<Multipole>(j_config);
                    pairpotential::to_json(j_config, multipole);
                    isotropic = false; // potential is now angular dependent
                    if (!have_dipole_self_energy) {
                        registerSelfEnergy(&multipole);
                        have_dipole_self_energy = true;
                    }
                    new_func = multipole;
                } else if (name == "hs-cigar") {
                    new_func = makePairPotential<HardSpheroCylinder>(j_config);
                } else if (name == "coswca-psc") {
                    new_func = makePairPotential<CigarCosAttractWCA>(j_config);
                }
                // place additional potentials here...
            } catch (std::exception& e) {
                usageTip.pick(name);
                throw ConfigurationError("{} -> {}", name, e.what());
            }
            if (new_func == nullptr) {
                throw ConfigurationError("potential '{}': unknown potential", name);
            }
            func = [func, new_func](const auto& a, const auto& b, auto r2, auto& r) {
                return func(a, b, r2, r) + new_func(a, b, r2, r);
            };
        }
    }
    return func;
}

void FunctorPotential::to_json(json &j) const {
    j["functor potential"] = backed_up_json_input;
    j["selfenergy"] = {{"monopole", have_monopole_self_energy}, {"dipole", have_dipole_self_energy}};
}

void FunctorPotential::from_json(const json &j) {
    have_monopole_self_energy = false;
    have_dipole_self_energy = false;
    backed_up_json_input = j;
    umatrix = decltype(umatrix)(atoms.size(), combinePairPotentials(backed_up_json_input.at("default")));
    for (const auto& [key, value] : backed_up_json_input.items()) {
        auto atompair = splitConvert<std::string>(key); // is this for a pair of atoms?
        if (atompair.size() == 2) {
            auto ids = names2ids(atoms, atompair);
            umatrix.set(ids[0], ids[1], combinePairPotentials(value));
        }
    }
}

FunctorPotential::FunctorPotential(const std::string& name)
    : PairPotential(name) {}

TEST_CASE("[Faunus] FunctorPotential") {
    using doctest::Approx;

    json j = R"({ "atomlist" : [
                 {"A": { "q":1.0,  "r":1.1, "eps":0.1 }},
                 {"B": { "q":-1.0, "r":2.0, "eps":0.05 }},
                 {"C": { "r":1.0, "mu":[2,0,0] }} ]})"_json;

    atoms = j["atomlist"].get<decltype(atoms)>();

    auto u = pairpotential::makePairPotential<FunctorPotential>(R"(
                { "default": [ { "coulomb" : {"epsr": 80.0, "type": "plain"} } ],
                  "A B" : [
                    { "coulomb" : {"epsr": 80.0, "type": "plain"} },
                    { "wca" : {"mixing": "LB"} }
                  ],
                  "C C" : [ { "hardsphere" : {} } ] })"_json);

    auto coulomb = pairpotential::makePairPotential<Coulomb>(R"({ "coulomb": {"epsr": 80.0} } )"_json);
    auto wca = pairpotential::makePairPotential<WeeksChandlerAndersen>(R"({ "wca" : {"mixing": "LB"} })"_json);

    Particle a = atoms[0];
    Particle b = atoms[1];
    Particle c = atoms[2];
    Point r = {2, 0, 0};
    auto r2 = r.squaredNorm();
    CHECK_EQ(u(a, a, r2, r), Approx(coulomb(a, a, r2, r)));
    CHECK_EQ(u(b, b, r2, r), Approx(coulomb(b, b, r2, r)));
    CHECK_EQ(u(a, b, r2, r), Approx(coulomb(a, b, r2, r) + wca(a, b, r2, r)));
    CHECK_EQ(u(c, c, (r * 1.01).squaredNorm(), r * 1.01), 0);
    CHECK_EQ(u(c, c, (r * 0.99).squaredNorm(), r * 0.99), pc::infty);

    SUBCASE("selfEnergy() - monopole") {
        // let's check that the self energy gets properly transferred to the functor potential
        const auto j = R"(
                {"default": [{ "coulomb" : {"epsr": 80.0, "type": "qpotential", "cutoff":20, "order":4} }]})"_json;

        auto functor = pairpotential::makePairPotential<FunctorPotential>(j);
        auto galore = pairpotential::makePairPotential<NewCoulombGalore>(j["default"].at(0));
        CHECK((functor.selfEnergy != nullptr));
        CHECK((galore.selfEnergy != nullptr));
        CHECK_EQ(functor.selfEnergy(a), Approx(galore.selfEnergy(a)));
    }

    SUBCASE("selfEnergy() - multipole") {
        // now test w. dipolar particles
        const auto j = R"(
                {"default": [{ "multipole" : {"epsr": 80.0, "type": "qpotential", "cutoff":20, "order":4} }]})"_json;

        auto functor = pairpotential::makePairPotential<FunctorPotential>(j);
        auto multipole = pairpotential::makePairPotential<Multipole>(j["default"].at(0));
        CHECK((functor.selfEnergy != nullptr));
        CHECK((multipole.selfEnergy != nullptr));
        CHECK_EQ(functor.selfEnergy(a), Approx(multipole.selfEnergy(a))); // q=1, mu=0
        CHECK_EQ(functor.selfEnergy(c), Approx(multipole.selfEnergy(c))); // q=0, mu=2
    }
}

// =============== SplinedPotential ===============

SplinedPotential::KnotData::KnotData(const base &b) : base(b) {}

/**
 * @param stream output stream
 * @param id1 fist atom id
 * @param id2 second atom id
 *
 * Stream splined and exact energy as a function of particle-particle separation to output stream
 */
void SplinedPotential::streamPairPotential(std::ostream& stream, const size_t id1, const size_t id2) {
    if (!stream || atoms.at(id1).implicit || atoms.at(id2).implicit) {
        return;
    }
    stream << "# r u_splined/kT u_exact/kT\n";
    const auto particle_1 = static_cast<Particle>(Faunus::atoms.at(id1));
    const auto particle_2 = static_cast<Particle>(Faunus::atoms.at(id2));
    const auto rmax = std::sqrt(matrix_of_knots(id1, id2).rmax2);
    for (auto r : arange(dr, rmax, dr)) {
        stream << fmt::format("{:.6E} {:.6E} {:.6E}\n", r, operator()(particle_1, particle_2, r* r, {r, 0, 0}),
                              FunctorPotential::operator()(particle_1, particle_2, r* r, {r, 0, 0}));
    }
}

/**
 * For each pair of atom types a file is created containing
 * the exact and splined pair potential as a function of distance
 */
void SplinedPotential::savePotentials() {
    for (size_t i = 0; i < Faunus::atoms.size(); ++i) { // loop over atom types
        for (size_t j = 0; j <= i; ++j) {               // and build matrix of spline data (knots) for each pair
            auto filename = fmt::format("{}-{}_tabulated.dat", Faunus::atoms.at(i).name, Faunus::atoms.at(j).name);
            if (auto stream = std::ofstream(filename); stream) {
                streamPairPotential(stream, i, j);
            }
        }
    }
}

/**
 * @param i Atom type index
 * @param j Atom type index
 * @param energy_threshold Absolute maximum energy used to determine min. distance for spline
 * @param rmin Initial starting point for rmin
 * @return Minimum distance for splining
 *
 * For repulsive potentials (at short sep.), the threshold is the maximum energy splined
 * For attactive potentials (at short sep.), the threshold is the minimum energy splined
 */
double SplinedPotential::findLowerDistance(int i, int j, double energy_threshold, double rmin) {
    assert(rmin > 0);
    Particle particle1 = Faunus::atoms.at(i);
    Particle particle2 = Faunus::atoms.at(j);
    int num_iterations = 0;
    while (rmin >= dr) {
        if (num_iterations++ == max_iterations) {
            throw std::runtime_error("Pair potential spline error: cannot determine minimum distance");
        }
        double u = std::fabs(FunctorPotential::operator()(particle1, particle2, rmin *rmin, {rmin, 0, 0}));
        if (u > energy_threshold * 1.1) {
            rmin += dr;
        } else if (u < energy_threshold / 1.1) {
            rmin -= dr;
        } else {
            break;
        }
    }
    assert(rmin >= 0);
    return rmin;
}

/**
 * @param i Atom type index
 * @param j Atom type index
 * @param energy_threshold
 * @param rmax Initial guess for max. distance
 * @return Maximum splining distance
 *
 * All pair potential are assumed to approach zero at large separations.
 * This function increases `rmin` until the absolute energy is lower than
 * `energy_threshold`.
 */
double SplinedPotential::findUpperDistance(int i, int j, double energy_threshold, double rmax) {
    assert(rmax > 0);
    Particle particle1 = Faunus::atoms.at(i);
    Particle particle2 = Faunus::atoms.at(j);
    int num_iterations = 0;
    while (rmax >= dr) {
        if (num_iterations++ == max_iterations) {
            throw std::runtime_error("Pair potential spline error: cannot determine maximum distance");
        }
        double u = FunctorPotential::operator()(particle1, particle2, rmax *rmax, {rmax, 0, 0});
        if (std::fabs(u) > energy_threshold) {
            rmax += dr;
        } else {
            break;
        }
    }
    return rmax;
}

void SplinedPotential::from_json(const json &js) {
    FunctorPotential::from_json(js);
    if (!isotropic) {
        throw std::runtime_error("Cannot spline anisotropic potentials");
    }
    spline.setTolerance(js.value("utol", 1e-3), js.value("ftol", 1e-2));
    hardsphere_repulsion = js.value("hardsphere", false);
    double energy_at_rmin = js.value("u_at_rmin", 20);
    double energy_at_rmax = js.value("u_at_rmax", 1e-6);

    faunus_logger->trace("Pair potential spline tolerance = {} kT", js.value("utol", 1e-5));

    for (size_t i = 0; i < Faunus::atoms.size(); ++i) { // loop over atom types
        for (size_t j = 0; j <= i; ++j) {               // and build matrix of spline data (knots) for each pair
            if (atoms[i].implicit || atoms[j].implicit) {
                continue;
            }
            double rmin = 0.5 * (Faunus::atoms[i].sigma + Faunus::atoms[j].sigma);
            double rmax = js.value("rmax", rmin * 10);
            if (auto it = js.find("cutoff_g2g"); it != js.end()) {
                if (it->is_number()) {
                    rmax = it->get<double>();
                } else if (it->is_object()) {
                    rmax = it->at("default").get<double>();
                }
            }
            rmin = findLowerDistance(i, j, energy_at_rmin, rmin);
            rmax = findUpperDistance(i, j, energy_at_rmax, rmax);
            assert(rmin < rmax);
            createKnots(i, j, rmin, rmax);
        }
    }
    if (js.value("to_disk", false)) {
        savePotentials();
    }
}

SplinedPotential::SplinedPotential(const std::string &name) : FunctorPotential(name) {}

/**
 * @param i Atom index
 * @param j Atom index
 * @param rmin Minimum splining distance
 * @param rmax Maximum splining distance
 */
void SplinedPotential::createKnots(int i, int j, double rmin, double rmax) {
    Particle particle1 = Faunus::atoms.at(i);
    Particle particle2 = Faunus::atoms.at(j);
    KnotData knotdata = spline.generate(
        [&](double r_squared) {
            return FunctorPotential::operator()(particle1, particle2, r_squared, {0, 0, 0});
        },
        rmin * rmin, rmax * rmax); // spline along r^2

    // if set, hard-sphere repulsion (infinity) is used IF the potential is repulsive below rmin
    knotdata.hardsphere_repulsion = hardsphere_repulsion;
    if (spline.eval(knotdata, knotdata.rmin2 + dr) < 0) { // disable hard sphere
        knotdata.hardsphere_repulsion = false;            // repulsion for attractive potentials
    }
    if (knotdata.hardsphere_repulsion) {
        faunus_logger->trace("Hardsphere repulsion enabled for {}-{} spline", Faunus::atoms.at(i).name,
                             Faunus::atoms.at(j).name);
    }
    matrix_of_knots.set(i, j, knotdata); // register knots for the pair

    double max_error = 0.0; // maximum absolute error of the spline along r
    for (const auto r : arange(rmin + dr, rmax, dr)) {
        auto error = std::fabs(operator()(particle1, particle2, r* r, {r, 0, 0}) -
                                 FunctorPotential::operator()(particle1, particle2, r* r, {r, 0, 0}));
        max_error = std::max(error, max_error);
    }
    faunus_logger->trace(
        "{}-{} interaction splined between [{:6.2f}:{:6.2f}] {} using {} knots w. maximum absolute error of {:.1E} kT",
        Faunus::atoms[i].name, Faunus::atoms[j].name, rmin, rmax, unicode::angstrom, knotdata.numKnots(), max_error);
}

// =============== NewCoulombGalore ===============

void NewCoulombGalore::setSelfEnergy() {
    selfEnergy = [lB = bjerrum_length, self_energy = pot.selfEnergyFunctor](const Particle &p) {
        return lB * self_energy({p.charge * p.charge, 0.0});
    }; // expose self-energy as a functor in potential base class
}

NewCoulombGalore::NewCoulombGalore(const std::string& name)
    : PairPotential(name) {
    setSelfEnergy();
}

TEST_CASE("[Faunus] NewCoulombGalore") {
    Particle a, b;
    a.charge = 1.0;
    b.charge = -1.0;
    a.pos = {0, 0, 0};
    b.pos = {0, 0, 7};
    Point b_towards_a = a.pos - b.pos;
    auto pot = makePairPotential<NewCoulombGalore>(R"({"epsr": 80, "type": "plain"})"_json);
    Point force_on_a = pot.force(a, b, b_towards_a.squaredNorm(), b_towards_a);
    CHECK_EQ(force_on_a.x(), doctest::Approx(0));
    CHECK_EQ(force_on_a.y(), doctest::Approx(0));
    CHECK_EQ(force_on_a.z(), doctest::Approx(0.1429734149)); // attraction -> positive direction expected
}

void NewCoulombGalore::from_json(const json &j) {
    using namespace ::CoulombGalore; // namespace for external CoulombGalore library
    const auto relative_dielectric_constant = j.at("epsr").get<double>();
    bjerrum_length = pc::bjerrumLength(relative_dielectric_constant);
    const auto electrolyte = Faunus::makeElectrolyte(j);
    pot.setTolerance(j.value("utol", 0.005 / bjerrum_length));
    const auto method = j.at("type").get<std::string>();
    if (method == "yukawa") {
        if (json _j(j); _j.value("shift", false)) { // zero energy and force at cutoff
            faunus_logger->debug("energy and force shifted yukawa uses the 'poisson' scheme with C=1 and D=1");
            _j["type"] = "poisson";
            _j["C"] = 1;
            _j["D"] = 1;
            _j["debyelength"] = electrolyte.value().debyeLength(bjerrum_length);
            pot.spline<::CoulombGalore::Poisson>(_j);
        } else { // non-shifted yukawa equals `plain` with exponential screening
            if (_j.contains("cutoff")) {
                throw ConfigurationError("unexpected 'cutoff' for non-shifted yukawa which is always infinity");
            }
            _j["type"] = "plain";
            _j["debyelength"] = electrolyte.value().debyeLength(bjerrum_length);
            pot.spline<::CoulombGalore::Plain>(_j);
        }
    } else if (method == "plain") {
        if (j.contains("cutoff")) {
            throw ConfigurationError("unexpected cutoff for plain: it's *always* infinity");
        }
        pot.spline<::CoulombGalore::Plain>(j);
    } else if (method == "qpotential") {
        pot.spline<::CoulombGalore::qPotential>(j);
    } else if (method == "wolf") {
        pot.spline<::CoulombGalore::Wolf>(j);
    } else if (method == "poisson") {
        pot.spline<::CoulombGalore::Poisson>(j);
    } else if (method == "fanourgakis") {
        pot.spline<::CoulombGalore::Fanourgakis>(j);
    } else if (method == "zahn") {
        pot.spline<::CoulombGalore::Zahn>(j);
    } else if (method == "fennell") {
        pot.spline<::CoulombGalore::Fennell>(j);
    } else if (method == "zerodipole") {
        pot.spline<::CoulombGalore::ZeroDipole>(j);
    } else if (method == "ewald") {
        auto _j(j);
        if (electrolyte) {
            _j["debyelength"] = electrolyte.value().debyeLength(bjerrum_length);
        }
        pot.spline<::CoulombGalore::Ewald>(_j);
    } else if (method == "reactionfield") {
        pot.spline<::CoulombGalore::ReactionField>(j);
    } else {
        throw ConfigurationError("unknown type '{}'", method);
    }
    faunus_logger->info("splitting function for '{}' splined with {} knots", method, pot.numKnots().at(0));
    setSelfEnergy();
}

void NewCoulombGalore::to_json(json &j) const {
    pot.to_json(j);
    j["lB"] = bjerrum_length;
}
const CoulombGalore::Splined& NewCoulombGalore::getCoulombGalore() const { return pot; }
[[maybe_unused]] double NewCoulombGalore::dielectric_constant(double M2V) { return pot.calc_dielectric(M2V); }

// =============== Multipole ===============

Multipole::Multipole(const std::string &name) : NewCoulombGalore(name) {
    isotropic = false; // this potential is angular dependent
    setSelfEnergy();
}

void Multipole::setSelfEnergy() {
    selfEnergy = [lB = bjerrum_length, self_energy = pot.selfEnergyFunctor](const Particle& p) {
        double mu_x_mu = 0;     // dipole-dipole product
        if (p.hasExtension()) { // only access dipole if the particle has extended properties
            mu_x_mu = p.getExt().mulen * p.getExt().mulen;
        }
        return lB * self_energy({p.charge * p.charge, mu_x_mu});
    }; // expose self-energy as a functor in potential base class
}

TEST_CASE("[Faunus] Dipole-dipole interactions") {
    using doctest::Approx;
    json j = R"({ "atomlist" : [
                 {"A": { "mu":[1.0,0.0,0.0], "mulen":3.0 }},
                 {"B": { "mu":[0.0,1.0,0.0], "mulen":3.0 }},
                 {"C": { "mu":[1.0,1.0,0.0] }} ]})"_json;
    atoms = j["atomlist"].get<decltype(atoms)>();

    auto u = pairpotential::makePairPotential<FunctorPotential>(R"(
                { "default": [ { "multipole" : {"epsr": 1.0, "type": "plain"} } ] } )"_json);

    auto dipoledipole = pairpotential::makePairPotential<Multipole>(R"({"epsr": 1.0, "type": "plain"})"_json);

    Particle a = atoms[0];
    Particle b = atoms[1];
    Particle c = atoms[2];
    Point r = {2, 0, 0};
    double r2 = r.squaredNorm();
    CHECK_EQ(u(a, a, r2, r),
          Approx(dipoledipole(a, a, r2,
                              r))); // interaction between two parallell dipoles, directed parallell to their seperation
    CHECK_EQ(u(b, b, r2, r),
          Approx(dipoledipole(
              b, b, r2, r))); // interaction between two parallell dipoles, directed perpendicular to their seperation
    CHECK_EQ(u(a, b, r2, r), Approx(dipoledipole(a, b, r2, r))); // interaction between two perpendicular dipoles
    CHECK_EQ(u(a, a, r2, r), -2.25 * dipoledipole.bjerrum_length);
    CHECK_EQ(u(b, b, r2, r), 1.125 * dipoledipole.bjerrum_length);
    CHECK_EQ(u(a, c, r2, r), -0.75 * dipoledipole.bjerrum_length);
    CHECK_EQ(u(b, c, r2, r), 0.375 * dipoledipole.bjerrum_length);
    CHECK_EQ(u(a, b, r2, r), 0);

    r = {3, 0, 0};
    r2 = 3 * 3;
    CHECK_EQ(u(a, a, r2, r),
          Approx(dipoledipole(a, a, r2,
                              r))); // interaction between two parallell dipoles, directed parallell to their seperation
    CHECK_EQ(u(b, b, r2, r),
          Approx(dipoledipole(
              b, b, r2, r))); // interaction between two parallell dipoles, directed perpendicular to their seperation
    CHECK_EQ(u(a, b, r2, r), Approx(dipoledipole(a, b, r2, r))); // interaction between two perpendicular dipoles
    CHECK_EQ(u(a, a, r2, r), Approx(-(2.0 / 3.0) * dipoledipole.bjerrum_length));
    CHECK_EQ(u(b, b, r2, r), (1.0 / 3.0) * dipoledipole.bjerrum_length);
    CHECK_EQ(u(a, c, r2, r), Approx(-2.0 / 9.0 * dipoledipole.bjerrum_length));
    CHECK_EQ(u(b, c, r2, r), 1.0 / 9.0 * dipoledipole.bjerrum_length);
    CHECK_EQ(u(a, b, r2, r), 0);
}

// =============== WeeksChandlerAndersen ===============

WeeksChandlerAndersen::WeeksChandlerAndersen(const std::string &name, const std::string &cite,
                                             CombinationRuleType combination_rule)
    : LennardJones(name, cite, combination_rule) {}

TEST_CASE("[Faunus] WeeksChandlerAndersen") {
    SUBCASE("JSON initilization") {
        atoms = R"([{"A": {"sigma":2.0, "eps":0.9}},
                    {"B": {"sigma":8.0, "eps":0.1}}])"_json.get<decltype(atoms)>();
        Particle a = atoms[0], b = atoms[1];

        CHECK_THROWS_AS(makePairPotential<WeeksChandlerAndersen>(R"({"mixing": "unknown"})"_json), std::runtime_error);
        CHECK_NOTHROW(makePairPotential<WeeksChandlerAndersen>(R"({})"_json));
        CHECK_NOTHROW(makePairPotential<WeeksChandlerAndersen>(
            R"({"mixing": "LB", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json));

        SUBCASE("Missing coefficient") {
            auto wca = makePairPotential<WeeksChandlerAndersen>(R"({"mixing": "LB", "sigma": "sigma_wca"})"_json);
            CHECK_EQ(std::isnan(wca(a, a, 10.0 * 10.0, {0, 0, 10.0})), true);
            CHECK_EQ(std::isnan(wca(a, a, 1.0 * 1.0, {0, 0, 1.0})), true);
        }
    }
    SUBCASE("JSON initilization custom coefficient names") {
        atoms = R"([{"A": {"sigma_wca":2.0, "eps":0.9}},
                    {"B": {"sigma_wca":8.0, "eps":0.1}}])"_json.get<decltype(atoms)>();
        Particle a = atoms[0], b = atoms[1];

        CHECK_NOTHROW(makePairPotential<WeeksChandlerAndersen>(R"({"mixing": "LB", "sigma": "sigma_wca"})"_json));
        // Shall throw after non-default potentials are created properly,
        // i.e., not needed pairs are not evaluated at all for the matrices
        // CHECK_THROWS_AS_MESSAGE(WeeksChandlerAndersen wca = R"({"mixing": "LB"})"_json, std::runtime_error,
        //                         "unknown atom property");
        // CHECK_THROWS_AS_MESSAGE(WeeksChandlerAndersen wca = R"({"mixing": "LB", "sigma": "unknown"})"_json,
        //                         std::runtime_error, "unknown atom property");
        // different atom and custom coefficient names are not allowed
        // CHECK_THROWS_AS_MESSAGE(
        //     WeeksChandlerAndersen wca =
        //         R"({"mixing": "LB", "sigma": "sigma_wca", "custom": [{"A B": {"eps": 0.5, "sigma": 8}}]})"_json,
        //     std::runtime_error, "unknown atom property");
    }
    SUBCASE("JSON serialization") {
        auto wca =
            pairpotential::makePairPotential<WeeksChandlerAndersen>(R"({"mixing": "LB", "sigma": "sigma_wca"})"_json);
        json j = wca;
        json &j_wca(j["wca"]);
        CHECK_EQ(j_wca["mixing"], "lorentz_berthelot");
        CHECK_EQ(j_wca["sigma"], "sigma_wca");
        CHECK_EQ(j_wca["eps"], "eps");
    }
}

PairPotentialException::PairPotentialException(const std::string& msg)
    : std::runtime_error(msg){}

void CosAttractMixed::initPairMatrices() {
    auto comb_distance = PairMixer::getCombinator(combination_rule, PairMixer::CoefficientType::SIGMA);
    auto comb_epsilon = PairMixer::getCombinator(combination_rule, PairMixer::CoefficientType::EPSILON);

    switching_distance = PairMixer(extract_rc, comb_distance, &PairMixer::modIdentity).createPairMatrix(atoms, *custom_pairs);
    switching_width = PairMixer(extract_wc, comb_distance, &PairMixer::modIdentity).createPairMatrix(atoms, *custom_pairs);
    epsilon = PairMixer(extract_eps, comb_epsilon, &PairMixer::modIdentity).createPairMatrix(atoms, *custom_pairs);

    faunus_logger->trace("Pair matrices for {} sigma ({}×{}) and epsilon ({}×{}) created using {} custom pairs.", name,
                         switching_distance->rows(), switching_distance->cols(), epsilon->rows(), epsilon->cols(), custom_pairs->size());
}

void CosAttractMixed::extractorsFromJson(const json& j) {
    auto name = j.value("rc", "rc");
    json_extra_params["rc"] = name;
    extract_rc = [=](auto& a) { return a.at(name) * 1.0_angstrom; };

    name = j.value("wc", "wc");
    json_extra_params["wc"] = name;
    extract_wc = [=](auto& a) { return a.at(name) * 1.0_angstrom; };

    name = j.value("eps", "eps");
    json_extra_params["eps"] = name;
    extract_eps = [=](auto& a) { return a.at(name) * 1.0_kJmol; };
}

CosAttractMixed::CosAttractMixed(const std::string& name, const std::string& cite, CombinationRuleType combination_rule)
    : MixerPairPotentialBase(name, cite, combination_rule) {}

double CosAttractMixed::cutOffSquared(AtomData::index_type id1, AtomData::index_type id2) const {
    const auto rc = (*switching_distance)(id1, id2);
    const auto wc = (*switching_width)(id1, id2);
    return (rc + wc) * (rc + wc);
}
} // namespace Faunus::pairpotential
