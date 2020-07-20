#include "potentials.h"
#include "multipole.h"
#include "units.h"
#include "spdlog/spdlog.h"
#include <coulombgalore.h>

namespace Faunus {
namespace Potential {

// =============== PairMixer ===============

TCombinatorFunc PairMixer::getCombinator(CombinationRuleType combination_rule, CoefficientType coefficient) {
    TCombinatorFunc combinator;
    switch (combination_rule) {
    case COMB_UNDEFINED:
        combinator = &combUndefined;
        break;
    case COMB_ARITHMETIC:
        combinator = &combArithmetic;
        break;
    case COMB_GEOMETRIC:
        combinator = &combGeometric;
        break;
    case COMB_LORENTZ_BERTHELOT:
        switch (coefficient) {
        case COEF_SIGMA:
            combinator = &combArithmetic;
            break;
        case COEF_EPSILON:
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

TPairMatrixPtr PairMixer::createPairMatrix(const std::vector<AtomData> &atoms) {
    size_t n = atoms.size(); // number of atom types
    TPairMatrixPtr matrix = std::make_shared<TPairMatrix>(n, n);
    for (auto &i : atoms) {
        for (auto &j : atoms) {
            if(i.implicit || j.implicit) {
                // implicit atoms are ignored as the missing properties, e.g., sigma and epsilon, might raise errors
                (*matrix)(i.id(), j.id()) = combUndefined();
            } else if (i.id() == j.id()) {
                // if the combinator is "undefined" the homogeneous interaction is still well defined
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
    for (auto &i : interactions) {
        if (i.atom_id[0] >= 0 && i.atom_id[1] >= 0 && i.atom_id[0] < dimension && i.atom_id[1] < dimension) {
            // interaction is always symmetric
            (*matrix)(i.atom_id[0], i.atom_id[1]) = (*matrix)(i.atom_id[1], i.atom_id[0]) =
                modifier(extractor(i.interaction));
        } else {
            throw std::range_error("atomtype index out of range");
        }
    }
    return matrix;
}

void from_json(const json &j, CustomInteractionData &c) {
    if (!j.is_object() || j.size() != 1) {
        throw ConfigurationError("invalid JSON for custom interaction parameters");
    }
    auto j_item = j.items().begin();
    try {
        const auto atom_names = words2vec<std::string>(j_item.key());
        const auto atom_ids = names2ids(atoms, atom_names);
        if (atom_ids.size() != c.atom_id.size()) {
            faunus_logger->error("Custom interaction parameters require exactly {} space-separated atoms: {}.",
                                 c.atom_id.size(), vec2words(atom_names));
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
    j = {{vec2words(atom_names), interaction.interaction}};
}

void from_json(const json &j, std::vector<CustomInteractionData> &interactions) {
    if(j.is_array()) {
        for (const auto j_pair: j) {
            interactions.push_back(j_pair);
        }
    } else if(j.is_object()) {
        for (const auto j_kv: j.items()) {
            interactions.push_back(json {{j_kv.key(), j_kv.value()}});
        }
    } else {
        throw ConfigurationError("invalid JSON for custom interaction parameters");
    }
}

// =============== PairPotentialBase ===============

PairPotentialBase::PairPotentialBase(const std::string &name, const std::string &cite, bool isotropic)
    : name(name), cite(cite), isotropic(isotropic) {}

Point PairPotentialBase::force(const Particle &, const Particle &, double, const Point &) const {
    assert(false && "We should never reach this point!");
    throw(std::logic_error("Force computation not implemented for this setup!"));
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
        if (j.count("mixing") == 1) {
            json mixing = j.at("mixing");
            combination_rule = mixing.get<CombinationRuleType>();
            if (combination_rule == COMB_UNDEFINED && mixing != "undefined") {
                // an ugly hack because the first pair in the json ↔ enum mapping is silently selected by default
                throw PairPotentialException("unknown combination rule " + mixing.get<std::string>());
            }
        }
        if (j.count("custom") == 1) {
            // *custom_pairs = j["custom"]; // does not work, perhaps as from_json is also a method (a namespace conflict)
            Potential::from_json(j["custom"], *custom_pairs);
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

void RepulsionR3::from_json(const json &j) {
    f = j.value("prefactor", 1.0);
    e = j.value("lj-prefactor", 1.0);
    s = j.value("sigma", 1.0);
}

void RepulsionR3::to_json(json &j) const { j = {{"prefactor", f}, {"lj-prefactor", e}, {"sigma", s}}; }

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

void Coulomb::to_json(json &j) const {
    j["epsr"] = pc::relativeDielectricFromBjerrumLength(lB);
    j["lB"] = lB;
}

void Coulomb::from_json(const json &j) { lB = pc::bjerrumLength(j.at("epsr")); }

void DipoleDipole::to_json(json &j) const {
    j["epsr"] = pc::relativeDielectricFromBjerrumLength(lB);
    j["lB"] = lB;
}

void DipoleDipole::from_json(const json &j) { lB = pc::bjerrumLength(j.at("epsr")); }

void FENE::from_json(const json &j) {
    k = j.at("stiffness");
    r02 = std::pow(double(j.at("maxsep")), 2);
    r02inv = 1 / r02;
}

void FENE::to_json(json &j) const { j = {{"stiffness", k}, {"maxsep", std::sqrt(r02)}}; }

void to_json(json &j, const PairPotentialBase &base) {
    base.name.empty() ? base.to_json(j) : base.to_json(j[base.name]);
}

void from_json(const json &j, PairPotentialBase &base) {
    try {
        if (not base.name.empty()) {
            if (j.count(base.name) == 1) {
                base.from_json(j.at(base.name));
                return;
            }
        }
        base.from_json(j);
    } catch (std::exception &e) {
        throw std::runtime_error(e.what() + usageTip[base.name]);
    }
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

double SASApotential::area(double R, double r, double d_squared) const {
    R += proberadius;
    r += proberadius;
    double area = 4 * pc::pi * (R * R + r * r); // full volume of both spheres
    double offset = (shift ? area : 0);
    if (d_squared > (R + r) * (R + r))
        return area - offset;
    if (r > R)
        std::swap(r, R);
    double d = sqrt(d_squared);
    if (d + r <= R)
        return 4 * pc::pi * R * R - offset;          // full volume of biggest sphere
    double h1 = (r - R + d) * (r + R - d) / (2 * d); // height of spherical caps
    double h2 = (R - r + d) * (R + r - d) / (2 * d); // comprising intersecting lens
    return area - 2 * pc::pi * (R * h1 + r * h2) - offset;
}

void CustomPairPotential::from_json(const json &j) {
    Rc2 = j.value("cutoff", pc::infty);
    Rc2 = Rc2 * Rc2;
    jin = j;
    auto &_j = jin["constants"];
    if (_j == nullptr)
        _j = json::object();
    _j["e0"] = pc::e0;
    _j["kB"] = pc::kB;
    _j["kT"] = pc::kT();
    _j["Nav"] = pc::Nav;
    _j["Rc"] = std::sqrt(Rc2);
    _j["T"] = pc::temperature;
    expr.set(jin, {{"r", &d->r}, {"q1", &d->q1}, {"q2", &d->q2}, {"s1", &d->s1}, {"s2", &d->s2}});
}

void CustomPairPotential::to_json(json &j) const {
    j = jin;
    if (std::isfinite(Rc2))
        j["cutoff"] = std::sqrt(Rc2);
}
CustomPairPotential::CustomPairPotential(const std::string &name)
    : PairPotentialBase(name), d(std::make_shared<Data>()) {}

// =============== Dummy ===============

Dummy::Dummy() { name = "dummy"; }
void Dummy::from_json(const json &) {}
void Dummy::to_json(json &) const {}

// =============== LennardJones ===============

void LennardJones::initPairMatrices() {
    const TCombinatorFunc comb_sigma = PairMixer::getCombinator(combination_rule, PairMixer::COEF_SIGMA);
    const TCombinatorFunc comb_epsilon = PairMixer::getCombinator(combination_rule, PairMixer::COEF_EPSILON);

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
    extract_sigma = [sigma_name](const InteractionData &a) -> double { return a.get(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const InteractionData &a) -> double { return a.get(epsilon_name) * 1.0_kJmol; };
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
    extract_sigma = [sigma_name](const InteractionData &a) -> double { return a.get(sigma_name) * 1.0_angstrom; };
}

// =============== Hertz ===============

void Hertz::initPairMatrices() {
    const TCombinatorFunc comb_diameter = PairMixer::getCombinator(combination_rule, PairMixer::COEF_SIGMA);
    const TCombinatorFunc comb_epsilon = PairMixer::getCombinator(combination_rule, PairMixer::COEF_EPSILON);

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
    extract_sigma = [sigma_name](const InteractionData &a) -> double { return a.get(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const InteractionData &a) -> double { return a.get(epsilon_name) * 1.0_kJmol; };
}

// =============== SquareWell ===============

void SquareWell::initPairMatrices() {
    const TCombinatorFunc comb_diameter = PairMixer::getCombinator(combination_rule, PairMixer::COEF_SIGMA);
    const TCombinatorFunc comb_depth = PairMixer::getCombinator(combination_rule, PairMixer::COEF_EPSILON);

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
    extract_sigma = [sigma_name](const InteractionData &a) -> double { return a.get(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const InteractionData &a) -> double { return a.get(epsilon_name) * 1.0_kJmol; };
}

// =============== Polarizability ===============

void Polarizability::from_json(const json &j) {
    epsr = j.at("epsr").get<double>();
    double lB = pc::bjerrumLength(epsr);
    for (auto &i : atoms) {
        for (auto &j : atoms) {
            m_neutral->set(i.id(), j.id(), -3 * i.alphax * pow(0.5 * i.sigma, 3) * j.alphax * pow(0.5 * j.sigma, 3));
            m_charged->set(i.id(), j.id(),
                           -lB / 2 *
                               (pow(i.charge, 2) * j.alphax * pow(0.5 * j.sigma, 3) +
                                pow(j.charge, 2) * i.alphax * pow(0.5 * i.sigma, 3)));
        }
    }
}

//----------------- FunctorPotential ---------------------

void FunctorPotential::registerSelfEnergy(PairPotentialBase *pot) {
    if (pot->selfEnergy) {
        if (not selfEnergy) // no self energy is defined
            selfEnergy = pot->selfEnergy;
        else // accumulate self energies
            selfEnergy = [pot = pot, &selfEnergy = selfEnergy](const Particle &p) {
                return pot->selfEnergy(p) + selfEnergy(p);
            };
        faunus_logger->debug("Added selfEnergy function from {} to {}", pot->name, name);
    } else
        faunus_logger->trace("Failed to register non-defined selfEnergy() for {}", pot->name);
}

FunctorPotential::uFunc FunctorPotential::combineFunc(json &j) {
    uFunc u = [](const Particle &, const Particle &, double, const Point &) { return 0.0; };
    if (j.is_array()) {
        for (auto &i : j) { // loop over all defined potentials in array
            if (i.is_object() and (i.size() == 1)) {
                for (auto it : i.items()) {
                    uFunc _u = nullptr;
                    try {
                        if (it.key() == "custom")
                            _u = CustomPairPotential() = it.value();

                        // add Coulomb potential and self-energy
                        // terms if not already added
                        else if (it.key() == "coulomb") { // temporary name
                            std::get<0>(potlist).from_json(it.value()); // initialize w. json object
                            std::get<0>(potlist).to_json(it.value());   // write back to json object with added values
                            _u = std::get<0>(potlist);
                            if (not have_monopole_self_energy) {
                                registerSelfEnergy(&std::get<0>(potlist));
                                have_monopole_self_energy = true;
                            }
                        } else if (it.key() == "cos2")
                            _u = std::get<1>(potlist) = i;
                        else if (it.key() == "polar")
                            _u = std::get<2>(potlist) = i;
                        else if (it.key() == "hardsphere")
                            _u = std::get<3>(potlist) = i;
                        else if (it.key() == "lennardjones")
                            _u = std::get<4>(potlist) = i;
                        else if (it.key() == "repulsionr3")
                            _u = std::get<5>(potlist) = i;
                        else if (it.key() == "sasa")
                            _u = std::get<6>(potlist) = i;
                        else if (it.key() == "wca")
                            _u = std::get<7>(potlist) = i;
                        else if (it.key() == "pm")
                            _u = std::get<8>(potlist) = it.value();
                        else if (it.key() == "pmwca")
                            _u = std::get<9>(potlist) = it.value();
                        else if (it.key() == "hertz")
                            _u = std::get<10>(potlist) = i;
                        else if (it.key() == "squarewell")
                            _u = std::get<11>(potlist) = i;
                        else if (it.key() == "dipoledipole") {
                            faunus_logger->error("'{}' is deprecated, use 'multipole' instead", it.key());
                        } else if (it.key() == "stockmayer") {
                            faunus_logger->error("'{}' is deprecated, use 'lennardjones'+'multipole' instead",
                                                 it.key());
                        } else if (it.key() == "multipole") {
                            std::get<12>(potlist).from_json(it.value()); // init from json
                            std::get<12>(potlist).to_json(it.value());   // write back added info to json
                            _u = std::get<12>(potlist);
                            isotropic = false;                         // potential is now angular dependent
                            if (not have_dipole_self_energy) {
                                registerSelfEnergy(&std::get<12>(potlist));
                                have_dipole_self_energy = true;
                            }
                        }
                        // place additional potentials here...
                    } catch (std::exception &e) {
                        throw std::runtime_error(it.key() + ": " + e.what() + usageTip[it.key()]);
                    }

                    if (_u != nullptr) // if found, sum them into new function object
                        u = [u, _u](const Particle &a, const Particle &b, double r2, const Point &r) {
                            return u(a, b, r2, r) + _u(a, b, r2, r);
                        };
                    else
                        throw std::runtime_error("unknown potential: " + it.key());
                }
            }
        }
    } else
        throw std::runtime_error("dictionary of potentials required");

    return u;
}

void FunctorPotential::to_json(json &j) const {
    j = _j;
    j["selfenergy"] = {{"monopole", have_monopole_self_energy}, {"dipole", have_dipole_self_energy}};
}

void FunctorPotential::from_json(const json &j) {
    have_monopole_self_energy = false;
    have_dipole_self_energy = false;
    _j = j;
    umatrix = decltype(umatrix)(atoms.size(), combineFunc(_j.at("default")));
    for (auto it = _j.begin(); it != _j.end(); ++it) {
        auto atompair = words2vec<std::string>(it.key()); // is this for a pair of atoms?
        if (atompair.size() == 2) {
            auto ids = names2ids(atoms, atompair);
            umatrix.set(ids[0], ids[1], combineFunc(it.value()));
        }
    }
}

FunctorPotential::FunctorPotential(const std::string &name) : PairPotentialBase(name) {}

SplinedPotential::KnotData::KnotData(const base &b) : base(b) {}

/**
 * @param stream output stream
 * @param i fist atom index
 * @param j second atom index
 *
 * Stream splined and exact energy as a function of particle-particle separation to output stream
 */
void SplinedPotential::stream_pair_potential(std::ostream &stream, int i, int j) {
    if (stream) {
        if (!atoms.at(i).implicit && !atoms.at(j).implicit) {
            Particle particle1 = Faunus::atoms.at(i);
            Particle particle2 = Faunus::atoms.at(j);
            stream << "# r u_splined/kT u_exact/kT\n";
            double rmax = sqrt(matrix_of_knots(i, j).rmax2);
            for (double r = dr; r < rmax; r += dr) {
                stream << fmt::format("{:.6E} {:.6E} {:.6E}\n", r, operator()(particle1, particle2, r *r, {r, 0, 0}),
                                      FunctorPotential::operator()(particle1, particle2, r *r, {r, 0, 0}));
            }
        }
    }
}

/**
 * For each pair of atom types a file is created containing
 * the exact and splined pair potential as a function of distance
 */
void SplinedPotential::save_potentials() {
    for (size_t i = 0; i < Faunus::atoms.size(); ++i) { // loop over atom types
        for (size_t j = 0; j <= i; ++j) {               // and build matrix of spline data (knots) for each pair
            auto filename = fmt::format("{}-{}_tabulated.dat", Faunus::atoms.at(i).name, Faunus::atoms.at(j).name);
            if (auto stream = std::ofstream(filename); stream) {
                stream_pair_potential(stream, i, j);
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
        save_potentials();
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
    for (double r = rmin + dr; r < rmax; r += dr) {
        double error = std::fabs(operator()(particle1, particle2, r *r, {r, 0, 0}) -
                                 FunctorPotential::operator()(particle1, particle2, r *r, {r, 0, 0}));
        max_error = std::max(error, max_error);
    }
    faunus_logger->debug(
        "{}-{} interaction splined between [{:6.2f}:{:6.2f}] {} using {} knots w. maximum absolute error of {:.1E} kT",
        Faunus::atoms[i].name, Faunus::atoms[j].name, rmin, rmax, u8::angstrom, knotdata.numKnots(), max_error);
}

void NewCoulombGalore::setSelfEnergy() {
    selfEnergy = [lB = bjerrum_length, self_energy = pot.selfEnergyFunctor](const Particle &p) {
        return lB * self_energy({p.charge * p.charge, 0.0});
    }; // expose self-energy as a functor in potential base class
}

NewCoulombGalore::NewCoulombGalore(const std::string &name) : PairPotentialBase(name) { setSelfEnergy(); }

Point NewCoulombGalore::force(const Particle &a, const Particle &b, double, const Point &r) const {
    return bjerrum_length * pot.ion_ion_force(a.charge, b.charge, r);
}

void NewCoulombGalore::from_json(const json &j) {
    using namespace ::CoulombGalore; // namespace for external CoulombGalore library
    double relative_dielectric_constant = j.at("epsr");
    bjerrum_length = pc::bjerrumLength(relative_dielectric_constant);
    pot.setTolerance(j.value("utol", 0.005 / bjerrum_length));
    std::string method = j.at("type");
    if (method == "yukawa") {
        if (json _j(j); _j.value("shift", false)) { // zero energy and force at cutoff
            faunus_logger->debug("energy and force shifted yukawa uses the 'poisson' scheme with C=1 and D=1");
            _j["type"] = "poisson";
            _j["C"] = 1;
            _j["D"] = 1;
            pot.spline<::CoulombGalore::Poisson>(_j);
        } else { // non-shifted yukawa equals `plain` with exponential screening
            if (_j.contains("cutoff")) {
                throw ConfigurationError("unexpected 'cutoff' for non-shifted yukawa which is always infinity");
            }
            _j["type"] = "plain";
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
        pot.spline<::CoulombGalore::Ewald>(j);
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

Multipole::Multipole(const std::string &name) : NewCoulombGalore(name) {
    isotropic = false; // this potential is angular dependent
    setSelfEnergy();
}

void Multipole::setSelfEnergy() {
    selfEnergy = [lB = bjerrum_length, self_energy = pot.selfEnergyFunctor](const Particle &p) {
        double mu_x_mu = 0;   // dipole-dipole product
        if (p.hasExtension()) // only access dipole if the particle has extended properties
            mu_x_mu = p.getExt().mulen * p.getExt().mulen;
        return lB * self_energy({p.charge * p.charge, mu_x_mu});
    }; // expose self-energy as a functor in potential base class
}

Point Multipole::force(const Faunus::Particle &a, const Faunus::Particle &b, double, const Faunus::Point &r) const {
    Point mua = a.getExt().mu * a.getExt().mulen;
    Point mub = b.getExt().mu * b.getExt().mulen;
    Point ionion = pot.ion_ion_force(a.charge, b.charge, r);
    Point iondip = pot.ion_dipole_force(a.charge, mub, r) + pot.ion_dipole_force(b.charge, mua, r);
    Point dipdip = pot.dipole_dipole_force(mua, mub, r);
    return bjerrum_length * (ionion + iondip + dipdip);
}

} // namespace Potential
} // namespace Faunus
