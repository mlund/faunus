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
                (*matrix)(i.id(), j.id()) = modifier(extractor(i));
            } else {
                (*matrix)(i.id(), j.id()) = modifier(combinator(extractor(i), extractor(j)));
            }
        }
    }
    return matrix;
}

TPairMatrixPtr PairMixer::createPairMatrix(const std::vector<AtomData> &atoms,
                                           const std::vector<InteractionData> &interactions) {
    TPairMatrixPtr matrix = PairMixer::createPairMatrix(atoms);
    auto dimension = std::min(matrix->rows(), matrix->cols());
    for (auto &i : interactions) {
        if (i.atom_id[0] >= 0 && i.atom_id[1] >= 0 && i.atom_id[0] < dimension && i.atom_id[1] < dimension) {
            // interaction is always symmetric
            (*matrix)(i.atom_id[0], i.atom_id[1]) = (*matrix)(i.atom_id[1], i.atom_id[0]) =
                modifier(extractor(i.interaction));
        } else {
            throw std::runtime_error("atomtype index out of range");
        }
    }
    return matrix;
}

void from_json(const json &j, std::vector<InteractionData> &interactions) {
    auto &custom_list = j.at("custom");
    if (!custom_list.is_object() && !custom_list.is_array()) {
        throw PairPotentialException("custom parameters syntax error");
    }

    for (auto custom_pair = custom_list.begin(); custom_pair != custom_list.end(); ++custom_pair) {
        // atomdata is an array with items ecapsulated as objects hence we emulate here
        AtomData a = custom_list.is_object() ? json::object({{custom_pair.key(), *custom_pair}}) : (*custom_pair);
        auto atoms_name = words2vec<std::string>(a.name);
        if (atoms_name.size() == 2) {
            auto atom0 = findName(atoms, atoms_name[0]);
            auto atom1 = findName(atoms, atoms_name[1]);
            if (atom0 == atoms.end() or atom1 == atoms.end()) {
                throw PairPotentialException(("unknown atom(s): ["s + atoms_name[0] + " " + atoms_name[1] + "]"));
            }
            interactions.push_back({{atom0->id(), atom1->id()}, a});
        } else {
            throw PairPotentialException("custom parameters require exactly two space-separated atoms");
        }
    }
}

void to_json(json &j, const std::vector<InteractionData> &interactions) {
    if (!interactions.empty()) {
        j = json::array();
        for (auto &i : interactions) {
            j.push_back(i.interaction);
        }
    }
}

// =============== PairPotentialBase ===============

PairPotentialBase::PairPotentialBase(const std::string &name, const std::string &cite, bool isotropic)
    : name(name), cite(cite), isotropic(isotropic) {}

Point PairPotentialBase::force(const Particle &, const Particle &, double, const Point &) const {
    assert(false && "We should never reach this point!");
    return {0, 0, 0};
}

// =============== MixerPairPotentialBase ===============

void MixerPairPotentialBase::init() {
    json j_combination_rule = combination_rule;
    faunus_logger->debug("Combination rule {} in effect for the {} potential.", j_combination_rule, name);
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
            *custom_pairs = j;
        }
    } catch (const PairPotentialException &e) {
        faunus_logger->error(std::string(e.what()) + " in potential " + name);
        throw std::runtime_error("error deserialising potential " + name + " from json");
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
    j["epsr"] = pc::lB2epsr(lB);
    j["lB"] = lB;
}

void Coulomb::from_json(const json &j) { lB = pc::lB(j.at("epsr")); }

void DipoleDipole::to_json(json &j) const {
    j["epsr"] = pc::lB2epsr(lB);
    j["lB"] = lB;
}

void DipoleDipole::from_json(const json &j) { lB = pc::lB(j.at("epsr")); }

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
    assertKeys(j, {"shift", "molarity", "radius"});
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
    extract_sigma = [sigma_name](const AtomData &a) -> double { return a.getProperty(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const AtomData &a) -> double { return a.getProperty(epsilon_name) * 1.0_kJmol; };
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
    extract_sigma = [sigma_name](const AtomData &a) -> double { return a.getProperty(sigma_name) * 1.0_angstrom; };
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
    extract_sigma = [sigma_name](const AtomData &a) -> double { return a.getProperty(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const AtomData &a) -> double { return a.getProperty(epsilon_name) * 1.0_kJmol; };
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
    extract_sigma = [sigma_name](const AtomData &a) -> double { return a.getProperty(sigma_name) * 1.0_angstrom; };
    auto epsilon_name = j.value("eps", "eps");
    json_extra_params["eps"] = epsilon_name;
    extract_epsilon = [epsilon_name](const AtomData &a) -> double { return a.getProperty(epsilon_name) * 1.0_kJmol; };
}

// =============== Polarizability ===============

void Polarizability::from_json(const json &j) {
    epsr = j.at("epsr").get<double>();
    double lB = pc::lB(epsr);
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
        self_energy_vector.push_back(pot->selfEnergy);
        assert(self_energy_vector.back());
    } else
        faunus_logger->trace("Failed to register non-defined selfEnergy() for {}", pot->name);
}

FunctorPotential::uFunc FunctorPotential::combineFunc(json &j) {
    uFunc u = [](const Particle &, const Particle &, const Point &) { return 0.0; };
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
                            _u = std::get<0>(potlist) = it.value();
                            std::get<0>(potlist).to_json(it.value());
                            if (not have_monopole_self_energy) {
                                registerSelfEnergy(&std::get<0>(potlist));
                                have_monopole_self_energy = true;
                                assert(not self_energy_vector.empty());
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
                            faunus_logger->error("'dipoledipole' is deprecated, use 'multipole' instead");
                        } else if (it.key() == "stockmayer") {
                            faunus_logger->error("'stockmayer' is deprecated, use 'lennardjones'+'multipole' instead");
                        } else if (it.key() == "multipole") {
                            _u = std::get<12>(potlist) = it.value();
                            std::get<12>(potlist).to_json(it.value()); // store json output
                            isotropic = false;                         // potential is now angular dependent
                            if (not have_dipole_self_energy) {
                                registerSelfEnergy(&std::get<12>(potlist));
                                have_dipole_self_energy = true;
                                assert(not self_energy_vector.empty());
                            }
                        }
                        // place additional potentials here...
                    } catch (std::exception &e) {
                        throw std::runtime_error(it.key() + ": " + e.what() + usageTip[it.key()]);
                    }

                    if (_u != nullptr) // if found, sum them into new function object
                        u = [u, _u](const Particle &a, const Particle &b, const Point &r) {
                            return u(a, b, r) + _u(a, b, r);
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
    _j = j;
    umatrix = decltype(umatrix)(atoms.size(), combineFunc(_j.at("default")));
    for (auto it = _j.begin(); it != _j.end(); ++it) {
        auto atompair = words2vec<std::string>(it.key()); // is this for a pair of atoms?
        if (atompair.size() == 2) {
            auto ids = names2ids(atoms, atompair);
            umatrix.set(ids[0], ids[1], combineFunc(it.value()));
        }
    }
    // at the moment we assume that self energies apply to all
    // pairs of particle types
    if (self_energy_vector.empty())
        selfEnergy = nullptr;
    else {
        // sum self energies from all registered pair-potentials
        selfEnergy = [&](const Particle &p) {
            double sum = 0;
            for (auto &func : self_energy_vector) {
                assert(func);
                sum += func(p);
            }
            return sum;
        };
        faunus_logger->debug("added {} selfEnergy terms to {}", self_energy_vector.size(), name);
    }
}

//---------------- TabulatedPotential ---------------------

void TabulatedPotential::from_json(const json &j) {
    FunctorPotential::from_json(j);

    // if user specifies an anisotropic potential, make sure to bail out
    if (not isotropic)
        throw std::runtime_error("cannot spline anisotropic potentials");

    tblt.setTolerance(j.value("utol", 1e-5), j.value("ftol", 1e-2));
    double u_at_rmin = j.value("u_at_rmin", 20);
    double u_at_rmax = j.value("u_at_rmax", 1e-6);
    hardsphere = j.value("hardsphere", false);

    // build matrix of spline data, each element corresponding
    // to a pair of atom types
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t k = 0; k <= i; ++k) {
            if (atoms[i].implicit == false and atoms[k].implicit == false) {
                Particle a = atoms.at(i);
                Particle b = atoms.at(k);
                double rmin2 = .5 * (atoms[i].sigma + atoms[k].sigma);
                rmin2 = rmin2 * rmin2;
                double rmax2 = rmin2 * 100;
                auto it = j.find("cutoff_g2g");
                if (j.count("rmax") == 1) {
                    rmax2 = std::pow(j.at("rmax").get<double>(), 2);
                } else if (it != j.end()) {
                    if (it->is_number())
                        rmax2 = std::pow(it->get<double>(), 2);
                    else if (it->is_object())
                        rmax2 = std::pow(it->at("default").get<double>(), 2);
                }

                // adjust lower splining distance to match
                // the given energy threshold (u_at_min2)
                double dr = 1e-2;
                while (rmin2 >= dr) {
                    double u = std::fabs(this->umatrix(i, k)(a, b, {0, 0, sqrt(rmin2)}));
                    if (u > u_at_rmin * 1.1)
                        rmin2 = rmin2 + dr;
                    else if (u < u_at_rmin / 1.1)
                        rmin2 = rmin2 - dr;
                    else
                        break;
                }

                assert(rmin2 >= 0);

                while (rmax2 >= dr) {
                    double u = std::fabs(this->umatrix(i, k)(a, b, {0, 0, sqrt(rmax2)}));
                    if (u > u_at_rmax)
                        rmax2 = rmax2 + dr;
                    else
                        break;
                }

                assert(rmin2 < rmax2);

                Ttable knotdata = tblt.generate(
                    [&](double r2) {
                        return this->umatrix(i, k)(a, b, {0, 0, sqrt(r2)});
                    },
                    rmin2, rmax2);

                // assert if potential is negative for r<rmin
                if (tblt.eval(knotdata, knotdata.rmin2 + dr) < 0)
                    knotdata.isNegativeBelowRmin = true;

                tmatrix.set(i, k, knotdata);
                if (j.value("to_disk", false)) {
                    std::ofstream f(atoms[i].name + "-" + atoms[k].name + "_tabulated.dat"); // output file
                    f << "# r splined exact\n";
                    Point r = {dr, 0, 0}; // variable distance vector between particle a and b
                    for (; r.x() < sqrt(rmax2); r.x() += dr)
                        f << r.x() << " " << operator()(a, b, r) << " " << this->umatrix(i, k)(a, b, r) << "\n";
                }
            }
        }
    }
}

NewCoulombGalore::NewCoulombGalore(const std::string &name) : PairPotentialBase(name) {
    selfEnergy = [&](const Particle &p) {
        return pot.self_energy({p.charge * p.charge, 0});
    }; // expose self-energy as a functor in potential base class
}

Point NewCoulombGalore::force(const Particle &a, const Particle &b, double, const Point &r) const {
    return lB * pot.ion_ion_force(a.charge, b.charge, r);
}

void NewCoulombGalore::from_json(const json &j) {
    using namespace ::CoulombGalore; // namespace for external CoulombGalore library
    double epsr = j.at("epsr");
    lB = pc::lB(epsr); // Bjerrum length
    std::string type = j.at("type");
    if (type == "yukawa")
        faunus_logger->error("'yukawa' is deprecated, use 'plain' with 'debyelength'");
    if (type == "plain")
        pot.spline<::CoulombGalore::Plain>(j);
    else if (type == "qpotential")
        pot.spline<::CoulombGalore::qPotential>(j);
    else if (type == "wolf")
        pot.spline<Wolf>(j);
    else if (type == "poisson")
        pot.spline<Poisson>(j);
    else if (type == "fanourgakis")
        pot.spline<Fanourgakis>(j);
    else if (type == "ewald")
        pot.spline<Ewald>(j);
    else
        throw std::runtime_error("unknown coulomb scheme");
}

void NewCoulombGalore::to_json(json &j) const {
    pot.to_json(j);
    j["lB"] = lB;
}

Multipole::Multipole(const std::string &name) : NewCoulombGalore(name) {
    isotropic = false; // this potential is angular dependent
    selfEnergy = [&](const Particle &p) {
        double mu_x_mu = 0;   // dipole-dipole product
        if (p.hasExtension()) // only access dipole of the particle has extended properties
            mu_x_mu = p.getExt().mulen * p.getExt().mulen;
        return pot.self_energy({p.charge * p.charge, mu_x_mu});
    }; // expose CoulombGalore self-energy as a functor in potential base class
}

Point Multipole::force(const Faunus::Particle &a, const Faunus::Particle &b, double, const Faunus::Point &r) const {
    Point mua = a.getExt().mu * a.getExt().mulen;
    Point mub = b.getExt().mu * b.getExt().mulen;
    Point ionion = pot.ion_ion_force(a.charge, b.charge, r);
    Point iondip = pot.ion_dipole_force(a.charge, mub, r) + pot.ion_dipole_force(b.charge, mua, r);
    Point dipdip = pot.dipole_dipole_force(mua, mub, r);
    return lB * (ionion + iondip + dipdip);
}

} // namespace Potential
} // namespace Faunus
