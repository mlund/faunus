#include "energy.h"
#include "penalty.h"
#include "potentials.h"
#include "externalpotential.h"

namespace Faunus {
namespace Energy {

void EwaldData::update(const Point &box) {
    L = box;
    int kcc = std::ceil(kc);
    check_k2_zero = 0.1 * std::pow(2 * pc::pi / L.maxCoeff(), 2);
    int kVectorsLength = (2 * kcc + 1) * (2 * kcc + 1) * (2 * kcc + 1) - 1;
    if (kVectorsLength == 0) {
        kVectors.resize(3, 1);
        Aks.resize(1);
        kVectors.col(0) = Point(1, 0, 0); // Just so it is not the zero-vector
        Aks[0] = 0;
        kVectorsInUse = 1;
        Qion.resize(1);
        Qdip.resize(1);
    } else {
        double kc2 = kc * kc;
        kVectors.resize(3, kVectorsLength);
        Aks.resize(kVectorsLength);
        kVectorsInUse = 0;
        kVectors.setZero();
        Aks.setZero();
        int startValue = 1 - int(ipbc);
        for (int kx = 0; kx <= kcc; kx++) {
            double dkx2 = double(kx * kx);
            for (int ky = -kcc * startValue; ky <= kcc; ky++) {
                double dky2 = double(ky * ky);
                for (int kz = -kcc * startValue; kz <= kcc; kz++) {
                    double factor = 1.0;
                    if (kx > 0) // optimization of PBC Ewald (and always the case for IPBC Ewald)
                        factor *= 2;
                    if (ky > 0 && ipbc) // only for IPBC Ewald
                        factor *= 2;
                    if (kz > 0 && ipbc) // only for IPBC Ewald
                        factor *= 2;
                    double dkz2 = double(kz * kz);
                    Point kv = 2 * pc::pi * Point(kx / L.x(), ky / L.y(), kz / L.z());
                    double k2 = kv.dot(kv) + kappa2; // last term is only for Yukawa-Ewald
                    if (k2 < check_k2_zero) // Check if k2 != 0
                        continue;
                    if (spherical_sum)
                        if ((dkx2 / kc2) + (dky2 / kc2) + (dkz2 / kc2) > 1)
                            continue;
                    kVectors.col(kVectorsInUse) = kv;
                    Aks[kVectorsInUse] = factor * std::exp(-k2 / (4 * alpha * alpha)) / k2;
                    kVectorsInUse++;
                }
            }
        }
        Qion.resize(kVectorsInUse);
        Qdip.resize(kVectorsInUse);
        Aks.conservativeResize(kVectorsInUse);
        kVectors.conservativeResize(3, kVectorsInUse);
    }
}

EwaldData::EwaldData(const json &j) {
    alpha = j.at("alpha"); // damping-parameter
    rc = j.at("cutoff");   // real space cut-off
    kc = j.at("ncutoff");  // reciprocal space cut-off
    ipbc = j.value("ipbc", false); // using PBC or IPBC?
    spherical_sum = j.value("spherical_sum", true); // Using spherical summation of k-vectors in reciprocal space?
    lB = pc::lB(j.at("epsr"));
    eps_surf = j.value("epss", 0.0); // dielectric constant of surrounding medium
    const_inf = (eps_surf < 1) ? 0 : 1; // if unphysical (<1) use epsr infinity for surrounding medium
    kappa = j.value("kappa", 0.0);
    kappa2 = kappa * kappa;
}

void to_json(json &j, const EwaldData &d) {
    j = {{"lB", d.lB},
         {"ipbc", d.ipbc},
         {"epss", d.eps_surf},
         {"alpha", d.alpha},
         {"cutoff", d.rc},
         {"ncutoff", d.kc},
         {"wavefunctions", d.kVectors.cols()},
         {"spherical_sum", d.spherical_sum},
         {"kappa", d.kappa}};
}

double Example2D::energy(Change &) {
    double s = 1 + std::sin(2 * pc::pi * i.x()) + std::cos(2 * pc::pi * i.y());
    if (i.x() >= -2.00 && i.x() <= -1.25)
        return 1 * s;
    if (i.x() >= -1.25 && i.x() <= -0.25)
        return 2 * s;
    if (i.x() >= -0.25 && i.x() <= 0.75)
        return 3 * s;
    if (i.x() >= 0.75 && i.x() <= 1.75)
        return 4 * s;
    if (i.x() >= 1.75 && i.x() <= 2.00)
        return 5 * s;
    return 1e10;
}
Example2D::Example2D(const json &, Space &spc) : i(spc.p.at(0).pos) { name = "Example2D"; }

double ContainerOverlap::energy(Change &change) {
    // if (spc.geo.type not_eq Geometry::CUBOID) // cuboid have PBC in all directions
    if (change) {
        // all groups have been updated
        if (change.dV or change.all) {
            for (auto &g : spc.groups) // loop over *all* groups in system
                for (auto &p : g)      // loop over *all* active particles in group
                    if (spc.geo.collision(p.pos))
                        return pc::infty;
            return 0;
        }

        // only a subset of groups have been updated
        for (auto &d : change.groups) {
            auto &g = spc.groups[d.index];
            // all atoms were updated
            if (d.all) {
                for (auto &p : g) // loop over *all* active particles in group
                    if (spc.geo.collision(p.pos))
                        return pc::infty;
            } else
                // only a subset of atoms were updated
                for (int i : d.atoms) // loop over specific atoms
                    if (i < g.size())
                        if (spc.geo.collision((g.begin() + i)->pos))
                            return pc::infty;
        }
    }
    return 0;
}

// ------------- Isobaric ---------------

Isobaric::Isobaric(const json &j, Space &spc) : spc(spc) {
    name = "isobaric";
    cite = "Frenkel & Smith 2nd Ed (Eq. 5.4.13)";
    P = j.value("P/mM", 0.0) * 1.0_mM;
    if (P < 1e-10) {
        P = j.value("P/Pa", 0.0) * 1.0_Pa;
        if (P < 1e-10)
            P = j.at("P/atm").get<double>() * 1.0_atm;
    }
}

double Isobaric::energy(Change &change) {
    if (change.dV || change.all || change.dN) {
        size_t N = 0;
        for (auto &g : spc.groups) {
            if (!g.empty()) {
                if (g.atomic)
                    N += g.size();
                else
                    N++;
            }
        }
        double V = spc.geo.getVolume();
        return P * V - (N+1)*std::log(V);
    } else
        return 0;
}
void Isobaric::to_json(json &j) const {
    j["P/atm"] = P / 1.0_atm;
    j["P/mM"] = P / 1.0_mM;
    j["P/Pa"] = P / 1.0_Pa;
    _roundjson(j, 5);
}

Constrain::Constrain(const json &j, Space &spc) {
    using namespace Faunus::ReactionCoordinate;
    name = "constrain";
    type = j.at("type").get<std::string>();
    rc = ReactionCoordinate::createReactionCoordinate({{type, j}}, spc);
}

double Constrain::energy(Change &change) {
    if (change) {
        double val = (*rc)();     // calculate reaction coordinate
        if (not rc->inRange(val)) // is it within allowed range?
            return pc::infty;     // if not, return infinite energy
    }
    return 0;
}
void Constrain::to_json(json &j) const {
    j = json(*rc).at(type);
    j.erase("resolution");
    j["type"] = type;
}
void Bonded::update_intra() {
    using namespace Potential;
    intra.clear();
    for (size_t i = 0; i < spc.groups.size(); i++) {
        auto &group = spc.groups.at(i);
        for (auto &bond : molecules.at(group.id).bonds) {
            intra[i].push_back<BondData>(bond->clone()); // deep copy BondData from MoleculeData
            intra[i].back()->shift(std::distance(spc.p.begin(), group.begin()));
            Potential::setBondEnergyFunction(intra[i].back(), spc.p);
        }
    }
}
double Bonded::sum_energy(const Bonded::BondVector &bonds) const {
    double energy = 0;
    for (auto &bond : bonds) {
        assert(bond->hasEnergyFunction());
        energy += bond->energy(spc.geo.getDistanceFunc());
    }
    return energy;
}
double Bonded::sum_energy(const Bonded::BondVector &bonds, const std::vector<int> &particles_ndx) const {
    double energy = 0;
    // outer loop over bonds to ensure that each bond is counted at most once
    for (auto &bond : bonds) {
        for (auto particle_ndx : particles_ndx) {
            if (std::find(bond->index.begin(), bond->index.end(), particle_ndx) != bond->index.end()) {
                assert(bond->hasEnergyFunction());
                energy += bond->energy(spc.geo.getDistanceFunc());
                break; // count each interaction at most once
            }
        }
    }
    return energy;
}
Bonded::Bonded(const json &j, Space &spc) : spc(spc) {
    name = "bonded";
    update_intra();
    if (j.is_object())
        if (j.count("bondlist") == 1)
            inter = j["bondlist"].get<BondVector>();
    for (auto &i : inter) // set all energy functions
        Potential::setBondEnergyFunction(i, spc.p);
}
void Bonded::to_json(json &j) const {
    if (!inter.empty())
        j["bondlist"] = inter;
    if (!intra.empty()) {
        json &_j = j["bondlist-intramolecular"];
        _j = json::array();
        for (auto &i : intra)
            for (auto &b : i.second)
                _j.push_back(b);
    }
}
double Bonded::energy(Change &change) {
    double energy = 0;
    if (change) {
        energy += sum_energy(inter); // energy of inter-molecular bonds

        if (change.all || change.dV) {              // compute all active groups
            for (auto &i : intra) {                 // energies of intra-molecular bonds
                if (!spc.groups[i.first].empty()) { // add only if group is active
                    energy += sum_energy(i.second);
                }
            }
        } else { // compute only the affected groups
            for (auto &group : change.groups) {
                auto &intra_group = intra[group.index];
                if (group.internal) {
                    if (group.all) { // all internal positions updated
                        if (not spc.groups[group.index].empty())
                            energy += sum_energy(intra_group);
                    } else { // only partial update of affected atoms
                        std::vector<int> atoms_ndx;
                        // an offset is the index of the first particle in the group
                        int offset = std::distance(spc.p.begin(), spc.groups[group.index].begin());
                        // add an offset to the group atom indices to get the absolute indices
                        std::transform(group.atoms.begin(), group.atoms.end(), std::back_inserter(atoms_ndx),
                                       [offset](int i) { return i + offset; });
                        energy += sum_energy(intra_group, atoms_ndx);
                    }
                }
            }
        } // for-loop over groups
    }
    return energy;
}

//---------- Hamiltonian ------------

void Hamiltonian::to_json(json &j) const {
    for (auto i : this->vec)
        j.push_back(*i);
}
void Hamiltonian::addEwald(const json &j, Space &spc) {
    // note this will currently not detect multipolar energies ore
    // deeply nested "coulomb" pair-potentials
    json _j;
    if (j.count("default") == 1) { // try to detect FunctorPotential
        for (auto &i : j["default"]) {
            if (i.count("coulomb") == 1) {
                _j = i["coulomb"];
                break;
            }
        }
    } else if (j.count("coulomb") == 1)
        _j = j["coulomb"];
    else
        return;

    if (_j.count("type"))
        if (_j.at("type") == "ewald")
            emplace_back<Energy::Ewald<>>(_j, spc);
}

Hamiltonian::Hamiltonian(Space &spc, const json &j) {
    using namespace Potential;

    typedef CombinedPairPotential<NewCoulombGalore, LennardJones> CoulombLJ; // temporary name
    typedef CombinedPairPotential<NewCoulombGalore, HardSphere> CoulombHS;
    typedef CombinedPairPotential<NewCoulombGalore, WeeksChandlerAndersen> CoulombWCA;
    typedef CombinedPairPotential<Coulomb, WeeksChandlerAndersen> PrimitiveModelWCA;
    typedef CombinedPairPotential<Coulomb, HardSphere> PrimitiveModel;

    if (not j.is_array())
        throw std::runtime_error("json array expected for energy");

    name = "hamiltonian";

    // add container overlap energy for non-cuboidal geometries
    if (spc.geo.type not_eq Geometry::CUBOID)
        emplace_back<Energy::ContainerOverlap>(spc);

    for (auto &m : j) { // loop over energy list
        size_t oldsize = vec.size();
        for (auto it : m.items()) {
            try {
                if (it.key() == "nonbonded_coulomblj")
                    emplace_back<Energy::Nonbonded<CoulombLJ>>(it.value(), spc, *this);

                else if (it.key() == "nonbonded_newcoulomblj")
                    emplace_back<Energy::Nonbonded<CoulombLJ>>(it.value(), spc, *this);

                else if (it.key() == "nonbonded_coulomblj_EM")
                    emplace_back<Energy::NonbondedCached<CoulombLJ>>(it.value(), spc, *this);

                else if (it.key() == "nonbonded_splined")
                    emplace_back<Energy::Nonbonded<TabulatedPotential>>(it.value(), spc, *this);

                else if (it.key() == "nonbonded" or it.key() == "nonbonded_exact")
                    emplace_back<Energy::Nonbonded<FunctorPotential>>(it.value(), spc, *this);

                else if (it.key() == "nonbonded_cached")
                    emplace_back<Energy::NonbondedCached<TabulatedPotential>>(it.value(), spc, *this);

                else if (it.key() == "nonbonded_coulombwca")
                    emplace_back<Energy::Nonbonded<CoulombWCA>>(it.value(), spc, *this);

                else if (it.key() == "nonbonded_pm" or it.key() == "nonbonded_coulombhs")
                    emplace_back<Energy::Nonbonded<PrimitiveModel>>(it.value(), spc, *this);

                else if (it.key() == "nonbonded_pmwca")
                    emplace_back<Energy::Nonbonded<PrimitiveModelWCA>>(it.value(), spc, *this);

                // this should be moved into `Nonbonded` and added when appropriate
                // Nonbonded now has access to Hamiltonian (*this) and can therefore
                // add energy terms
                addEwald(it.value(), spc); // add reciprocal Ewald terms if appropriate

                if (it.key() == "bonded")
                    emplace_back<Energy::Bonded>(it.value(), spc);

                else if (it.key() == "customexternal")
                    emplace_back<Energy::CustomExternal>(it.value(), spc);

                else if (it.key() == "akesson")
                    emplace_back<Energy::ExternalAkesson>(it.value(), spc);

                else if (it.key() == "confine")
                    emplace_back<Energy::Confine>(it.value(), spc);

                else if (it.key() == "constrain")
                    emplace_back<Energy::Constrain>(it.value(), spc);

                else if (it.key() == "example2d")
                    emplace_back<Energy::Example2D>(it.value(), spc);

                else if (it.key() == "isobaric")
                    emplace_back<Energy::Isobaric>(it.value(), spc);

                else if (it.key() == "penalty")
#ifdef ENABLE_MPI
                    emplace_back<Energy::PenaltyMPI>(it.value(), spc);
#else
                    emplace_back<Energy::Penalty>(it.value(), spc);
#endif
#if defined ENABLE_FREESASA
                else if (it.key() == "sasa")
                    emplace_back<Energy::SASAEnergy>(it.value(), spc);
#endif
                // additional energies go here...

                else if (it.key() == "maxenergy") {
                    maxenergy = it.value().get<double>();
                    continue;
                }

                if (vec.size() == oldsize)
                    throw std::runtime_error("unknown term");

            } catch (std::exception &e) {
                throw std::runtime_error("energy '" + it.key() + "': " + e.what() + usageTip[it.key()]);
            }
        } // end of loop over energy input terms
    }
    // Check if there are molecules with bonds and warn
    // if "bonded" has not been added
    for (auto &a : Faunus::molecules)
        if (not a.bonds.empty() and this->find<Energy::Bonded>().empty())
            faunus_logger->warn(a.name + " bonds specified in topology but missing in energy");
}
double Hamiltonian::energy(Change &change) {
    double du = 0;
    for (auto i : this->vec) { // loop over terms in Hamiltonian
        i->key = key;
        i->timer.start(); // time each term
        du += i->energy(change);
        i->timer.stop();
        if (du >= maxenergy)
            break; // stop summing energies
    }
    return du;
}
void Hamiltonian::init() {
    for (auto i : this->vec)
        i->init();
}
void Hamiltonian::sync(Energybase *basePtr, Change &change) {
    auto other = dynamic_cast<decltype(this)>(basePtr);
    if (other)
        if (other->size() == size()) {
            for (size_t i = 0; i < size(); i++)
                this->vec[i]->sync(other->vec[i].get(), change);
            return;
        }
    throw std::runtime_error("hamiltonian mismatch");
}

#ifdef ENABLE_FREESASA

SASAEnergy::SASAEnergy(Space &spc, double cosolute_concentration, double probe_radius)
    : spc(spc), cosolute_concentration(cosolute_concentration)
{
    name = "sasa"; // todo predecessor constructor
    cite = "doi:10.12688/f1000research.7931.1"; // todo predecessor constructor
    parameters = freesasa_default_parameters;
    parameters.probe_radius = probe_radius;
    init();
}

SASAEnergy::SASAEnergy(const json &j, Space &spc)
    : SASAEnergy(spc, j.value("molarity", 0.0) * 1.0_molar, j.value("radius", 1.4) * 1.0_angstrom) {}

void SASAEnergy::updatePositions([[gnu::unused]] const ParticleVector &p) {
    assert(p.size() == spc.positions().size());
    positions.clear();
    for(auto pos: spc.positions()) {
        auto xyz = pos.data();
        positions.insert(positions.end(), xyz, xyz+3);
    }
}

void SASAEnergy::updateRadii(const ParticleVector &p) {
    radii.resize(p.size());
    std::transform(p.begin(), p.end(), radii.begin(),
                   [](auto &a) { return atoms[a.id].sigma * 0.5; });
}

void SASAEnergy::updateSASA(const ParticleVector &p, const Change &) {
    updateRadii(p);
    updatePositions(p);
    auto result = freesasa_calc_coord(positions.data(), radii.data(), p.size(), &parameters);
    if(result) {
        assert(result->n_atoms == p.size());
        sasa.clear();
        sasa.insert(sasa.begin(), result->sasa, result->sasa + result->n_atoms); // copy
        assert(sasa.size() == p.size());
        freesasa_result_free(result);
    } else {
        throw std::runtime_error("FreeSASA failed");
    }
}

void SASAEnergy::init() {
    auto box = spc.geo.getLength();
    auto box_pbc = box;
    spc.geo.boundary(box_pbc);
    if(box_pbc != box) {
        faunus_logger->error("PBC applied, but PBC not implemented for FreeSASA. Expect unphysical results.");
    }
    Change change;
    change.all = true;
    updateSASA(spc.p, change);
}

double SASAEnergy::energy(Change &change) {
    double u = 0, A = 0;
    updateSASA(spc.p, change); // ideally we want
    for (size_t i = 0; i < spc.p.size(); ++i) {
        auto &a = atoms[spc.p[i].id];
        u += sasa[i] * (a.tension + cosolute_concentration * a.tfe);
        A += sasa[i];
    }
    avgArea += A; // sample average area for accepted confs.
    return u;
}

void SASAEnergy::sync(Energybase *basePtr, Change &c) {
    auto other = dynamic_cast<decltype(this)>(basePtr);
    if (other) {
        if (c.all || c.dV) {
            radii = other->radii;
            positions = other->positions;
            sasa = other->sasa;
        } else {
            for (auto &d : c.groups) {
                int offset = std::distance(spc.p.begin(), spc.groups.at(d.index).begin());
                for (int j : d.atoms) {
                    int i = j + offset;
                    radii[i] = other->radii[i];
                    sasa[i] = other->sasa[i];
                    for(size_t k = 0; k < 3; ++k) {
                        positions[3*i + k] = spc.positions()[i][k];
                    }
                }
            }
        }
    }
}

void SASAEnergy::to_json(json &j) const {
    using namespace u8;
    j["molarity"] = cosolute_concentration / 1.0_molar;
    j["radius"] = parameters.probe_radius / 1.0_angstrom;
    j[bracket("SASA") + "/" + angstrom + squared] = avgArea.avg() / 1.0_angstrom;
    _roundjson(j, 5); // set json output precision
}
#endif
} // end of namespace Energy
} // end of namespace Faunus
