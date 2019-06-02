
#include "mpi.h"
#include "energy.h"
#include "penalty.h"
#include "potentials.h"
#include "externalpotential.h"

namespace Faunus {
namespace Energy {

void Energybase::to_json(json &) const {}

void Energybase::sync(Energybase *, Change &) {}

void Energybase::init() {}

void to_json(json &j, const Energybase &base) {
    assert(not base.name.empty());
    if (base.timer)
        j["relative time"] = base.timer.result();
    if (not base.cite.empty())
        j[base.name]["reference"] = base.cite;
    base.to_json(j[base.name]);
}

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
                    if (kx > 0)
                        factor *= 2;
                    if (ky > 0 && ipbc)
                        factor *= 2;
                    if (kz > 0 && ipbc)
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

void from_json(const json &j, EwaldData &d) {
    d.alpha = j.at("alpha");
    d.rc = j.at("cutoff");
    d.kc = j.at("kcutoff");
    d.ipbc = j.value("ipbc", false);
    d.spherical_sum = j.value("spherical_sum", true);
    d.lB = pc::lB(j.at("epsr"));
    d.eps_surf = j.value("epss", 0.0);
    d.const_inf = (d.eps_surf < 1) ? 0 : 1; // if unphysical (<1) use epsr infinity for surrounding medium
    d.kappa = j.value("kappa", 0.0);
    d.kappa2 = d.kappa*d.kappa;
}

void to_json(json &j, const EwaldData &d) {
    j = {{"lB", d.lB},
         {"ipbc", d.ipbc},
         {"epss", d.eps_surf},
         {"alpha", d.alpha},
         {"cutoff", d.rc},
         {"kcutoff", d.kc},
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
                    if (spc.geo.collision((g.begin() + i)->pos))
                        return pc::infty;
        }
    }
    return 0;
}

// Remove?
/*
SelfEnergy::SelfEnergy(const json &j, Space &spc) : spc(spc) {
    name = "selfenergy";
    type = j.at("type");
    rc = j.at("cutoff");
    epsr = j.at("epsr");
    lB = pc::lB(epsr);

    selfenergy_ion_prefactor = 0.0;
    selfenergy_dipole_prefactor = 0.0;
    if (type == "reactionfield") {
        epsrf = j.at("epsrf");
        selfenergy_ion_prefactor = 1.5 * epsrf / (2.0 * epsrf + epsr); // Correct?!, see Eq.14 in DOI: 10.1021/jp510612w
        selfenergy_dipole_prefactor = 2.0*(epsr - epsrf)/(2.0*epsrf + epsr); // Preliminary, needs to be checked!
    }
    if (type == "fanourgakis") {
        selfenergy_ion_prefactor = -0.875;
        selfenergy_dipole_prefactor = 0.0;
    }
    if (type == "poisson") {
        C = j.at("C");
        D = j.at("D");
        selfenergy_ion_prefactor = -double(C+D)/double(C);
        selfenergy_dipole_prefactor = 0.0; // check this!
    }
    if (type == "yukawapoisson") {
        C = j.at("C");
        D = j.at("D");
        kappa = j.at("kappa");
        selfenergy_ion_prefactor = -double(C+D)/double(C);
        selfenergy_dipole_prefactor = 0.0; // check this!
    }
    if (type == "yukawa") {
        kappa = j.at("kappa");
        selfenergy_ion_prefactor = 0.0; // check this!
        selfenergy_dipole_prefactor = 0.0; // check this!
    }
    if (type == "qpotential" || type == "q2potential") {
        selfenergy_ion_prefactor = -1.0;
        selfenergy_dipole_prefactor = -1.0;
    }
    if (type == "fennell") {
        alpha = j.at("alpha");
        selfenergy_ion_prefactor = -(erfc(alpha*rc) + alpha*rc / sqrt(pc::pi) * (1.0 + exp(-alpha*alpha*rc*rc)));
        selfenergy_dipole_prefactor = -0.5*( erfc(alpha*rc) + 2.0*alpha*rc/sqrt(pc::pi)*exp(-alpha*alpha*rc*rc) +
(4.0/3.0)*pow(alpha*rc,3.0)/sqrt(pc::pi) );
    }
    if (type == "wolf") {
        alpha = j.at("alpha");
        selfenergy_ion_prefactor = -0.5*(erfc(alpha*rc) + 2.0*alpha*rc / sqrt(pc::pi));
        selfenergy_dipole_prefactor = -0.5*( erfc(alpha*rc) + 2.0*alpha*rc/sqrt(pc::pi)*exp(-alpha*alpha*rc*rc) +
(4.0/3.0)*pow(alpha*rc,3.0)/sqrt(pc::pi) );
    }
    if (type == "ewald") {
        alpha = j.at("alpha");
        selfenergy_ion_prefactor = -alpha*rc/std::sqrt(pc::pi);
        selfenergy_dipole_prefactor = -2.0*pow(alpha*rc,3.0)/3.0/std::sqrt(pc::pi);
    }
}

double SelfEnergy::energy(Change &change) {
    double Eq = 0;
    double Emu = 0;
    if (change.dN)
        for (auto cg : change.groups) {
            auto g = spc.groups.at(cg.index);
            for (auto i : cg.atoms)
                if (i < g.size()) {
                    Eq += std::pow((g.begin() + i)->charge, 2);
                    Emu += std::pow((g.begin() + i)->getExt().mulen, 2);
                }
        }
    else if (change.all and not change.dV)
        for (auto g : spc.groups)
            for (auto i : g) {
                Eq += i.charge * i.charge;
                Emu += i.getExt().mulen * i.getExt().mulen;
            }
    return ( selfenergy_ion_prefactor * Eq / rc + selfenergy_dipole_prefactor*Emu/pow(rc,3.0) )*lB;
}*/

// ------------- ParticleSelfEnergy ---------------

ParticleSelfEnergy::ParticleSelfEnergy(Space &spc, std::function<double(Particle &)> selfEnergy)
    : spc(spc), selfEnergy(selfEnergy) {
    name = "selfenergy";
#ifndef NDEBUG
    // test if self energy can be called
    Particle myparticle;
    if (this->selfEnergy)
        this->selfEnergy(myparticle);
#endif
}

double ParticleSelfEnergy::energy(Change &change) {
    double u = 0;
    if (selfEnergy) {
        if (change.dN)
            for (auto &cg : change.groups) {
                auto &g = spc.groups.at(cg.index);
                for (auto i : cg.atoms)
                    if (i < g.size())
                        u += selfEnergy(*(g.begin() + i));
            }
        else if (change.all and not change.dV)
            for (auto &g : spc.groups)
                for (auto &i : g)
                    u += selfEnergy(i);
    }
    return u;
}

void ParticleSelfEnergy::sync(Energybase *basePtr, Change &change) {
    if (change.all) {
        auto other = dynamic_cast<decltype(this)>(basePtr);
        assert(other);
        selfEnergy = other->selfEnergy;
    }
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
        double V = spc.geo.getVolume();
        size_t N = 0;
        for (auto &g : spc.groups)
            if (!g.empty()) {
                if (g.atomic)
                    N += g.size();
                else
                    N++;
            }
        return P * V - (N + 1) * std::log(V);
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
    try {
        if (type == "atom")
            rc = std::make_shared<AtomProperty>(j, spc);
        else if (type == "molecule")
            rc = std::make_shared<MoleculeProperty>(j, spc);
        else if (type == "system")
            rc = std::make_shared<SystemProperty>(j, spc);
        else if (type == "cmcm")
            rc = std::make_shared<MassCenterSeparation>(j, spc);
        if (rc == nullptr)
            throw std::runtime_error("unknown coordinate type");

    } catch (std::exception &e) {
        throw std::runtime_error("error for reaction coordinate '" + type + "': " + e.what() +
                                 usageTip["coords=[" + type + "]"]);
    }
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
    j = *rc;
    j["type"] = type;
    j.erase("resolution");
}
void Bonded::update_intra() {
    using namespace Potential;
    intra.clear();
    for (size_t i = 0; i < spc.groups.size(); i++) {
        auto &group = spc.groups.at(i);
        for (auto &bond : molecules.at(group.id).bonds) {
            intra[i].push_back(bond->clone()); // deep copy BondData from MoleculeData
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
    // note this will not find deeper placed coulomb potentials
    // in FunctorPotential etc. Nor dipolar energies
    if (j.count("coulomb") == 1)
        if (j["coulomb"].count("type") == 1)
            if (j["coulomb"].at("type") == "ewald")
                push_back<Energy::Ewald<>>(j["coulomb"], spc);
}

Hamiltonian::Hamiltonian(Space &spc, const json &j) {
    using namespace Potential;

    typedef CombinedPairPotential<CoulombGalore, LennardJones> CoulombLJ;
    typedef CombinedPairPotential<CoulombGalore, HardSphere> CoulombHS;
    typedef CombinedPairPotential<CoulombGalore, WeeksChandlerAndersen> CoulombWCA;
    typedef CombinedPairPotential<Coulomb, WeeksChandlerAndersen> PrimitiveModelWCA;
    typedef CombinedPairPotential<Coulomb, HardSphere> PrimitiveModel;

    name = "hamiltonian";

    // add container overlap energy for non-cuboidal geometries
    if (spc.geo.type not_eq Geometry::CUBOID)
        push_back<Energy::ContainerOverlap>(spc);

    for (auto &m : j.at("energy")) { // loop over move list
        size_t oldsize = vec.size();
        for (auto it = m.begin(); it != m.end(); ++it) {
            try {
                if (it.key() == "nonbonded_coulomblj")
                    push_back<Energy::Nonbonded<CoulombLJ>>(it.value(), spc, *this);

                if (it.key() == "nonbonded_coulomblj_EM")
                    push_back<Energy::NonbondedCached<CoulombLJ>>(it.value(), spc, *this);

                if (it.key() == "nonbonded")
                    push_back<Energy::Nonbonded<TabulatedPotential>>(it.value(), spc, *this);

                if (it.key() == "nonbonded_exact")
                    push_back<Energy::Nonbonded<FunctorPotential>>(it.value(), spc, *this);

                if (it.key() == "nonbonded_cached")
                    push_back<Energy::NonbondedCached<TabulatedPotential>>(it.value(), spc, *this);

                if (it.key() == "nonbonded_coulombwca")
                    push_back<Energy::Nonbonded<CoulombWCA>>(it.value(), spc, *this);

                if (it.key() == "nonbonded_pm" or it.key() == "nonbonded_coulombhs")
                    push_back<Energy::Nonbonded<PrimitiveModel>>(it.value(), spc, *this);

                if (it.key() == "nonbonded_pmwca")
                    push_back<Energy::Nonbonded<PrimitiveModelWCA>>(it.value(), spc, *this);

                // this should be moved into `Nonbonded` and added when appropriate
                // Nonbonded now has access to Hamiltonian (*this) and can therefore
                // add energy terms
                addEwald(it.value(), spc); // add reciprocal Ewald terms if appropriate

                if (it.key() == "bonded")
                    push_back<Energy::Bonded>(it.value(), spc);

                if (it.key() == "customexternal")
                    push_back<Energy::CustomExternal>(it.value(), spc);

                if (it.key() == "akesson")
                    push_back<Energy::ExternalAkesson>(it.value(), spc);

                if (it.key() == "confine")
                    push_back<Energy::Confine>(it.value(), spc);

                if (it.key() == "constrain")
                    push_back<Energy::Constrain>(it.value(), spc);

                if (it.key() == "example2d")
                    push_back<Energy::Example2D>(it.value(), spc);

                if (it.key() == "isobaric")
                    push_back<Energy::Isobaric>(it.value(), spc);

                if (it.key() == "penalty")
#ifdef ENABLE_MPI
                    push_back<Energy::PenaltyMPI>(it.value(), spc);
#else
                    push_back<Energy::Penalty>(it.value(), spc);
#endif
#ifdef ENABLE_POWERSASA
                if (it.key() == "sasa")
                    push_back<Energy::SASAEnergy>(it.value(), spc);
#endif
                // additional energies go here...

                if (it.key() == "maxenergy") {
                    maxenergy = it.value().get<double>();
                    continue;
                }

                if (vec.size() == oldsize)
                    throw std::runtime_error("unknown term");

            } catch (std::exception &e) {
                throw std::runtime_error("Error adding energy '" + it.key() + "': " + e.what() + usageTip[it.key()]);
            }
        } // end of loop over energy input terms

        // if (have_nonbonded) {
        //    push_back<Energy::ParticleSelfEnergy>(spc);
        //}
    }
}
double Hamiltonian::energy(Change &change) {
    double du = 0;
    for (auto i : this->vec) {
        i->key = key;
        i->timer.start();
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
void SASAEnergy::updateSASA(const ParticleVector &p) {
    assert(ps != nullptr);
    radii.resize(p.size());
    std::transform(p.begin(), p.end(), radii.begin(),
                   [this](auto &a) { return atoms[a.id].sigma * 0.5 + this->probe; });

    ps->update_coords(spc.positions(), radii); // slowest step!

    for (size_t i = 0; i < p.size(); i++) {
        auto &a = atoms[p[i].id];
        if (std::fabs(a.tfe) > 1e-9 || std::fabs(a.tension) > 1e-9)
            ps->calc_sasa_single(i);
    }
    sasa = ps->getSasa();
    assert(sasa.size() == p.size());
}
void SASAEnergy::to_json(json &j) const {
    using namespace u8;
    j["molarity"] = conc / 1.0_molar;
    j["radius"] = probe / 1.0_angstrom;
    j[bracket("SASA") + "/" + angstrom + squared] = avgArea.avg() / 1.0_angstrom;
    _roundjson(j, 5);
}
void SASAEnergy::sync(Energybase *basePtr, Change &c) {
    auto other = dynamic_cast<decltype(this)>(basePtr);
    if (other) {
        if (c.all || c.dV) {
            radii = other->radii;
            sasa = other->sasa;
        } else {
            for (auto &d : c.groups) {
                int offset = std::distance(spc.p.begin(), spc.groups.at(d.index).begin());
                for (int j : d.atoms) {
                    int i = j + offset;
                    radii[i] = other->radii[i];
                    sasa[i] = other->sasa[i];
                }
            }
        }
    }
}
SASAEnergy::SASAEnergy(const json &j, Space &spc) : spc(spc) {
    name = "sasa";
    cite = "doi:10.1002/jcc.21844";
    probe = j.value("radius", 1.4) * 1.0_angstrom;
    conc = j.at("molarity").get<double>() * 1.0_molar;
    init();
}
void SASAEnergy::init() {
    radii.resize(spc.p.size());
    std::transform(spc.p.begin(), spc.p.end(), radii.begin(),
                   [this](auto &a) { return atoms[a.id].sigma * 0.5 + this->probe; });

    if (ps == nullptr)
        ps = std::make_shared<POWERSASA::PowerSasa<float, Point>>(spc.positions(), radii);
    updateSASA(spc.p);
}
double SASAEnergy::energy(Change &) {
    double u = 0, A = 0;
    /*
     * ideally we want to call `update` only if `key==NEW` but
     * syncronising the PowerSasa object is difficult since it's
     * non-copyable.
     */
    updateSASA(spc.p); // ideally we want
    for (size_t i = 0; i < spc.p.size(); ++i) {
        auto &a = atoms[spc.p[i].id];
        u += sasa[i] * (a.tension + conc * a.tfe);
        A += sasa[i];
    }
    avgArea += A; // sample average area for accepted confs. only
    return u;
}
} // end of namespace Energy
} // end of namespace Faunus
