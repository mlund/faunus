#pragma once

#include "space.h"
#include "bonds.h"
#include "auxiliary.h"
#include <range/v3/view.hpp>
#include <Eigen/Dense>

#ifdef ENABLE_POWERSASA
#include <power_sasa.h>
#endif

namespace Faunus {

namespace ReactionCoordinate {
struct ReactionCoordinateBase;
}

namespace Potential {
struct PairPotentialBase;
}

namespace Energy {

class Energybase {
  public:
    enum keys { OLD, NEW, NONE };
    keys key = NONE;
    std::string name;
    std::string cite;
    TimeRelativeOfTotal<std::chrono::microseconds> timer;
    virtual double energy(Change &) = 0; //!< energy due to change
    virtual void to_json(json &j) const; //!< json output
    virtual void sync(Energybase *, Change &);
    virtual void init();                               //!< reset and initialize
    virtual inline void force(std::vector<Point> &){}; // update forces on all particles
    inline virtual ~Energybase(){};
};

void to_json(json &j, const Energybase &base); //!< Converts any energy class to json object

/**
 * @brief Check for overlap between atoms and the simulation container
 *
 * If found infinite energy is returned. Not needed for cuboidal geometry
 * as there's nover any overlap due to PBC.
 */
struct ContainerOverlap : public Energybase {
    const Space &spc;
    ContainerOverlap(const Space &spc) : spc(spc) { name = "ContainerOverlap"; }
    double energy(Change &change) override;
};

/**
 * This holds Ewald setup and must *not* depend on particle type, nor depend on Space
 */
struct EwaldData {
    typedef std::complex<double> Tcomplex;
    Eigen::Matrix3Xd kVectors;   // k-vectors, 3xK
    Eigen::VectorXd Aks;         // 1xK, to minimize computational effort (Eq.24,DOI:10.1063/1.481216)
    Eigen::VectorXcd Qion, Qdip; // 1xK
    double alpha, rc, kc, check_k2_zero, lB;
    double const_inf, eps_surf, kappa, kappa2;
    bool spherical_sum = true;
    bool ipbc = false;
    int kVectorsInUse = 0;
    Point L; //!< Box dimensions

    void update(const Point &box);
};

void from_json(const json &j, EwaldData &d);

void to_json(json &j, const EwaldData &d);

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] Ewald - EwaldData") {
    using doctest::Approx;

    EwaldData data = R"({
                "ipbc": false, "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "kcutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json;

    data.update(Point(10, 10, 10));

    CHECK(data.ipbc == false);
    CHECK(data.const_inf == 1);
    CHECK(data.alpha == 0.894427190999916);
    CHECK(data.kVectors.cols() == 2975);
    CHECK(data.Qion.size() == data.kVectors.cols());

    data.ipbc = true;
    data.update(Point(10, 10, 10));
    CHECK(data.kVectors.cols() == 846);
    CHECK(data.Qion.size() == data.kVectors.cols());
}
#endif

/** @brief recipe or policies for ion-ion ewald */
template <bool eigenopt = false /** use Eigen matrix ops where possible */> struct PolicyIonIon {
    typedef typename ParticleVector::iterator iter;
    Space *spc;
    Space *old = nullptr; // set only if key==NEW at first call to `sync()`

    PolicyIonIon(Space &spc) : spc(&spc) {}

    void updateComplex(EwaldData &data) const {
        auto active = spc->activeParticles();
        if (eigenopt)
            if (data.ipbc == false) {
                auto pos = asEigenMatrix(active.begin().base(), active.end().base(), &Space::Tparticle::pos); //  Nx3
                auto charge =
                    asEigenVector(active.begin().base(), active.end().base(), &Space::Tparticle::charge); // Nx1
                Eigen::MatrixXd kr = pos.matrix() * data.kVectors; // Nx3 * 3xK = NxK
                data.Qion.real() = (kr.array().cos().colwise() * charge).colwise().sum();
                data.Qion.imag() = kr.array().sin().colwise().sum();
                return;
            }
        for (int k = 0; k < data.kVectors.cols(); k++) {
            const Point &kv = data.kVectors.col(k);
            EwaldData::Tcomplex Q(0, 0);
            if (data.ipbc)
                for (auto &i : active)
                    Q += kv.cwiseProduct(i.pos).array().cos().prod() * i.charge;
            else
                for (auto &i : active) {
                    double dot = kv.dot(i.pos);
                    Q += i.charge * EwaldData::Tcomplex(std::cos(dot), std::sin(dot));
                }
            data.Qion[k] = Q;
        }
    } //!< Update all k vectors

    void updateComplex(EwaldData &data, Change &change) const {
        assert(old != nullptr);
        assert(spc->p.size() == old->p.size());
        for (int k = 0; k < data.kVectors.cols(); k++) {
            auto &Q = data.Qion[k];
            Point q = data.kVectors.col(k);
            if (data.ipbc)
                for (auto cg : change.groups) {
                    auto g_new = spc->groups.at(cg.index);
                    auto g_old = old->groups.at(cg.index);
                    for (auto i : cg.atoms) {
                        if (i < g_new.size())
                            Q += q.cwiseProduct((g_new.begin() + i)->pos).array().cos().prod() *
                                 (g_new.begin() + i)->charge;
                        if (i < g_old.size())
                            Q -= q.cwiseProduct((g_old.begin() + i)->pos).array().cos().prod() *
                                 (g_old.begin() + i)->charge;
                    }
                }
            else
                for (auto cg : change.groups) {
                    auto g_new = spc->groups.at(cg.index);
                    auto g_old = old->groups.at(cg.index);
                    for (auto i : cg.atoms) {
                        if (i < g_new.size()) {
                            double _new = q.dot((g_new.begin() + i)->pos);
                            Q += (g_new.begin() + i)->charge * EwaldData::Tcomplex(std::cos(_new), std::sin(_new));
                        }
                        if (i < g_old.size()) {
                            double _old = q.dot((g_old.begin() + i)->pos);
                            Q -= (g_old.begin() + i)->charge * EwaldData::Tcomplex(std::cos(_old), std::sin(_old));
                        }
                    }
                }
        }
    } //!< Optimized update of k subset. Require access to old positions through `old` pointer

    double selfEnergy(const EwaldData &d, Change &change) {
        double Eq = 0;
        if (change.dN)
            for (auto cg : change.groups) {
                auto g = spc->groups.at(cg.index);
                for (auto i : cg.atoms)
                    if (i < g.size())
                        Eq += std::pow((g.begin() + i)->charge, 2);
            }
        else if (change.all and not change.dV)
            for (auto g : spc->groups)
                for (auto i : g)
                    Eq += i.charge * i.charge;
        return -d.alpha * Eq / std::sqrt(pc::pi) * d.lB;
    }

    double surfaceEnergy(const EwaldData &d, Change &change) {
        if (d.const_inf < 0.5)
            return 0;
        Point qr(0, 0, 0);
        if (change.all or change.dV)
            for (auto g : spc->groups)
                for (auto i : g)
                    qr += i.charge * i.pos;
        else if (change.groups.size() > 0)
            for (auto cg : change.groups) {
                auto g = spc->groups.at(cg.index);
                for (auto i : cg.atoms)
                    if (i < g.size())
                        qr += (g.begin() + i)->charge * (g.begin() + i)->pos;
            }
        return d.const_inf * 2 * pc::pi / ((2 * d.eps_surf + 1) * spc->geo.getVolume()) * qr.dot(qr) * d.lB;
    }

    double reciprocalEnergy(const EwaldData &d) {
        double E = 0;
        if (eigenopt) // known at compile time
            E = d.Aks.cwiseProduct(d.Qion.cwiseAbs2()).sum();
        else
            for (int k = 0; k < d.Qion.size(); k++)
                E += d.Aks[k] * std::norm(d.Qion[k]);
        return 2 * pc::pi / spc->geo.getVolume() * E * d.lB;
    }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] Ewald - IonIonPolicy") {
    using doctest::Approx;
    Space spc;
    spc.p.resize(2);
    spc.geo = R"( {"type": "cuboid", "length": 10} )"_json;
    spc.p[0] = R"( {"pos": [0,0,0], "q": 1.0} )"_json;
    spc.p[1] = R"( {"pos": [1,0,0], "q": -1.0} )"_json;
    Group<Particle> g(spc.p.begin(), spc.p.end());
    spc.groups.push_back(g);

    PolicyIonIon<> ionion(spc);
    EwaldData data = R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "kcutoff": 11.0, "spherical_sum": true, "cutoff": 5.0})"_json;
    Change c;
    c.all = true;
    data.ipbc = false; // PBC Ewald (http://dx.doi.org/10.1063/1.481216)
    data.update(spc.geo.getLength());
    ionion.updateComplex(data);
    CHECK(ionion.selfEnergy(data, c) == Approx(-1.0092530088080642 * data.lB));
    CHECK(ionion.surfaceEnergy(data, c) == Approx(0.0020943951023931952 * data.lB));
    CHECK(ionion.reciprocalEnergy(data) == Approx(0.21303063979675319 * data.lB));

    data.ipbc = true; // IPBC Ewald
    data.update(spc.geo.getLength());
    ionion.updateComplex(data);
    CHECK(ionion.selfEnergy(data, c) == Approx(-1.0092530088080642 * data.lB));
    CHECK(ionion.surfaceEnergy(data, c) == Approx(0.0020943951023931952 * data.lB));
    CHECK(ionion.reciprocalEnergy(data) == Approx(0.0865107467 * data.lB));
}
#endif

/** @brief Ewald summation reciprocal energy */
template <class Policy = PolicyIonIon<>> class Ewald : public Energybase {
  private:
    EwaldData data;
    Policy policy;
    Space &spc;

  public:
    Ewald(const json &j, Space &spc) : policy(spc), spc(spc) {
        name = "ewald";
        data = j;
        init();
    }

    void init() override {
        data.update(spc.geo.getLength());
        policy.updateComplex(data); // brute force. todo: be selective
    }

    double energy(Change &change) override {
        double u = 0;
        if (change) {
            // If the state is NEW (trial state), then update all k-vectors
            if (key == NEW) {
                if (change.all || change.dV) { // everything changes
                    data.update(spc.geo.getLength());
                    policy.updateComplex(data); // update all (expensive!)
                } else {
                    if (change.groups.size() > 0)
                        policy.updateComplex(data, change);
                }
            }
            u = policy.surfaceEnergy(data, change) + policy.reciprocalEnergy(data) + policy.selfEnergy(data, change);
        }
        return u;
    }

    void sync(Energybase *basePtr, Change &) override {
        auto other = dynamic_cast<decltype(this)>(basePtr);
        assert(other);
        if (other->key == OLD)
            policy.old = &(other->spc); // give NEW access to OLD space for optimized updates
        data = other->data;             // copy everything!

    } //!< Called after a move is rejected/accepted as well as before simulation

    void to_json(json &j) const override { j = data; }
};

/** @brief Self-energy term of electrostatic potentials  */
class SelfEnergy : public Energybase {
  private:
    std::string type;
    double selfenergy_ion_prefactor, selfenergy_dipole_prefactor, epsr, lB, rc, alpha, kappa, epsrf;
    int C, D;
    Space &spc;

  public:
    SelfEnergy(const json &j, Space &spc);
    double energy(Change &change) override;
};

class ParticleSelfEnergy : public Energybase {
  private:
    Space &spc;
    Potential::PairPotentialBase &pairpot;

  public:
    ParticleSelfEnergy(Space &, Potential::PairPotentialBase &);
    double energy(Change &change) override;
};

class Isobaric : public Energybase {
  private:
    Space &spc;
    double P; // P/kT
  public:
    Isobaric(const json &j, Space &spc);
    double energy(Change &change) override;
    void to_json(json &j) const override;
};

/**
 * @brief Constrain system using reaction coordinates
 *
 * If outside specified `range`, infinity energy is returned, causing rejection.
 */
class Constrain : public Energybase {
  private:
    std::string type;
    std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> rc = nullptr;

  public:
    Constrain(const json &j, Space &spc);
    double energy(Change &change) override;
    void to_json(json &j) const override;
};

/*
 * The keys of the `intra` map are group index and the values
 * is a vector of `BondData`. For bonds between groups, fill
 * in `inter` which is evaluated for every update of call to
 * `energy`.
 *
 * @todo Optimize.
 */
class Bonded : public Energybase {
  private:
    Space &spc;
    typedef std::vector<std::shared_ptr<Potential::BondData>> BondVector;
    BondVector inter;                // inter-molecular bonds
    std::map<int, BondVector> intra; // intra-molecular bonds

  private:
    void update_intra();                              // finds and adds all intra-molecular bonds of active molecules
    double sum_energy(const BondVector &bonds) const; // sum energy in vector of BondData
    double sum_energy(const BondVector &bonds, const std::vector<int> &particles_ndx)
        const; // sum energy in vector of BondData for matching particle indices

  public:
    Bonded(const json &j, Space &spc);
    void to_json(json &j) const override;
    double energy(Change &change) override; // brute force -- refine this!
};

/**
 * @brief Nonbonded energy using a pair-potential
 */
template <typename Tpairpot> class Nonbonded : public Energybase {
  private:
    double g2gcnt = 0, g2gskip = 0;
    PairMatrix<double> cutoff2; // matrix w. group-to-group cutoff

  protected:
    typedef typename Space::Tgroup Tgroup;
    double Rc2_g2g = pc::infty;

    // control of when OpenMP should be used
    bool omp_enable = false;
    bool omp_i2all = false;
    bool omp_g2g = false;
    bool omp_p2p = false;

    void to_json(json &j) const override {
        j["pairpot"] = pairpot;
        if (omp_enable) {
            json _a = json::array();
            if (omp_p2p)
                _a.push_back("p2p");
            if (omp_g2g)
                _a.push_back("g2g");
            if (omp_i2all)
                _a.push_back("i2all");
            j["openmp"] = _a;
        }
        j["cutoff_g2g"] = json::object();
        auto &_j = j["cutoff_g2g"];
        for (auto &a : Faunus::molecules)
            for (auto &b : Faunus::molecules)
                if (a.id() >= b.id())
                    _j[a.name + " " + b.name] = sqrt(cutoff2(a.id(), b.id()));
    }

    template <typename T> inline bool cut(const T &g1, const T &g2) {
        g2gcnt++;
        if (g1.atomic || g2.atomic)
            return false;
        if (spc.geo.sqdist(g1.cm, g2.cm) < cutoff2(g1.id, g2.id))
            return false;
        g2gskip++;
        return true;
    } //!< true if group<->group interaction can be skipped

    template <typename T> inline double i2i(const T &a, const T &b) {
        assert(&a != &b && "a and b cannot be the same particle");
        return pairpot(a, b, spc.geo.vdist(a.pos, b.pos));
    }

    /*
     * Internal energy in group, calculating all with all or, if `index`
     * is given, only a subset. Index specifies the internal index (starting
     * from zero) of changed particles within the group.
     */
    double g_internal(const Tgroup &g, const std::vector<int> &index = std::vector<int>()) {
        using namespace ranges;
        double u = 0;
        if (index.empty() and not molecules.at(g.id).rigid) // assume that all atoms have changed
            for (auto i = g.begin(); i != g.end(); ++i)
                for (auto j = i; ++j != g.end();)
                    u += i2i(*i, *j);
        else { // only a subset has changed
            auto fixed = view::ints(0, int(g.size())) |
                         view::remove_if([&index](int i) { return std::binary_search(index.begin(), index.end(), i); });
            for (int i : index) { // moved<->static
                for (int j : fixed) {
                    u += i2i(*(g.begin() + i), *(g.begin() + j));
                }
            }
            for (int i : index) // moved<->moved
                for (int j : index)
                    if (j > i)
                        u += i2i(*(g.begin() + i), *(g.begin() + j));
        }
        return u;
    }

    /*
     * Calculates the interaction energy of a particle, `i`,
     * and checks (1) if it is already part of Space, or (2)
     * external to space.
     */
    double i2all(const typename Space::Tparticle &i) {
        double u = 0;
        auto it = spc.findGroupContaining(i); // iterator to group
        if (it != spc.groups.end()) {         // check if i belongs to group in space
#pragma omp parallel for reduction(+ : u) if (omp_enable and omp_i2all)
            for (size_t ig = 0; ig < spc.groups.size(); ig++) {
                auto &g = spc.groups[ig];
                if (&g != &(*it))         // avoid self-interaction
                    if (not cut(g, *it))  // check g2g cut-off
                        for (auto &j : g) // loop over particles in other group
                            u += i2i(i, j);
            }
            for (auto &j : *it) // i with all particles in own group
                if (&j != &i)
                    u += i2i(i, j);
        } else                         // particle does not belong to any group
            for (auto &g : spc.groups) // i with all other *active* particles
                for (auto &j : g)      // (this will include only active particles)
                    u += i2i(i, j);
        return u;
    }

    /*
     * Group-to-group energy. A subset of `g1` can be given with `index` which refers
     * to the internal index (starting at zero) of the first group, `g1
     * NOTE: the interpretation of this function is extended to also consider the mutual interactions
     * of a subset of each group and in such case returns sub1 <-> 2 and !sub1<->sub2,
     * hence excluding !sub1 <-> !sub2 in comparision to calling onconstrained g2g. In absence
     * of sub1 any sub2 is ignored.
     */
    virtual double g2g(const Tgroup &g1, const Tgroup &g2, const std::vector<int> &index = std::vector<int>(),
                       const std::vector<int> &jndex = std::vector<int>()) {
        using namespace ranges;
        double u = 0;
        if (not cut(g1, g2)) {
            if (index.empty() && jndex.empty()) // if index is empty, assume all in g1 have changed
#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (omp_enable and omp_p2p)
                for (size_t i = 0; i < g1.size(); i++)
                    for (size_t j = 0; j < g2.size(); j++)
                        u += i2i(*(g1.begin() + i), *(g2.begin() + j));
            else { // only a subset of g1
                for (auto i : index)
                    for (auto j = g2.begin(); j != g2.end(); ++j)
                        u += i2i(*(g1.begin() + i), *j);
                if (not jndex.empty()) {
                    auto fixed = view::ints(0, int(g1.size())) | view::remove_if([&index](int i) {
                                     return std::binary_search(index.begin(), index.end(), i);
                                 });
                    for (auto i : jndex)     // moved2        <-|
                        for (auto j : fixed) // static1   <-|
                            u += i2i(*(g2.begin() + i), *(g1.begin() + j));
                }
            }
        }
        return u;
    }

  public:
    Space &spc;       //!< Space to operate on
    Tpairpot pairpot; //!< Pair potential

    Nonbonded(const json &j, Space &spc) : spc(spc) {
        name = "nonbonded";
        pairpot = j;

        // controls for OpenMP
        auto it = j.find("openmp");
        if (it != j.end())
            if (it->is_array())
                if (it->size() > 0) {
                    omp_enable = true;
                    for (const std::string &k : *it)
                        if (k == "g2g")
                            omp_g2g = true;
                        else if (k == "p2p")
                            omp_p2p = true;
                        else if (k == "i2all")
                            omp_i2all = true;
#ifndef _OPENMP
                    std::cerr << "warning: nonbonded requests unavailable OpenMP." << endl;
#endif
                }

        // disable all group-to-group cutoffs by setting infinity
        for (auto &i : Faunus::molecules)
            for (auto &j : Faunus::molecules)
                cutoff2.set(i.id(), j.id(), pc::infty);

        it = j.find("cutoff_g2g");
        if (it != j.end()) {
            // old style input w. only a single cutoff
            if (it->is_number()) {
                Rc2_g2g = std::pow(it->get<double>(), 2);
                for (auto &i : Faunus::molecules)
                    for (auto &j : Faunus::molecules)
                        cutoff2.set(i.id(), j.id(), Rc2_g2g);
            }
            // new style input w. multiple cutoffs between molecules
            else if (it->is_object()) {
                // ensure that there is a default, fallback cutoff
                Rc2_g2g = std::pow(it->at("default").get<double>(), 2);
                for (auto &i : Faunus::molecules)
                    for (auto &j : Faunus::molecules)
                        cutoff2.set(i.id(), j.id(), Rc2_g2g);
                // loop for space separated molecule pairs in keys
                for (auto &i : it->items()) {
                    auto v = words2vec<std::string>(i.key());
                    if (v.size() == 2) {
                        int id1 = (*findName(Faunus::molecules, v[0])).id();
                        int id2 = (*findName(Faunus::molecules, v[1])).id();
                        cutoff2.set(id1, id2, std::pow(i.value().get<double>(), 2));
                    }
                }
            }
        }
    }

    void force(std::vector<Point> &forces) override {
        auto &p = spc.p; // alias to particle vector (reference)
        assert(forces.size() == p.size() && "the forces size must match the particle size");
        for (size_t i = 0; i < p.size() - 1; i++)
            for (size_t j = i + 1; j < p.size(); j++) {
                // Point r = spc.geo.vdist(p[i].pos, p[j].pos); // minimum distance vector
                Point f; //= pairpot.force( p[i], p[j], r.squaredNorm(), r );
                forces[i] += f;
                forces[j] -= f;
            }
    }

    double energy(Change &change) override {
        using namespace ranges;
        double u = 0;

        if (change) {

            if (change.dV) {
#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (omp_enable and omp_g2g)
                for (auto i = spc.groups.begin(); i < spc.groups.end(); ++i) {
                    for (auto j = i; ++j != spc.groups.end();)
                        u += g2g(*i, *j);
                    if (i->atomic)
                        u += g_internal(*i);
                }
                return u;
            }

            // did everything change?
            if (change.all) {
#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (omp_enable and omp_g2g)
                for (auto i = spc.groups.begin(); i < spc.groups.end(); ++i) {
                    for (auto j = i; ++j != spc.groups.end();)
                        u += g2g(*i, *j);
                    u += g_internal(*i);
                }
                // more todo here...
                return u;
            }

            // if exactly ONE molecule is changed
            if (change.groups.size() == 1 && not change.dN) {
                auto &d = change.groups[0];
                auto gindex = spc.groups.at(d.index).to_index(spc.p.begin()).first;

                // exactly one atom has move
                if (d.atoms.size() == 1)
                    return i2all(spc.p.at(gindex + d.atoms[0]));

                // more atoms moved
                auto &g1 = spc.groups.at(d.index);
                // for (auto &g2 : spc.groups)
#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (omp_enable and omp_g2g)
                for (size_t i = 0; i < spc.groups.size(); i++) {
                    auto &g2 = spc.groups[i];
                    if (&g1 != &g2)
                        u += g2g(g1, g2, d.atoms);
                }
                if (d.internal)
                    u += g_internal(g1, d.atoms);
                return u;
            }

            auto moved = change.touchedGroupIndex(); // index of moved groups
            auto fixed = view::ints(0, int(spc.groups.size())) | view::remove_if([&moved](int i) {
                             return std::binary_search(moved.begin(), moved.end(), i);
                         }); // index of static groups

            if (change.dN) {
                /*auto moved = change.touchedGroupIndex(); // index of moved groups
                std::vector<int> Moved;
                for (auto i: moved) {
                    Moved.push_back(i);
                }
                std::sort( Moved.begin(), Moved.end() );
                auto fixed = view::ints( 0, int(spc.groups.size()) )
                    | view::remove_if(
                            [&Moved](int i){return std::binary_search(Moved.begin(), Moved.end(), i);}
                            ); // index of static groups*/
                for (auto cg1 = change.groups.begin(); cg1 < change.groups.end();
                     ++cg1) {                              // Loop over all changed groups
                    std::vector<int> ifiltered, jfiltered; // Active atoms
                    auto g1 = &spc.groups.at(cg1->index);
                    for (auto i : cg1->atoms) {
                        if (i < g1->size())
                            ifiltered.push_back(i);
                    }
                    // Skip if the group is empty
                    if (not ifiltered.empty())
                        for (auto j : fixed)
                            u += g2g(*g1, spc.groups[j], ifiltered, jfiltered);

                    for (auto cg2 = cg1; ++cg2 != change.groups.end();) {
                        for (auto i : cg2->atoms)
                            if (i < spc.groups.at(cg2->index).size())
                                jfiltered.push_back(i);
                        // Skip if both groups are empty
                        if (not(ifiltered.empty() && jfiltered.empty()))
                            u += g2g(*g1, spc.groups.at(cg2->index), ifiltered, jfiltered);
                        jfiltered.clear();
                    }
                    if (not ifiltered.empty() and not molecules.at(g1->id).rigid) {
                        if (cg1->all) {
                            u += g_internal(*g1);
                        } else
                            u += g_internal(*g1, ifiltered);
                    }
                }
                return u;
            }

            // moved<->moved
            if (change.moved2moved) {
                for (auto i = moved.begin(); i != moved.end(); ++i)
                    for (auto j = i; ++j != moved.end();)
                        u += g2g(spc.groups[*i], spc.groups[*j]);
            }

            // moved<->static
            if (omp_enable and omp_g2g) {
                std::vector<std::pair<int, int>> pairs(size(moved) * size(fixed));
                size_t cnt = 0;
                for (auto i : moved)
                    for (auto j : fixed)
                        pairs[cnt++] = {i, j};
#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (omp_enable and omp_g2g)
                for (size_t i = 0; i < pairs.size(); i++)
                    u += g2g(spc.groups[pairs[i].first], spc.groups[pairs[i].second]);
            } else
                for (auto i : moved)
                    for (auto j : fixed)
                        u += g2g(spc.groups[i], spc.groups[j]);

            // more todo!
        }
        return u;
    }

}; //!< Nonbonded, pair-wise additive energy term

template <typename Tpairpot> class NonbondedCached : public Nonbonded<Tpairpot> {
  private:
    typedef Nonbonded<Tpairpot> base;
    typedef typename Space::Tgroup Tgroup;
    Eigen::MatrixXf cache;
    Space &spc;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
    double g2g(const Tgroup &g1, const Tgroup &g2, const std::vector<int> &index = std::vector<int>(),
               const std::vector<int> &jndex = std::vector<int>()) override {
#pragma GCC diagnostic pop
        int i = &g1 - &base::spc.groups.front();
        int j = &g2 - &base::spc.groups.front();
        if (j < i)
            std::swap(i, j);
        if (base::key == Energybase::NEW) { // if this is from the trial system,
            double u = 0;
            if (not base::cut(g1, g2)) {
                for (auto &i : g1)
                    for (auto &j : g2)
                        u += base::i2i(i, j);
            }
            cache(i, j) = u;
        }
        return cache(i, j); // return (cached) value
    }

  public:
    NonbondedCached(const json &j, Space &spc) : base(j, spc), spc(spc) {
        base::name += "EM";
        init();
    }

    void init() override {
        cache.resize(spc.groups.size(), spc.groups.size());
        cache.setZero();
        for (auto i = base::spc.groups.begin(); i < base::spc.groups.end(); ++i) {
            for (auto j = i; ++j != base::spc.groups.end();) {
                int k = &(*i) - &base::spc.groups.front();
                int l = &(*j) - &base::spc.groups.front();
                if (l < k)
                    std::swap(k, l);
                double u = 0;
                if (!base::cut(*i, *j)) {
                    for (auto &k : *i)
                        for (auto &l : *j)
                            u += base::i2i(k, l);
                }
                cache(k, l) = u;
            }
        }
    } //!< Cache pair interactions in matrix

    double energy(Change &change) override {
        using namespace ranges;
        double u = 0;

        if (change) {

            if (change.all || change.dV) {
#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (this->omp_enable)
                for (auto i = base::spc.groups.begin(); i < base::spc.groups.end(); ++i) {
                    for (auto j = i; ++j != base::spc.groups.end();)
                        u += g2g(*i, *j);
                }
                return u;
            }

            // if exactly ONE molecule is changed
            if (change.groups.size() == 1) {
                auto &d = change.groups[0];
                auto &g1 = base::spc.groups.at(d.index);

#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (this->omp_enable and this->omp_g2g)
                for (size_t i = 0; i < spc.groups.size(); i++) {
                    auto &g2 = spc.groups[i];
                    if (&g1 != &g2)
                        u += g2g(g1, g2, d.atoms);
                }
                return u;
            }

            auto moved = change.touchedGroupIndex(); // index of moved groups
            auto fixed = view::ints(0, int(base::spc.groups.size())) | view::remove_if([&moved](int i) {
                             return std::binary_search(moved.begin(), moved.end(), i);
                         }); // index of static groups

            // moved<->moved
            if (change.moved2moved)
                for (auto i = moved.begin(); i != moved.end(); ++i)
                    for (auto j = i; ++j != moved.end();)
                        u += g2g(base::spc.groups[*i], base::spc.groups[*j]);
            // moved<->static
            if (this->omp_enable and this->omp_g2g) {
                std::vector<std::pair<int, int>> pairs(size(moved) * size(fixed));
                size_t cnt = 0;
                for (auto i : moved)
                    for (auto j : fixed)
                        pairs[cnt++] = {i, j};
#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (this->omp_enable and this->omp_g2g)
                for (size_t i = 0; i < pairs.size(); i++)
                    u += g2g(spc.groups[pairs[i].first], spc.groups[pairs[i].second]);
            } else
                for (auto i : moved)
                    for (auto j : fixed)
                        u += g2g(base::spc.groups[i], base::spc.groups[j]);

            // more todo!
        }
        return u;
    }

    void sync(Energybase *basePtr, Change &change) override {
        auto other = dynamic_cast<decltype(this)>(basePtr);
        assert(other);
        if (change.all || change.dV)
            cache.triangularView<Eigen::StrictlyUpper>() =
                (other->cache).template triangularView<Eigen::StrictlyUpper>();
        else
            for (auto &d : change.groups) {
                for (int i = 0; i < d.index; i++)
                    cache(i, d.index) = other->cache(i, d.index);
                for (size_t i = d.index + 1; i < base::spc.groups.size(); i++)
                    cache(d.index, i) = other->cache(d.index, i);
            }
    } //!< Copy energy matrix from other
};    //!< Nonbonded with cached energies (Energy Matrix)

#ifdef ENABLE_POWERSASA
/*
 * @todo:
 * - can only a subset of sasa be calculated? Note that it's the
 *   `update_coord()` function that takes up most time.
 * - delegate to GPU? In the PowerSasa paper this is mentioned
 */
class SASAEnergy : public Energybase {
  public:
    std::vector<float> sasa, radii;

  private:
    Space &spc;
    double probe;            // sasa probe radius (angstrom)
    double conc = 0;         // co-solute concentration (mol/l)
    Average<double> avgArea; // average surface area
    std::shared_ptr<POWERSASA::PowerSasa<float, Point>> ps = nullptr;

    void updateSASA(const ParticleVector &p);
    void to_json(json &j) const override;

    /*
     * @note
     * This is not enough as the PowerSasa object contains data
     * that also need syncing. It works due to the `update` (expensive!)
     * call in `energy`.
     */
    void sync(Energybase *basePtr, Change &c) override;

  public:
    SASAEnergy(const json &j, Space &spc);
    void init() override;
    double energy(Change &) override;
}; //!< SASA energy from transfer free energies
#endif

struct Example2D : public Energybase {
    Point &i; // reference to 1st particle in the system
    Example2D(const json &, Space &spc);
    double energy(Change &change) override;
};

class Hamiltonian : public Energybase, public BasePointerVector<Energybase> {
  protected:
    double maxenergy = pc::infty; //!< Maximum allowed energy change
    void to_json(json &j) const override;
    void addEwald(const json &j, Space &spc); //!< Adds an instance of reciprocal space Ewald energies (if appropriate)
    void
    addSelfEnergy(const json &j,
                  Space &spc); //!< Adds an instance of the self term of the electrostatic potential (if appropriate)

  public:
    Hamiltonian(Space &spc, const json &j);
    double energy(Change &change) override; //!< Energy due to changes
    void init() override;
    void sync(Energybase *basePtr, Change &change) override;
}; //!< Aggregates and sum energy terms

} // namespace Energy
} // namespace Faunus
