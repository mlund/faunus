#pragma once

#include "bonds.h"
#include "externalpotential.h" // Energybase implemented here
#include "space.h"
#include "aux/iteratorsupport.h"
#include <range/v3/view.hpp>
#include <Eigen/Dense>
#include "spdlog/spdlog.h"

#ifdef ENABLE_FREESASA
#include <freesasa.h>
#endif

namespace Faunus {

namespace ReactionCoordinate {
struct ReactionCoordinateBase;
}

namespace Potential {
struct PairPotentialBase;
}

namespace Energy {

class Hamiltonian;

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
 * @brief Data class for Ewald k-space calculations
 *
 * Currently, the Eigen policies map to the non-eigen
 * variants, e.g. `PBCEigen == PBC`.
 *
 * Related reading:
 * - PBC Ewald (DOI:10.1063/1.481216)
 * - IPBC Ewald (DOI:10/css8)
 * - Update optimization (DOI:10.1063/1.481216, Eq. 24)
 */
struct EwaldData {
    typedef std::complex<double> Tcomplex;
    Eigen::Matrix3Xd k_vectors;             //!< k-vectors, 3xK
    Eigen::VectorXd Aks;                    //!< 1xK for update optimization (see Eq.24, DOI:10.1063/1.481216)
    Eigen::VectorXcd Q_ion, Q_dipole;       //!< Complex 1xK vectors
    double r_cutoff = 0;                    //!< Real-space cutoff
    double n_cutoff = 0;                    //!< Inverse space cutoff
    double surface_dielectric_constant = 0; //!< Surface dielectric constant;
    double bjerrum_length = 0;              //!< Bjerrum length
    double kappa = 0;                       //!< Inverse Debye screening length
    double kappa_squared = 0;               //!< Squared inverse Debye screening length
    double alpha = 0;
    double const_inf = 0;
    double check_k2_zero = 0;
    bool use_spherical_sum = true;
    int num_kvectors = 0;
    Point box_length = {0.0, 0.0, 0.0};                        //!< Box dimensions
    enum Policies { PBC, PBCEigen, IPBC, IPBCEigen, INVALID }; //!< Possible k-space updating schemes
    Policies policy = PBC;                                     //!< Policy for updating k-space
    EwaldData(const json &);                                   //!< Initialize from json
};

NLOHMANN_JSON_SERIALIZE_ENUM(EwaldData::Policies, {
                                                      {EwaldData::INVALID, nullptr},
                                                      {EwaldData::PBC, "PBC"},
                                                      {EwaldData::PBCEigen, "PBCEigen"},
                                                      {EwaldData::IPBC, "IPBC"},
                                                      {EwaldData::IPBCEigen, "IPBCEigen"},
                                                  })

void to_json(json &, const EwaldData &);

/**
 * @brief Base class for Ewald k-space updates policies
 */
class EwaldPolicyBase {
  public:
    std::string cite; //!< Optional reference, preferably DOI, to further information
    virtual ~EwaldPolicyBase() = default;
    virtual void updateBox(EwaldData &, const Point &) const = 0; //!< Prepare k-vectors according to given box vector
    virtual void updateComplex(EwaldData &,
                               Space::Tgvec &) const = 0; //!< Update all k vectors
    virtual void updateComplex(EwaldData &, Change &, Space::Tgvec &,
                               Space::Tgvec &) const = 0; //!< Update subset of k vectors. Require `old` pointer
    virtual double selfEnergy(const EwaldData &, Change &,
                              Space::Tgvec &) = 0; //!< Self energy contribution due to a change
    virtual double surfaceEnergy(const EwaldData &, Change &,
                                 Space::Tgvec &) = 0;       //!< Surface energy contribution due to a change
    virtual double reciprocalEnergy(const EwaldData &) = 0; //!< Total reciprocal energy

    /**
     * @brief Represent charges and positions using an Eigen facade (Map)
     *
     * Requires that all groups are fully active, i.e. does not work for GCMC.
     *
     * @param groups Vector of groups to represent
     * @return tuple with positions, charges
     */
    auto mapGroupsToEigen(Space::Tgvec &groups) const {
        for (auto &g : groups)
            if (g.size() != g.capacity())
                throw std::runtime_error("Eigen optimized Ewald not available with inactive groups");
        auto first_particle = groups.front().begin();
        auto last_particle = groups.back().end();
        auto pos = asEigenMatrix(first_particle, last_particle,
                                 &Space::Tparticle::pos); // N x 3
        auto charge = asEigenVector(first_particle, last_particle,
                                    &Space::Tparticle::charge); // N x 1
        return std::make_tuple(pos, charge);
    }

    static std::shared_ptr<EwaldPolicyBase> makePolicy(EwaldData::Policies); //!< Policy factory
};

/**
 * @brief Ion-Ion Ewald using periodic boundary conditions (PBC)
 */
struct PolicyIonIon : public EwaldPolicyBase {
    PolicyIonIon();
    void updateBox(EwaldData &, const Point &) const override;
    void updateComplex(EwaldData &, Space::Tgvec &) const override;
    void updateComplex(EwaldData &, Change &, Space::Tgvec &, Space::Tgvec &) const override;
    double selfEnergy(const EwaldData &, Change &, Space::Tgvec &) override;
    double surfaceEnergy(const EwaldData &, Change &, Space::Tgvec &) override;
    double reciprocalEnergy(const EwaldData &) override;
};

/**
 * @brief Ion-Ion Ewald with periodic boundary conditions (PBC) using Eigen
 * operations
 * @warning Will not work with Space with inactive particles (GCMC, for example)
 *
 * For compilers that offer good vectorization (gcc on linux) this brings a 4-5
 * fold speed increase.
 * Status on February, 2020:
 * - Clang9: Eigen version is slower than generic version (macos/ubuntu)
 * - GCC9: Eigen is 4-5 times faster on x86 linux; ~1.5 times *lower on macos.
 */
struct PolicyIonIonEigen : public PolicyIonIon {
    void updateComplex(EwaldData &, Space::Tgvec &) const override;
    double reciprocalEnergy(const EwaldData &) override;
};

/**
 * @brief Ion-Ion Ewald with isotropic periodic boundary conditions (IPBC)
 */
struct PolicyIonIonIPBC : public PolicyIonIon {
    PolicyIonIonIPBC();
    void updateBox(EwaldData &, const Point &) const override;
    void updateComplex(EwaldData &, Space::Tgvec &) const override;
    void updateComplex(EwaldData &, Change &, Space::Tgvec &, Space::Tgvec &) const override;
};

/**
 * @brief Ion-Ion Ewald with isotropic periodic boundary conditions (IPBC) using Eigen operations
 * @warning Incomplete and under construction
 */
struct PolicyIonIonIPBCEigen : public PolicyIonIonIPBC {
    void updateComplex(EwaldData &, Space::Tgvec &) const override;
};

/** @brief Ewald summation reciprocal energy */
class Ewald : public Energybase {
  private:
    EwaldData data;
    std::shared_ptr<EwaldPolicyBase> policy; //!< Policy for updating k-space
    Space &spc;
    Space::Tgvec *old_groups = nullptr;

  public:
    Ewald(const json &, Space &);
    void init() override;
    double energy(Change &) override;
    void sync(Energybase *,
              Change &) override; //!< Called after a move is rejected/accepted
                                  //! as well as before simulation
    void to_json(json &) const override;
};

class Isobaric : public Energybase {
  private:
    Space &spc;
    double P; // P/kT
  public:
    Isobaric(const json &, Space &);
    double energy(Change &) override;
    void to_json(json &) const override;
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
    Constrain(const json &, Space &);
    double energy(Change &) override;
    void to_json(json &) const override;
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
    typedef BasePointerVector<Potential::BondData> BondVector;
    BondVector inter;                // inter-molecular bonds
    std::map<int, BondVector> intra; // intra-molecular bonds

  private:
    void update_intra();                              // finds and adds all intra-molecular bonds of active molecules
    double sum_energy(const BondVector &) const;      // sum energy in vector of BondData
    double sum_energy(const BondVector &,
                      const std::vector<int> &) const; // sum energy in vector of BondData for matching particle indices

  public:
    Bonded(const json &, Space &);
    void to_json(json &) const override;
    double energy(Change &) override; // brute force -- refine this!
};

/**
 * @brief Determine if two groups are separated beyond the cutoff distance.
 *
 * The cutoff distance can be specified independently for each group pair to override the default value.
 *
 * @see PairEnergy
 */
class Cutoff {
    double default_cutoff_squared = pc::infty;
    PairMatrix<double> cutoff_squared;  //!< matrix with group-to-group cutoff distances squared in angstrom squared
    double total_cnt = 0, skip_cnt = 0; //!< statistics
    Space &spc;
    friend void from_json(const json &j, Cutoff &c);
    friend void to_json(json &j, const Cutoff &c);

  public:
    /**
     * @return true if group<->group interaction is beyond the cutoff distance, i.e., it can be skipped, false otherwise
     */
    template <typename T> inline bool cut(const T &group1, const T &group2) {
        bool result = false;
        ++total_cnt;
        if (!group1.atomic && !group2.atomic // atomic groups have no meaningful cm
            && spc.geo.sqdist(group1.cm, group2.cm) >= cutoff_squared(group1.id, group2.id)) {
            result = true;
            ++skip_cnt;
        }
        return result;
    }

    /**
     * @brief Functor alias for @see cut
     */
    template <typename... Args> inline auto operator()(Args &&... args) { return cut(std::forward<Args>(args)...); }

    Cutoff(Space &spc) : spc(spc) {}
};

void from_json(const json &j, Cutoff &c);
void to_json(json &j, const Cutoff &c);

/**
 * @brief Nonbonded energy using a pair-potential
 * @tparam Tpairpot Pair potential to inject
 * @tparam allow_anisotropic_pair_potential If true, a distance vector is passed to the pair potential
 */
template <typename Tpairpot, bool allow_anisotropic_pair_potential = true> class Nonbonded : public Energybase {
  private:
    std::vector<const Particle *> i_interact_with_these;

    TimeRelativeOfTotal<> timer_g_internal;

  protected:
    typedef typename Space::Tgroup Tgroup;
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
        j["timings"] = {{"g_internal", timer_g_internal.result()}};
        Energy::to_json(j, cut);
    }

    Cutoff cut;

    template <typename T> inline double i2i(const T &a, const T &b) {
        assert(&a != &b); // a and b cannot be the same particle
        if constexpr (allow_anisotropic_pair_potential) {
            Point r = spc.geo.vdist(a.pos, b.pos);
            return pairpot(a, b, r.squaredNorm(), r);
        } else {
            return pairpot(a, b, spc.geo.sqdist(a.pos, b.pos), {0, 0, 0});
        }
    }

    /**
     * Internal energy in group, calculating all with all or, if `index`
     * is given, only a subset. Index specifies the internal index (starting
     * from zero) of changed particles within the group.
     *
     * @todo investigate overhead of `fixed_index` filtering
     */
    double g_internal(const Tgroup &g, const std::vector<int> &moved_index = std::vector<int>()) {
        timer_g_internal.start(); // measure time spent in this subroutine
        using namespace ranges;
        double u = 0;
        auto &moldata = g.traits();
        if (moldata.rigid)
            return u;

        switch (moved_index.size()) {

        // if zero, assume all particles have changed
        case 0:
            for (int i = 0; i < (int)g.size() - 1; i++)
                for (int j = i + 1; j < (int)g.size(); j++)
                    if (not moldata.isPairExcluded(i, j))
                        u += i2i(g[i], g[j]);
            break; // exit switch

        // a single particle has changed in an atomic group; no exclusion check
        case 1:
            if (g.atomic) {
                size_t i = moved_index[0];
                for (size_t j = 0; j < i; j++)
                    u += i2i(g[i], g[j]);
                for (size_t j = i + 1; j < g.size(); j++)
                    u += i2i(g[i], g[j]);
                break; // exit switch only if atomic group, otherwise continue to default:
            }
            [[fallthrough]];

        // one or more particle(s) have changed
        default:
            auto fixed_index = ranges::views::ints(0, int(g.size())) | ranges::views::remove_if([&moved_index](int i) {
                                   return std::binary_search(moved_index.begin(), moved_index.end(), i);
                               }); // index of all static particles
            for (int i : moved_index) {
                for (int j : fixed_index) // moved <-> static
                    if (not moldata.isPairExcluded(i, j))
                        u += i2i(g[i], g[j]);
                for (int j : moved_index) // moved <-> moved
                    if (j > i)
                        if (not moldata.isPairExcluded(i, j))
                            u += i2i(g[i], g[j]);
            }
        }
        timer_g_internal.stop();
        return u;
    }

    /*
     * Calculates the interaction energy of a particle, `i`,
     * and checks (1) if it is already part of Space, or (2)
     * external to space.
     */
    double i2all(const typename Space::Tparticle &i) {
        if (omp_enable and omp_i2all) {
            return i2all_parallel(i);
        }
        double u = 0;
        auto it = spc.findGroupContaining(i); // iterator to group
        if (it != spc.groups.end()) {         // check if i belongs to group in space
            for (size_t ig = 0; ig < spc.groups.size(); ig++) {
                auto &g = spc.groups[ig];
                if (&g != &(*it))         // avoid self-interaction
                    if (not cut(g, *it))  // check g2g cut-off
                        for (auto &j : g) // loop over particles in other group
                            u += i2i(i, j);
            }
            std::ptrdiff_t i_ndx = &i - &(*(it->begin()));   // fixme c++ style
            u += g_internal(*it, {static_cast<int>(i_ndx)}); // only int indices are used internally
        } else {                          // particle does not belong to any group
            for (auto &g : spc.groups) {  // i with all other *active* particles
                for (auto &j : g) {       // (this will include only active particles)
                    u += i2i(i, j);
                }
            }
        }
        return u;
    }

    double i2all_parallel(const Particle &p_i) {
        using ranges::cpp20::views::filter;
        using ranges::cpp20::views::transform;

        i_interact_with_these.clear(); // vector of particle pointers to interact with
        double u = 0;
        auto own_group = spc.findGroupContaining(p_i); // iterator to group containing p_i
        if (own_group != spc.groups.end()) {           // check if p_i belongs to group in space
            for (auto &g : spc.groups) {
                if (&g != &(*own_group))        // avoid self-interaction
                    if (not cut(g, *own_group)) // check g2g cut-off
                        std::transform(g.begin(), g.end(), std::back_inserter(i_interact_with_these),
                                       [](auto &p) { return &p; });
            }
            // p_i with all particles in own group - avoid self interaction
            auto f = *own_group | filter([&](auto &j) { return &j != &p_i; }) | transform([](auto &j) { return &j; });
            i_interact_with_these.insert(i_interact_with_these.end(), f.begin(), f.end());
        } else {                                         // particle does not belong to any group
            for (auto &g : spc.groups)                   // i with all other *active* particles
                for (auto &j : g)                        // (this will include only active particles)
                    i_interact_with_these.push_back(&j); // u += i2i(i, j);
        }

#pragma omp parallel for reduction(+ : u) if (omp_enable and omp_i2all)
        for (size_t k = 0; k < i_interact_with_these.size(); k++)
            u += i2i(p_i, *i_interact_with_these[k]);
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
                        u += i2i(g1[i], g2[j]);
            else { // only a subset of g1
                for (auto i : index)
                    for (auto j = g2.begin(); j != g2.end(); ++j)
                        u += i2i(g1[i], *j);
                if (not jndex.empty()) {
                    auto fixed = ranges::views::ints(0, int(g1.size())) | ranges::views::remove_if([&index](int i) {
                                     return std::binary_search(index.begin(), index.end(), i);
                                 });
                    for (auto i : jndex)     // moved2        <-|
                        for (auto j : fixed) // static1   <-|
                            u += i2i(g2[i], g1[j]);
                }
            }
        }
        return u;
    }

    // add self energy term to Hamiltonian if appropriate
    void addPairPotentialSelfEnergy() {
        if (pairpot.selfEnergy) { // only add if self energy is defined
            faunus_logger->debug("{} is adding self-energy from {} to Hamiltonian", name, pairpot.name);
            pot.emplace_back<Energy::ParticleSelfEnergy>(spc, pairpot.selfEnergy);
        }
    }

    void configureOpenMP(const json &j) {
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
                    faunus_logger->warn("{}: requested openmp unavailable", name);
#endif
                }
    }

  public:
    Space &spc; //!< Space to operate on
    BasePointerVector<Energybase> &pot;
    Tpairpot pairpot; //!< Pair potential

    Nonbonded(const json &j, Space &spc, BasePointerVector<Energybase> &pot) : cut(spc), spc(spc), pot(pot) {
        name = "nonbonded";
        pairpot.from_json(j);

        if (!pairpot.isotropic and !allow_anisotropic_pair_potential)
            throw std::runtime_error("Only isotropic pair potentials are allowed");

        // some pair-potentials give rise to self-energies (Wolf etc.)
        // which are added here if needed
        addPairPotentialSelfEnergy();

        configureOpenMP(j);
        from_json(j, cut);
    }

    /**
     * Calculates the force on all particles
     * @todo Change to reflect only active particle, see Space::activeParticles()
     */
    void force(std::vector<Point> &forces) override {
        auto &p = spc.p; // alias to particle vector (reference)
        assert(forces.size() == p.size() && "the forces size must match the particle size");
        for (size_t i = 0; i < p.size() - 1; i++) {
            for (size_t j = i + 1; j < p.size(); j++) {
                Point r = spc.geo.vdist(p[i].pos, p[j].pos); // minimum distance vector
                Point f = pairpot.force(p[i], p[j], r.squaredNorm(), r);
                forces[i] += f;
                forces[j] -= f;
            }
        }
    }

    double energy(Change &change) override {
        double u = 0;

        if (change) {
            // there's a change in system volume
            if (change.dV) {
#pragma omp parallel for reduction(+ : u) schedule(dynamic) if (omp_enable and omp_g2g)
                for (auto i = spc.groups.begin(); i < spc.groups.end(); ++i) {
                    for (auto j = i; ++j != spc.groups.end();)
                        u += g2g(*i, *j);
                    if (i->atomic or i->compressible)
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

                // exactly one atom has moved
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
            auto fixed = ranges::views::ints(0, int(spc.groups.size())) | ranges::views::remove_if([&moved](int i) {
                             return std::binary_search(moved.begin(), moved.end(), i);
                         }); // index of static groups

            if (change.dN) {
                /*auto moved = change.touchedGroupIndex(); // index of moved groups
                std::vector<int> Moved;
                for (auto i: moved) {
                    Moved.push_back(i);
                }
                std::sort( Moved.begin(), Moved.end() );
                auto fixed = ranges::views::ints( 0, int(spc.groups.size()) )
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
                std::vector<std::pair<int, int>> pairs(size(moved) * range_size(fixed));
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
        //assert(index.empty() && "unimplemented");
        //assert(jndex.empty() && "unimplemented");

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
    NonbondedCached(const json &j, Space &spc, BasePointerVector<Energybase> &pot) : base(j, spc, pot), spc(spc) {
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
            auto fixed =
                ranges::views::ints(0, int(base::spc.groups.size())) | ranges::views::remove_if([&moved](int i) {
                    return std::binary_search(moved.begin(), moved.end(), i);
                }); // index of static groups

            // moved<->moved
            if (change.moved2moved)
                for (auto i = moved.begin(); i != moved.end(); ++i)
                    for (auto j = i; ++j != moved.end();)
                        u += g2g(base::spc.groups[*i], base::spc.groups[*j]);
            // moved<->static
            if (this->omp_enable and this->omp_g2g) {
                std::vector<std::pair<int, int>> pairs(size(moved) * range_size(fixed));
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

#ifdef ENABLE_FREESASA
/**
 * @brief Interface to the FreeSASA C-library. Experimental and unoptimized.
 *
 * https://freesasa.github.io/
 */
class SASAEnergy : public Energybase {
  public:
    std::vector<double> sasa, radii, positions;

  private:
    Space &spc;
    double cosolute_concentration;             // co-solute concentration (mol/l)
    freesasa_parameters parameters;
    Average<double> avgArea; // average surface area

    void updatePositions(const ParticleVector &p);
    void updateRadii(const ParticleVector &p);

    void updateSASA(const ParticleVector &p, const Change &change);
    void to_json(json &j) const override;
    void sync(Energybase *basePtr, Change &c) override;

  public:
    /**
     * @param spc
     * @param cosolute_concentration in particles per angstrom cubed
     * @param probe_radius in angstrom
     */
    SASAEnergy(Space &spc, double cosolute_concentration = 0.0, double probe_radius = 1.4);
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
  public:
    Hamiltonian(Space &spc, const json &j);
    double energy(Change &change) override; //!< Energy due to changes
    void init() override;
    void sync(Energybase *basePtr, Change &change) override;
}; //!< Aggregates and sum energy terms

} // namespace Energy
} // namespace Faunus
