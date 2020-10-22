#pragma once

#include "bonds.h"
#include "externalpotential.h" // Energybase implemented here
#include "space.h"
#include "aux/iteratorsupport.h"
#include "aux/pairmatrix.h"
#include <range/v3/view.hpp>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>
#include <numeric>

#ifdef ENABLE_FREESASA
#include <freesasa.h>
#endif

#ifdef __cpp_lib_parallel_algorithm
#include <execution>
#endif

namespace Faunus {

namespace ReactionCoordinate {
struct ReactionCoordinateBase;
}

namespace Potential {
struct PairPotentialBase;
}

/**
 *  @par Non-bonded energy
 *
 *  Several classes (class templates) are used together to allow computation in change of the non-bonded energy upon
 *  a MC move.
 *
 *  The energy change is calculated by the Nonbonded class. It internally uses one of the pairing policies
 *  to efficiently get all pair interactions affected by the MC move (as described by the Change object).
 *
 *  Pairing policies allow efficient summation of pair energies over the whole system, between groups, inside a group,
 *  etc. The pairing policy is optimized for performance in a different execution environment, e.g., sequential or
 *  OMP parallelism.
 *
 *  Policies have direct access to the pair interaction energy functor represented by a simple PairEnergy template.
 *  Furthermore, the GroupCutoff object is provided to limit free energy computation using a cutoff distance between
 *  respective groups.
 *
 *  @see Nonbonded, PairingBasePolicy, PairEnergy, GroupCutoff
 */
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
    using PolicyIonIon::updateComplex;
    void updateComplex(EwaldData &, Space::Tgvec &) const override;
    double reciprocalEnergy(const EwaldData &) override;
};

/**
 * @brief Ion-Ion Ewald with isotropic periodic boundary conditions (IPBC)
 */
struct PolicyIonIonIPBC : public PolicyIonIon {
    using PolicyIonIon::updateComplex;
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
    using PolicyIonIonIPBC::updateComplex;
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
    void force(std::vector<Point> &) override; // update forces on all particles
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
    std::map<int, BondVector> intra; // intra-molecular bonds; key is group index

  private:
    void update_intra();                              // finds and adds all intra-molecular bonds of active molecules
    double sum_energy(const BondVector &) const;      // sum energy in vector of BondData

    /**
     * @brief Sum energy in vector of BondData for matching particle indices
     * @param bonds List of bonds
     * @param indices_of_particles Particle index
     *
     * To speed up the bond search, the given indices must be ordered which allows
     * for binary search which on large systems provides superior performance compared
     * to simplistic search which scales with number of bonds as O(N^2).
     */
    template <class RangeOfIndex>
    double sum_energy(const Bonded::BondVector &bonds, const RangeOfIndex &indices_of_particles) const {
        assert(std::is_sorted(indices_of_particles));

        auto bond_filter = [&](const auto &bond_ptr) { // determine if bond is part of indices of particles
            for (auto index : bond_ptr->index) {
                if (std::binary_search(indices_of_particles.begin(), indices_of_particles.end(), index)) {
                    return true;
                }
            }
            return false;
        };
        auto affected_bonds = bonds | ranges::cpp20::views::filter(bond_filter);

#if (defined(__clang__) && __clang_major__ >= 10) || (defined(__GNUC__) && __GNUC__ >= 9 && __GNUC_MINOR__ >= 3)
        return std::transform_reduce(affected_bonds.begin(), affected_bonds.end(), 0.0, std::plus<>(),
                                     [&](const auto &bond_ptr) {
                                         assert(bond_ptr->hasEnergyFunction());
                                         return bond_ptr->energyFunc(spc.geo.getDistanceFunc());
                                     });
#else
        double energy = 0.0;
        for (const auto &bond_ptr : affected_bonds) {
            assert(bond_ptr->hasEnergyFunction());
            energy += bond_ptr->energyFunc(spc.geo.getDistanceFunc());
        }
        return energy;
#endif
    }

  public:
    Bonded(const json &, Space &);
    void to_json(json &) const override;
    double energy(Change &) override;          //!< brute force -- refine this!
    void force(std::vector<Point> &) override; //!< Calculates the forces on all particles
};

/**
 * @brief Provides a complementary set of ints with respect to the iota set of a given size.
 * @remark It is used as a helper function for pair interactions.
 *
 * @tparam TSize  a number castable to int
 * @tparam TSet  a finite iterable container on ints
 * @param size  the iota superset contains all integers in the range [0, size)
 * @param set  an original set of integers
 * @return a set of ints complementary to the original set
 */
template <typename TSize, typename TSet> inline auto indexComplement(const TSize size, const TSet &set) {
    assert(size <= std::numeric_limits<int>::max());
    return ranges::views::ints(0, static_cast<int>(size)) | ranges::views::remove_if([&set](TSize i) {
      return std::binary_search(set.begin(), set.end(), i);
    });
}

/**
 * @brief Determines if two groups are separated beyond the cutoff distance.
 *
 * The distance between centers of mass is considered. The cutoff distance can be specified independently for each
 * group pair to override the default value.
 *
 * @see PairEnergy
 */
class GroupCutoff {
    double default_cutoff_squared = pc::max_value;
    PairMatrix<double> cutoff_squared;  //!< matrix with group-to-group cutoff distances squared in angstrom squared
    double total_cnt = 0, skip_cnt = 0; //!< statistics
    Space::Tgeometry &geometry;         //!< geometry to compute the inter group distance with
    friend void from_json(const json&, GroupCutoff &);
    friend void to_json(json&, const GroupCutoff &);

  public:
    /**
     * @brief Determines if two groups are separated beyond the cutoff distance.
     * @return true if the group-to-group distance is beyond the cutoff distance, false otherwise
     */
    template <typename TGroup> inline bool cut(const TGroup &group1, const TGroup &group2) {
        bool result = false;
        ++total_cnt;
        if (!group1.atomic && !group2.atomic // atomic groups have no meaningful cm
            && geometry.sqdist(group1.cm, group2.cm) >= cutoff_squared(group1.id, group2.id)) {
            result = true;
            ++skip_cnt;
        }
        return result;
    }

    /**
     * @brief A functor alias for cut().
     * @see cut()
     */
    template <typename... Args> inline auto operator()(Args &&... args) { return cut(std::forward<Args>(args)...); }

    /**
     * @brief Sets the geometry.
     * @param geometry  geometry to compute the inter group distance with
     */
    GroupCutoff(Space::Tgeometry &geometry);
};

void from_json(const json&, GroupCutoff &);
void to_json(json&, const GroupCutoff &);

/**
 * @brief Provides a fast inlineable interface for non-bonded pair potential energy computation.
 *
 * @tparam TPairPotential  a pair potential to compute with
 * @tparam allow_anisotropic_pair_potential  pass also a distance vector to the pair potential, slower
 */
template <typename TPairPotential, bool allow_anisotropic_pair_potential = true> class PairEnergy {
    Space::Tgeometry &geometry;                //!< geometry to operate with
    TPairPotential pair_potential;             //!< pair potential function/functor
    Space &spc;                                //!< space to init ParticleSelfEnergy with @see addPairPotentialSelfEnergy
    BasePointerVector<Energybase> &potentials; //!< registered non-bonded potentials @see addPairPotentialSelfEnergy
  public:
    /**
     * @param spc
     * @param potentials  registered non-bonded potentials
     */
    PairEnergy(Space &spc, BasePointerVector<Energybase> &potentials) : geometry(spc.geo), spc(spc), potentials(potentials) {}

    /**
     * @brief Computes pair potential energy.
     *
     * @param a  particle
     * @param b  particle
     * @return pair potential energy between particles a and b
     */
    template <typename T> inline double potential(const T &a, const T &b) const {
        assert(&a != &b); // a and b cannot be the same particle
        if constexpr (allow_anisotropic_pair_potential) {
            const Point r = geometry.vdist(a.pos, b.pos);
            return pair_potential(a, b, r.squaredNorm(), r);
        } else {
            return pair_potential(a, b, geometry.sqdist(a.pos, b.pos), {0, 0, 0});
        }
    }

    // just a temporary placement until PairForce class template will be implemented
    template <typename T> inline Point force(const T &a, const T &b) const {
        assert(&a != &b); // a and b cannot be the same particle
        const Point r = geometry.vdist(a.pos, b.pos);
        return pair_potential.force(a, b, r.squaredNorm(), r);
    }

    /**
     * @brief A functor alias for potential().
     * @see potential()
     */
    template <typename... Args> inline auto operator()(Args &&... args) {
        return potential(std::forward<Args>(args)...);
    }

    /**
     * @brief Registers the potential self-energy to hamiltonian if needed.
     * @see Hamiltonian::Hamiltonian
     */
    void addPairPotentialSelfEnergy() {
        if (pair_potential.selfEnergy) { // only add if self energy is defined
            faunus_logger->debug("Adding self-energy from {} to hamiltonian", pair_potential.name);
            potentials.emplace_back<Energy::ParticleSelfEnergy>(spc, pair_potential.selfEnergy);
        }
    }

    void from_json(const json &j) {
        pair_potential.from_json(j);
        if (!pair_potential.isotropic && !allow_anisotropic_pair_potential) {
            throw std::logic_error("Only isotropic pair potentials are allowed.");
        }
        addPairPotentialSelfEnergy();
    }

    void to_json(json &j) const { pair_potential.to_json(j); }
};

/**
 * @brief Particle pairing to calculate non-bonded pair potential energies.
 *
 * A complete basic (serial) implementation of particle-particle pairing. The class shall not be used directly.
 * The derived template class PairingPolicy shall be used instead.
 *
 * @remark Method arguments are generally not checked for correctness because of performance reasons.
 *
 * @tparam TPairEnergy  a functor to compute non-bonded energy between two particles
 * @tparam TCutoff  a cutoff scheme between groups
 * @see PairEnergy, PairingPolicy
 */
template <typename TPairEnergy, typename TCutoff> class PairingBasePolicy {
  protected:
    Space &spc;              //!< a space to operate on
    TPairEnergy pair_energy; //!< a functor to compute non-bonded energy between two particles @see PairEnergy
    GroupCutoff cut;         //!< a cutoff functor that determines if energy between two groups can be ignored

  public:
    /**
     * @param spc
     * @param potentials  registered non-bonded potentials
     */
    PairingBasePolicy(Space &spc, BasePointerVector<Energybase> &potentials)
        : spc(spc), pair_energy(spc, potentials), cut(spc.geo) {}

    void from_json(const json &j) {
        Energy::from_json(j, cut);
        pair_energy.from_json(j);
    }

    void to_json(json &j) const {
        pair_energy.to_json(j);
        Energy::to_json(j, cut);
    }

    template <typename T> inline double particle2particle(const T &a, const T &b) const {
        return pair_energy.potential(a, b);
    }

    /**
     * @brief Internal energy of a group.
     *
     * All non-bonded pair interaction within the group are summed up. The pair exclusions defined in the molecule
     * topology are honoured.
     *
     * @param group
     * @return energy sum between particle pairs
     */
    template <typename TGroup> double groupInternal(const TGroup &group) {
        double u = 0;
        auto &moldata = group.traits();
        if (!moldata.rigid) {
            const int group_size = group.size();
            for (int i = 0; i < group_size - 1; ++i) {
                for (int j = i + 1; j < group_size; ++j) {
                    // This compound condition is faster than an outer atomic condition;
                    // tested on bulk example in GCC 9.2.
                    if (group.atomic || !moldata.isPairExcluded(i, j)) {
                        u += particle2particle(group[i], group[j]);
                    }
                }
            }
        }
        return u;
    }

    /**
     * @brief Partial internal energy of a group limited to interactions of a single particle within the group.
     *
     * The pair exclusions defined in the molecule topology are honoured.
     *
     * @param group
     * @param index  internal index of the selected particle within the group
     * @return energy sum between particle pairs
     */
    template <typename TGroup> double groupInternal(const TGroup &group, const int index) {
        double u = 0;
        auto &moldata = group.traits();
        if (!moldata.rigid) {
            if (group.atomic) {
                // speed optimization: non-bonded interaction exclusions do not need to be checked for atomic groups
                for (int i = 0; i < index; ++i) {
                    u += particle2particle(group[index], group[i]);
                }
                for (int i = index + 1; i < group.size(); ++i) {
                    u += particle2particle(group[index], group[i]);
                }
            } else {
                // molecular group
                for (int i = 0; i < index; ++i) {
                    if (!moldata.isPairExcluded(index, i)) {
                        u += particle2particle(group[index], group[i]);
                    }
                }
                for (int i = index + 1; i < group.size(); ++i) {
                    if (!moldata.isPairExcluded(index, i)) {
                        u += particle2particle(group[index], group[i]);
                    }
                }
            }
        }
        return u;
    }

    /**
     * @brief Partial internal energy of the group limited to the particles present in the index.
     *
     * Only such non-bonded pair interactions within the group are considered if at least one particle is present
     * in the index. The pair exclusions defined in the molecule topology are honoured.
     *
     * @param group
     * @param index  internal indices of particles within the group
     * @return energy sum between particle pairs
     */
    template <typename TGroup, typename TIndex> double groupInternal(const TGroup &group, const TIndex &index) {
        double u = 0;
        auto &moldata = group.traits();
        if (!moldata.rigid) {
            if (index.size() == 1) {
                u = groupInternal(group, index[0]);
            } else {
                // TODO investigate overhead of `index_complement` filtering;
                // TODO perhaps allow different strategies based on the index-size/group-size ratio
                auto index_complement = indexComplement(group.size(), index);
                // moved <-> static
                for (int i : index) {
                    for (int j : index_complement) {
                        if (!moldata.isPairExcluded(i, j)) {
                            u += particle2particle(group[i], group[j]);
                        }
                    }
                }
                // moved <-> moved
                for (auto i_it = index.begin(); i_it < index.end(); ++i_it) {
                    for (auto j_it = std::next(i_it); j_it < index.end(); ++j_it) {
                        if (!moldata.isPairExcluded(*i_it, *j_it)) {
                            u += particle2particle(group[*i_it], group[*j_it]);
                        }
                    }
                }
            }
        }
        return u;
    }

    /**
     * @brief Complete cartesian pairing of particles in two groups.
     *
     * group1 × group2
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, no calculation is performed.
     * The group intersection must be an empty set, i.e., no particle is included in both groups. This is not verified
     * for performance reason.
     *
     * @param group1
     * @param group2
     * @return energy sum between particle pairs
     */
    template <typename TGroup> double group2group(const TGroup &group1, const TGroup &group2) {
        double u = 0;
        if (!cut(group1, group2)) {
            for (auto &particle1 : group1) {
                for (auto &particle2 : group2) {
                    u += particle2particle(particle1, particle2);
                }
            }
        }
        return u;
    }

    /**
     * @brief Cross pairing of particles in two groups. Only a cartesian subset of the complete cartesian product is
     * considered as the particles in the first group must be also present in the index. The aim is to capture only
     * interactions that involve changing (indexed) particles.
     *
     * ⊕group1 × group2, where ⊕ denotes a filter by an index
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, no calculation is performed.
     * The group intersection must be an empty set, i.e., no particle is included in both groups. This is not verified
     * for performance reason.

     * @param group1
     * @param group2
     * @param index1  list of particle indices in group1 relative to the group beginning
     * @return energy sum between particle pairs
     */
    template <typename TGroup>
    double group2group(const TGroup &group1, const TGroup &group2, const std::vector<int> &index1) {
        double u = 0;
        if (!cut(group1, group2)) {
            for (auto particle1_ndx : index1) {
                for (auto &particle2 : group2) {
                    u += particle2particle(*(group1.begin() + particle1_ndx), particle2);
                }
            }
        }
        return u;
    }

    /**
     * @brief Cross pairing of particles in two groups. Only a non-cartesian subset of the complete cartesian product
     * is considered as at least one particles in the pair must be also present in the respective index. The aim is
     * to capture only interactions that involve changing (indexed) particles, i.e., to avoid pairs containing only
     * non-indexed particles.
     *
     * (⊕group1 × ∁⊕group2) + (∁⊕group1 × ⊕group2) + (⊕group1 × ⊕group2) =
     * = group1 × group2 − (∁⊕group2 × ∁⊕group2), where ⊕ denotes a filter by an index and ∁ a complement
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, no calculation is performed.
     * The group intersection must be an empty set, i.e., no particle is included in both groups. This is not verified
     * for performance reason.
     *
     * @param group1
     * @param group2
     * @param index1  list of particle indices in group1 relative to the group beginning
     * @param index2  list of particle indices in group2 relative to the group beginning
     * @return energy sum between particle pairs
     */
    template <typename TGroup>
    double group2group(const TGroup &group1, const TGroup &group2, const std::vector<int> &index1,
                       const std::vector<int> &index2) {
        double u = 0;
        if (!cut(group1, group2)) {
            if (!index2.empty()) {
                // (∁⊕group1 × ⊕group2) + (⊕group1 × ⊕group2) = group1 × ⊕group2
                u += group2group(group2, group1, index2);
                // + (⊕group1 × ∁⊕group2)
                auto index2_complement = indexComplement(group2.size(), index2);
                for (auto particle1_ndx : index1) {
                    for (auto particle2_ndx : index2_complement) {
                        u += particle2particle(group2[particle2_ndx], group1[particle1_ndx]);
                    }
                }
            } else if (!index1.empty()) {
                // (⊕group1 × ∁⊕group2) + (⊕group1 × ⊕group2) = ⊕group1 × group2
                u += group2group(group1, group2, index1);
                // + (∁⊕group1 × ⊕group2) = Ø as ⊕group2 is empty
            } else {
                // both indices empty hence nothing to do
            }
        }
        return u;
    }

    /**
     * @brief Complete cartesian pairing between particles in a group and a union of groups.
     *
     * group × (∪ groups)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing between them
     * is skipped. The internal energy of the group is not computed even if the group is also present in the union
     * of groups.
     *
     * @param group
     * @param groups
     * @return energy sum between particle pairs
     */
    template <typename TGroup, typename TGroups> double group2groups(const TGroup &group, const TGroups &groups) {
        double u = 0;
        for (auto &other_group : groups) {
            if (&other_group != &group) {
                u += group2group(group, other_group);
            }
        }
        return u;
    }

    /**
     * @brief Cross pairing of particles in a group and a union of groups. Only a cartesian subset of the complete
     * cartesian product is considered as the particles of the first group must be also present in the index. The aim
     * is to capture only interactions that involve changing (indexed) particles.
     *
     * ⊕group × (∪ groups), where ⊕ denotes a filter by an index
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped. The internal energy of the group is not computed even if the group is also present
     * in the union of groups.
     *
     * @param group
     * @param group_index  groups as indices in Space::groups
     * @param index  list of particle indices in the group relative to the group beginning
     * @return energy sum between particle pairs
     */
    template <typename TGroup, typename TGroups>
    double group2groups(const TGroup &group, const TGroups &group_index, const std::vector<int> &index) {
        double u = 0;
        for (auto other_group_ndx : group_index) {
            const auto &other_group = spc.groups[other_group_ndx];
            if (&other_group != &group) {
                u += group2group(group, other_group, index);
            }
        }
        return u;
    }

    /**
     * @brief Complete cartesian pairing between particles in a group and particles in other groups in space.
     *
     * group × (space ∖ group)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @param group
     * @return energy sum between particle pairs
     */
    template <typename Tgroup> double group2all(const Tgroup &group) {
        double u = 0;
        for (auto &other_group : spc.groups) {
            if (&other_group != &group) {
                u += group2group(group, other_group);
            }
        }
        return u;
    }

    /**
     * @brief Complete cartesian pairing between a single particle in a group and particles in other groups in space.
     *
     * ⊕group × (space ∖ group), where ⊕ denotes a filter by an index (here a single particle)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped. This method is performance-optimized version of the multiple indices method.
     *
     * @param group
     * @param index  a particle index relative to the group beginning
     * @return energy sum between particle pairs
     */
    template <typename TGroup> double group2all(const TGroup &group, const int index) {
        double u = 0;
        const auto &particle = group[index];
        for (auto &other_group : spc.groups) {
            if (&other_group != &group) {                      // avoid self-interaction
                if (!cut(other_group, group)) {                // check g2g cut-off
                    for (auto &other_particle : other_group) { // loop over particles in other group
                        u += particle2particle(particle, other_particle);
                    }
                }
            }
        }
        return u;
    }

    /**
     * @brief Complete cartesian pairing between selected particles in a group and particles in other groups in space.
     *
     * ⊕group × (space ∖ group), where ⊕ denotes a filter by an index
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @param group
     * @param index  list of particle indices in the group relative to the group beginning
     * @return energy sum between particle pairs
     */
    template <typename Tgroup> double group2all(const Tgroup &group, const std::vector<int> &index) {
        double u = 0;
        if (index.size() == 1) {
            u = group2all(group, index[0]);
        } else {
            for (auto &other_group : spc.groups) {
                if (&other_group != &group) {
                    u += group2group(group, other_group, index);
                }
            }
        }
        return u;
    }

    /**
     * @brief Cross pairing of particles among a union of groups. No internal pairs within any group are considered.
     *
     * If the distance between any two groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @param group_index  list of groups
     * @return energy sum between particle pairs
     */
    template <typename T> double groups2self(const T &group_index) {
        double u = 0;
        for (auto group1_ndx_it = group_index.begin(); group1_ndx_it < group_index.end(); ++group1_ndx_it) {
            // u += groupInternal(spc.groups[*group_it]);
            for (auto group2_ndx_it = std::next(group1_ndx_it); group2_ndx_it < group_index.end(); group2_ndx_it++) {
                u += group2group(spc.groups[*group1_ndx_it], spc.groups[*group2_ndx_it]);
            }
        }
        return u;
    }

    /**
     * @brief Cross pairing of particles between a union of groups and its complement in space.
     *
     * If the distance between any two groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @param group_index  list of groups
     * @return energy sum between particle pairs
     */
    template <typename T> double groups2all(const T &group_index) {
        double u = 0;
        u += groups2self(group_index);
        auto index_complement = indexComplement(spc.groups.size(), group_index);
        for (auto group1_ndx : group_index) {
            for (auto group2_ndx : index_complement) {
                u += group2group(spc.groups[group1_ndx], spc.groups[group2_ndx]);
            }
        }
        return u;
    }

    /**
     * @brief Cross pairing between all particles in the space.
     *
     * If the distance between particles' groups is greater or equal to the group cutoff distance, no calculation is
     * performed.
     *
     * @return energy sum between particle pairs
     */
    double all() {
        double u = 0;
        for (auto group_it = spc.groups.begin(); group_it < spc.groups.end(); ++group_it) {
            u += groupInternal(*group_it);
            for (auto other_group_it = std::next(group_it); other_group_it < spc.groups.end(); other_group_it++) {
                u += group2group(*group_it, *other_group_it);
            }
        }
        return u;
    }

    /**
     * @brief Cross pairing between all particles in the space.
     *
     * If the distance between particles' groups is greater or equal to the group cutoff distance, no calculation is
     * performed.
     *
     * @tparam TCondition  a function returning bool and having a group as an argument
     * @param condition  a group filter if internal energy of the group shall be added
     * @return energy sum between particle pairs
     */
    template <typename TCondition> double all(TCondition condition) {
        double u = 0;
        for (auto group_it = spc.groups.begin(); group_it < spc.groups.end(); ++group_it) {
            if (condition(*group_it)) {
                u += groupInternal(*group_it);
            }
            for (auto other_group_it = std::next(group_it); other_group_it < spc.groups.end(); other_group_it++) {
                u += group2group(*group_it, *other_group_it);
            }
        }
        return u;
    }

    void force(std::vector<Point> &forces) {
        // just a temporary hack; perhaps better to allow PairForce instead of the PairEnergy template
        assert(forces.size() == spc.p.size() && "the forces size must match the particle size");
        for (size_t i = 0; i < spc.p.size() - 1; ++i) {
            for (size_t j = i + 1; j < spc.p.size(); ++j) {
                const Point f = pair_energy.force(spc.p[i], spc.p[j]);
                forces[i] += f;
                forces[j] -= f;
            }
        }
    }
};

template <typename TPairEnergy, typename TCutoff, bool parallel = false>
class PairingPolicy : public PairingBasePolicy<TPairEnergy, TCutoff> {
  public:
    using PairingBasePolicy<TPairEnergy, TCutoff>::PairingBasePolicy;
};

/**
 * @brief Computes change in the non-bonded energy, assuming pair-wise additive energy terms.
 *
 * @tparam TPairingPolicy  pairing policy to effectively sum up the pair-wise additive non-bonded energy
 */
template <typename TPairingPolicy> class Nonbonded : public Energybase {
  protected:
    Space &spc;             //!< space to operate on
    TPairingPolicy pairing; //!< pairing policy to effectively sum up the pair-wise additive non-bonded energy

    /**
     * @brief Computes non-bonded energy contribution if only a single group has changed.
     *
     * @param change
     * @return energy sum between particle pairs
     */
    double energyGroup(Change &change) {
        double u = 0;
        const auto &change_data = change.groups[0];
        const auto &group = spc.groups.at(change_data.index);
        if (change_data.atoms.size() == 1) {
            // faster algorithm if only a single particle moves
            u = pairing.group2all(group, change_data.atoms[0]);
            if (change_data.internal) {
                u += pairing.groupInternal(group, change_data.atoms[0]);
            }
        } else {
            const bool change_all = change_data.atoms.empty(); // all particles or only their subset?
            u = change_all ? pairing.group2all(group) : pairing.group2all(group, change_data.atoms);
            if (change_data.internal) {
                u += change_all ? pairing.groupInternal(group) : pairing.groupInternal(group, change_data.atoms);
            }
        }
        return u;
    }

    /**
     * @brief Computes non-bonded energy difference if the number of particles have changed.
     *
     * Particles have to be explicitly enumerated in the atom indices of the changed group. Implicit addition of atoms
     * with a group is not supported yet. Note that we do not have to care about missing (removed) particles at all.
     * They are taken into account in the original (old) space where they are present.
     *
     * @param change
     * @return energy sum between particle pairs
     */
    double energySpeciation(Change &change) {
        assert(change.dN);
        double u = 0;
        const auto &moved = change.touchedGroupIndex(); // index of moved groups
        const auto &fixed = indexComplement(int(spc.groups.size()), moved) | ranges::to<std::vector>; // index of static groups
        auto filter_active = [](int size){return ranges::views::filter([size](const auto i) { return i < size; }); };

        // loop over all changed groups
        for (auto change_group1_it = change.groups.begin(); change_group1_it < change.groups.end(); ++change_group1_it) {
            auto &group1 = spc.groups.at(change_group1_it->index);
            // filter only active particles
            const std::vector<int> index1 = change_group1_it->atoms | filter_active(group1.size()) | ranges::to<std::vector>;
            if (!index1.empty()) {
                // particles added into the group: compute (changed group) <-> (static group)
                u += pairing.group2groups(group1, fixed, index1);
            }
            // loop over successor changed groups (hence avoid double counting group1×group2 and group2×group1)
            for (auto change_group2_it = std::next(change_group1_it); change_group2_it < change.groups.end(); ++change_group2_it) {
                auto &group2 = spc.groups.at(change_group2_it->index);
                const std::vector<int> index2 = change_group2_it->atoms | filter_active(group2.size()) | ranges::to<std::vector>;
                if (!index1.empty() || !index2.empty()) {
                    // particles added into one or other group: compute (changed group) <-> (changed group)
                    u += pairing.group2group(group1, group2, index1, index2);
                }
            }
            if (!index1.empty() && !molecules.at(group1.id).rigid) {
                // compute internal energy in the changed group
                u += change_group1_it->all ? pairing.groupInternal(group1) : pairing.groupInternal(group1, index1);
            }
        }
        return u;
    }

  public:
    Nonbonded(const json &j, Space &spc, BasePointerVector<Energybase> &pot) : spc(spc), pairing(spc, pot) {
        name = "nonbonded";
        pairing.from_json(j);
    }

    void to_json(json &j) const override { pairing.to_json(j); }

    /**
     * @brief Calculates the force on all particles.
     *
     * @todo A stub. Change to reflect only active particle, see Space::activeParticles().
     */
    void force(std::vector<Point> &forces) override { pairing.force(forces); }

    /**
     * @brief Computes non-bonded energy contribution from changed particles.
     *
     * The internal energy contribution, i.e., the contribution from the intra group interactions, is added
     * only if a single group is changed or if all changed.
     *
     * @param change
     * @return energy sum between particle pairs
     */
    double energy(Change &change) override {
        assert(std::is_sorted(change.groups.begin(), change.groups.end()));
        double u = 0;
        if (change.all) {
            u = pairing.all();
        } else if (change.dV) {
            // sum all interaction energies except the internal energies of incompressible molecules
            u = pairing.all([](auto &group) { return group.atomic || group.compressible; });
        } else if (!change.dN) {
            if (change.groups.size() == 1) {
                // if only a single group changes use faster algorithm and optionally add the internal energy
                u = energyGroup(change);
            } else {
                // if multiple groups move, no internal energies are computed
                const auto &moved = change.touchedGroupIndex(); // index of moved groups
                u = pairing.groups2all(moved);
            }
        } else { // change.dN
            u = energySpeciation(change);
        }
        return u;
    }
};


/**
 * @brief Computes non-bonded energy contribution from changed particles. Cache group2group energy once calculated,
 * until a new trial configuration is provided. Not for general use as only partially implemented!
 *
 * Original implementation, only refurbished. Generally suboptimal as only PairingPolicy::group2group method
 * may be called.
 * No internal energy is ever computed. Cannot deal with particle count changes. And other unmentioned constrains.
 *
 * @tparam Tpairpot
 */
template <typename Tpairpot> class NonbondedCached : public Nonbonded<PairingPolicy<PairEnergy<Tpairpot>, GroupCutoff>> {
    typedef Nonbonded<PairingPolicy<PairEnergy<Tpairpot>, GroupCutoff>> base;
    typedef typename Space::Tgroup Tgroup;
    Eigen::MatrixXf cache;
    using base::spc;

    double g2g(const Tgroup &g1, const Tgroup &g2) {
        int i = &g1 - spc.groups.data();
        int j = &g2 - spc.groups.data();
        if (j < i) {
            std::swap(i, j);
        }
        if (base::key == Energybase::TRIAL_MONTE_CARLO_STATE) { // if this is from the trial system,
            cache(i, j) = base::pairing.group2group(g1, g2); // update the cache
        }
        return cache(i, j); // return (cached) value
    }

    double g2g(const Tgroup &g1, const Tgroup &g2, [[maybe_unused]] const std::vector<int> &index) {
        // index not implemented
        return g2g(g1, g2);
    }

  public:
    NonbondedCached(const json &j, Space &spc, BasePointerVector<Energybase> &pot) : base(j, spc, pot) {
        base::name += "EM";
        init();
    }

  /**
   * @brief Cache pair interactions in matrix.
   */
    void init() override {
        const auto groups_size = spc.groups.size();
        cache.resize(groups_size, groups_size);
        cache.setZero();
        for (auto i = 0; i < groups_size - 1; ++i) {
            for (auto j = i + 1; j < groups_size; ++j) {
                cache(i, j) = base::pairing.group2group(spc.groups[i], spc.groups[j]);
            }
        }
    }

    double energy(Change &change) override {
        // Only g2g may be called there to compute (and cache) energy!
        double u = 0;
        if (change) {
            if (change.all || change.dV) {
                for (auto i = spc.groups.begin(); i < spc.groups.end(); ++i) {
                    for (auto j = std::next(i); j < base::spc.groups.end(); ++j) {
                        u += g2g(*i, *j);
                    }
                }
            } else {
                if (change.groups.size() == 1) { // if exactly ONE molecule is changed
                    auto &d = change.groups[0];
                    auto &g1 = spc.groups.at(d.index);
                    for (auto g2_it = spc.groups.begin(); g2_it < spc.groups.end(); ++g2_it) {
                        if (&g1 != &(*g2_it)) {
                            u += g2g(g1, *g2_it, d.atoms);
                        }
                    }
                } else {                                     // many molecules are changed
                    auto moved = change.touchedGroupIndex(); // index of moved groups
                    // moved<->moved
                    if (change.moved2moved) {
                        for (auto i = moved.begin(); i < moved.end(); ++i) {
                            for (auto j = std::next(i); j < moved.end(); ++j) {
                                u += g2g(spc.groups[*i], spc.groups[*j]);
                            }
                        }
                    }
                    // moved<->static
#if true
                    // classic version
                    auto fixed = indexComplement(spc.groups.size(), moved); // index of static groups
                    for (auto i : moved) {
                        for (auto j : fixed) {
                            u += g2g(spc.groups[i], spc.groups[j]);
                        }
                    }
#else
                    // OMP-ready version
                    auto fixed =
                        indexComplement(spc.groups.size(), moved) | ranges::to<std::vector>; // index of static groups
                    const size_t moved_size = moved.size();
                    const size_t fixed_size = fixed.size();
                    for (auto i = 0; i < moved_size; ++i) {
                        for (auto j = 0; j < fixed_size; ++j) {
                            u += g2g(spc.groups[moved[i]], spc.groups[fixed[j]]);
                        }
                    }
#endif
                }
            }
            // more todo!
        }
        return u;
    }

    /**
     * @brief Copy energy matrix from other
     * @param base_ptr
     * @param change
     */
    void sync(Energybase *base_ptr, Change &change) override {
        auto other = dynamic_cast<decltype(this)>(base_ptr);
        assert(other);
        if (change.all || change.dV) {
            cache.triangularView<Eigen::StrictlyUpper>() =
                (other->cache).template triangularView<Eigen::StrictlyUpper>();
        } else {
            for (auto &d : change.groups) {
                for (int i = 0; i < d.index; i++) {
                    cache(i, d.index) = other->cache(i, d.index);
                }
                for (size_t i = d.index + 1; i < spc.groups.size(); i++) {
                    cache(d.index, i) = other->cache(d.index, i);
                }
            }
        }
    }
};

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
    void to_json(json &) const override;
    void addEwald(const json &, Space &); //!< Adds an instance of reciprocal space Ewald energies (if appropriate)
    void force(PointVector &) override;
  public:
    Hamiltonian(Space &spc, const json &j);
    double energy(Change &change) override; //!< Energy due to changes
    void init() override;
    void sync(Energybase *basePtr, Change &change) override;
}; //!< Aggregates and sum energy terms

} // namespace Energy
} // namespace Faunus
