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
#include <algorithm>

#ifdef ENABLE_FREESASA
#include <freesasa.h>
#endif

#if defined(__cpp_lib_parallel_algorithm) && __has_include(<tbb/tbb.h>)
#include <execution>
#endif

namespace Faunus {

namespace ReactionCoordinate {
class ReactionCoordinateBase;
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
     * to simplistic search which scales as number_of_bonds x number_of_moved_particles
     */
    template <class RangeOfIndex>
    double sum_energy(const Bonded::BondVector &bonds, const RangeOfIndex &indices_of_particles) const {
        assert(std::is_sorted(indices_of_particles.begin(), indices_of_particles.end()));

        auto bond_filter = [&](const auto &bond_ptr) { // determine if bond is part of indices of particles
            for (auto index : bond_ptr->index) {
                if (std::binary_search(indices_of_particles.begin(), indices_of_particles.end(), index)) {
                    return true;
                }
            }
            return false;
        };
        auto affected_bonds = bonds | ranges::cpp20::views::filter(bond_filter);

        auto bond_energy = [&](const auto &bond_ptr) { return bond_ptr->energyFunc(spc.geo.getDistanceFunc()); };

#if (defined(__clang__) && __clang_major__ >= 10) || (defined(__GNUC__) && __GNUC__ >= 10)
        return std::transform_reduce(affected_bonds.begin(), affected_bonds.end(), 0.0, std::plus<>(), bond_energy);
#else
        double energy = 0.0;
        for (const auto &bond_ptr : affected_bonds) {
            energy += bond_energy(bond_ptr);
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
 * @brief A basic accumulator which immediately computes and adds energy of a pair of particles upon addition using
 * the PairEnergy templated class.
 *
 * Generally this is the original way how the pairwise nonbonded energy has been computed in Faunus. Due to compiler
 * optimization, templated class method 'PairEnergy.potential' may be inlined to significantly improve performance.
 *
 * @tparam PairEnergy  pair energy implementing a potential(a, b) method for particles a and b
 */
template <typename PairEnergy>
struct BasicEnergyAccumulator {
  protected:
    PairEnergy &pair_energy; //!< recipe to compute non-bonded energy between two particles, see PairEnergy
    double value = 0.0;      //!< accumulated energy

  public:
    typedef const std::reference_wrapper<const Space::Tparticle> ParticleRef;

    BasicEnergyAccumulator(PairEnergy &pair_energy, const double value = 0.0) : pair_energy(pair_energy), value(value) {}

    inline BasicEnergyAccumulator &operator=(const double new_value) {
        value = new_value;
        return *this;
    }

    inline BasicEnergyAccumulator &operator+=(const double new_value) {
        value += new_value;
        return *this;
    }

    inline BasicEnergyAccumulator &operator+=(const std::pair<ParticleRef, ParticleRef> &pair) {
        // keep this short to get inlined
        value += pair_energy.potential(pair.first.get(), pair.second.get());
        return *this;
    }

    template <typename TOtherAccumulator>
    inline BasicEnergyAccumulator &operator+=(const TOtherAccumulator &acc) {
        value += static_cast<double>(acc);
        return *this;
    }

    inline explicit operator double() const { return value; }
};


/**
 * @brief Determines if two groups are separated beyond the cutoff distance.
 *
 * The distance between centers of mass is considered. The cutoff distance can be specified independently for each
 * group pair to override the default value.
 *
 * @see GroupPairingPolicy
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
    inline bool cut(const Space::Tgroup &group1, const Space::Tgroup &group2) {
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
    template <typename... Args>
    inline auto operator()(Args &&... args) { return cut(std::forward<Args>(args)...); }

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
    Space &spc;                                //!< space to init ParticleSelfEnergy with addPairPotentialSelfEnergy
    BasePointerVector<Energybase> &potentials; //!< registered non-bonded potentials, see addPairPotentialSelfEnergy
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
    template <typename T>
    inline double potential(const T &a, const T &b) const {
        assert(&a != &b); // a and b cannot be the same particle
        if constexpr (allow_anisotropic_pair_potential) {
            const Point r = geometry.vdist(a.pos, b.pos);
            return pair_potential(a, b, r.squaredNorm(), r);
        } else {
            return pair_potential(a, b, geometry.sqdist(a.pos, b.pos), {0, 0, 0});
        }
    }

    // just a temporary placement until PairForce class template will be implemented
    template <typename T>
    inline Point force(const T &a, const T &b) const {
        assert(&a != &b); // a and b cannot be the same particle
        const Point r = geometry.vdist(a.pos, b.pos);
        return pair_potential.force(a, b, r.squaredNorm(), r);
    }

    /**
     * @brief A functor alias for potential().
     * @see potential()
     */
    template <typename... Args>
    inline auto operator()(Args &&... args) {
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
 * @brief Particle pairing to calculate pairẃise interaction using particles' groups internally. Depending on
 * the accumulator provided, raw particle pairs, energy sum, etc. can be obtained.
 *
 * Accumulator is used as the first argument in all methods. Accumulator shall overload '+=' operator to accept a pair
 * of particle references as used in particle2particle method.
 *
 * @remark Method arguments are generally not checked for correctness because of performance reasons.
 *
 * @tparam TCutoff  a cutoff scheme between groups
 * @see BasicEnergyAccumulator, GroupCutoff
 */
template <typename TCutoff>
class GroupPairingPolicy {
  protected:
    const Space &spc; //!< a space to operate on
    TCutoff cut;      //!< a cutoff functor that determines if energy between two groups can be ignored

  public:
    /**
     * @param spc
     */
    GroupPairingPolicy(Space &spc)
        : spc(spc), cut(spc.geo) {}

    void from_json(const json &j) {
        Energy::from_json(j, cut);
    }

    void to_json(json &j) const {
        Energy::to_json(j, cut);
    }

    /**
     * @brief Add two interacting particles to the accumulator.
     *
     * Due to compiler optimization, the '+=' operator and this function itself may be inlined to significantly
     * improve performance.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam T  an interacting particle
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param a  first particle
     * @param b  second particle
     */
    template <typename TAccumulator, typename T>
    inline void particle2particle(TAccumulator &pair_accumulator, const T &a, const T &b) const {
        pair_accumulator += {std::cref(a), std::cref(b)};
    }

    /**
     * @brief All pairings within a group.
     *
     * All pair interaction within the group are accumulated. The pair exclusions defined in the molecule
     * topology are honoured.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param group
     * @param pair_accumulator  accumulator of interacting pairs of particles
     */
    template <typename TAccumulator, typename TGroup>
    void groupInternal(TAccumulator &pair_accumulator, const TGroup &group) {
        const auto &moldata = group.traits();
        if (!moldata.rigid) {
            const int group_size = group.size();
            for (int i = 0; i < group_size - 1; ++i) {
                for (int j = i + 1; j < group_size; ++j) {
                    // This compound condition is faster than an outer atomic condition;
                    // tested on bulk example in GCC 9.2.
                    if (group.atomic || !moldata.isPairExcluded(i, j)) {
                        particle2particle(pair_accumulator, group[i], group[j]);
                    }
                }
            }
        }
    }

    /**
     * @brief Pairings of a single particle within the group.
     *
     * The pair exclusions defined in the molecule topology are honoured.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param index  internal index of the selected particle within the group
     */
    template <typename TAccumulator, typename TGroup>
    void groupInternal(TAccumulator &pair_accumulator, const TGroup &group, const int index) {
        const auto &moldata = group.traits();
        if (!moldata.rigid) {
            if (group.atomic) {
                // speed optimization: non-bonded interaction exclusions do not need to be checked for atomic groups
                for (int i = 0; i < index; ++i) {
                    particle2particle(pair_accumulator, group[index], group[i]);
                }
                for (int i = index + 1; i < group.size(); ++i) {
                    particle2particle(pair_accumulator, group[index], group[i]);
                }
            } else {
                // molecular group
                for (int i = 0; i < index; ++i) {
                    if (!moldata.isPairExcluded(index, i)) {
                        particle2particle(pair_accumulator, group[index], group[i]);
                    }
                }
                for (int i = index + 1; i < group.size(); ++i) {
                    if (!moldata.isPairExcluded(index, i)) {
                        particle2particle(pair_accumulator, group[index], group[i]);
                    }
                }
            }
        }
    }

    /**
     * @brief Pairing in the group involving only the particles present in the index.
     *
     * Only such non-bonded pair interactions within the group are considered if at least one particle is present
     * in the index. The pair exclusions defined in the molecule topology are honoured.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @tparam TIndex
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param index  internal indices of particles within the group
     */
    template <typename TAccumulator, typename TGroup, typename TIndex>
    void groupInternal(TAccumulator &pair_accumulator, const TGroup &group, const TIndex &index) {
        auto &moldata = group.traits();
        if (!moldata.rigid) {
            if (index.size() == 1) {
                groupInternal(pair_accumulator, group, index[0]);
            } else {
                // TODO investigate overhead of `index_complement` filtering;
                // TODO perhaps allow different strategies based on the index-size/group-size ratio
                auto index_complement = indexComplement(group.size(), index);
                // moved <-> static
                for (int i : index) {
                    for (int j : index_complement) {
                        if (!moldata.isPairExcluded(i, j)) {
                            particle2particle(pair_accumulator, group[i], group[j]);
                        }
                    }
                }
                // moved <-> moved
                for (auto i_it = index.begin(); i_it < index.end(); ++i_it) {
                    for (auto j_it = std::next(i_it); j_it < index.end(); ++j_it) {
                        if (!moldata.isPairExcluded(*i_it, *j_it)) {
                            particle2particle(pair_accumulator, group[*i_it], group[*j_it]);
                        }
                    }
                }
            }
        }
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
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group1
     * @param group2
     */
    template <typename TAccumulator, typename TGroup>
    void group2group(TAccumulator &pair_accumulator, const TGroup &group1, const TGroup &group2) {
        if (!cut(group1, group2)) {
            for (auto &particle1 : group1) {
                for (auto &particle2 : group2) {
                    particle2particle(pair_accumulator, particle1, particle2);
                }
            }
        }
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

     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group1
     * @param group2
     * @param index1  list of particle indices in group1 relative to the group beginning
     */
    template <typename TAccumulator, typename TGroup>
    void group2group(TAccumulator &pair_accumulator, const TGroup &group1, const TGroup &group2,
                     const std::vector<int> &index1) {
        if (!cut(group1, group2)) {
            for (auto particle1_ndx : index1) {
                for (auto &particle2 : group2) {
                    particle2particle(pair_accumulator, *(group1.begin() + particle1_ndx), particle2);
                }
            }
        }
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
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group1
     * @param group2
     * @param index1  list of particle indices in group1 relative to the group beginning
     * @param index2  list of particle indices in group2 relative to the group beginning
     */
    template <typename TAccumulator, typename TGroup>
    void group2group(TAccumulator &pair_accumulator, const TGroup &group1, const TGroup &group2,
                     const std::vector<int> &index1, const std::vector<int> &index2) {
        if (!cut(group1, group2)) {
            if (!index2.empty()) {
                // (∁⊕group1 × ⊕group2) + (⊕group1 × ⊕group2) = group1 × ⊕group2
                group2group(pair_accumulator, group2, group1, index2);
                // + (⊕group1 × ∁⊕group2)
                auto index2_complement = indexComplement(group2.size(), index2);
                for (auto particle1_ndx : index1) {
                    for (auto particle2_ndx : index2_complement) {
                        particle2particle(pair_accumulator, group2[particle2_ndx], group1[particle1_ndx]);
                    }
                }
            } else if (!index1.empty()) {
                // (⊕group1 × ∁⊕group2) + (⊕group1 × ⊕group2) = ⊕group1 × group2
                group2group(pair_accumulator, group1, group2, index1);
                // + (∁⊕group1 × ⊕group2) = Ø as ⊕group2 is empty
            } else {
                // both indices empty hence nothing to do
            }
        }
    }

    /**
     * @brief Complete cartesian pairing between particles in a group and a union of groups.
     *
     * group × (∪ groups)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing between
     * them is skipped. The internal energy of the group is not computed even if the group is also present in the union
     * of groups.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @tparam TGroups
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param groups
     */
    template <typename TAccumulator, typename TGroup, typename TGroups>
    void group2groups(TAccumulator &pair_accumulator, const TGroup &group, const TGroups &groups) {
        for (auto &other_group : groups) {
            if (&other_group != &group) {
                group2group(pair_accumulator, group, other_group);
            }
        }
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
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @tparam TGroups
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param group_index  groups as indices in Space::groups
     * @param index  list of particle indices in the group relative to the group beginning
     */
    template <typename TAccumulator, typename TGroup, typename TGroups>
    void group2groups(TAccumulator &pair_accumulator, const TGroup &group, const TGroups &group_index,
                      const std::vector<int> &index) {
        for (auto other_group_ndx : group_index) {
            const auto &other_group = spc.groups[other_group_ndx];
            if (&other_group != &group) {
                group2group(pair_accumulator, group, other_group, index);
            }
        }
    }

    /**
     * @brief Complete cartesian pairing between particles in a group and particles in other groups in space.
     *
     * group × (space ∖ group)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     */
    template <typename TAccumulator, typename Tgroup>
    void group2all(TAccumulator &pair_accumulator, const Tgroup &group) {
        for (auto &other_group : spc.groups) {
            if (&other_group != &group) {
                group2group(pair_accumulator, group, other_group);
            }
        }
    }

    /**
     * @brief Complete cartesian pairing between a single particle in a group and particles in other groups in space.
     *
     * ⊕group × (space ∖ group), where ⊕ denotes a filter by an index (here a single particle)
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped. This method is performance-optimized version of the multiple indices method.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TGroup
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param index  a particle index relative to the group beginning
     */
    template <typename TAccumulator, typename TGroup>
    void group2all(TAccumulator &pair_accumulator, const TGroup &group, const int index) {
        const auto &particle = group[index];
        for (auto &other_group : spc.groups) {
            if (&other_group != &group) {                      // avoid self-interaction
                if (!cut(other_group, group)) {                // check g2g cut-off
                    for (auto &other_particle : other_group) { // loop over particles in other group
                        particle2particle(pair_accumulator, particle, other_particle);
                    }
                }
            }
        }
    }

    /**
     * @brief Complete cartesian pairing between selected particles in a group and particles in other groups in space.
     *
     * ⊕group × (space ∖ group), where ⊕ denotes a filter by an index
     *
     * If the distance between the groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group
     * @param index  list of particle indices in the group relative to the group beginning
     */
    template <typename TAccumulator, typename Tgroup>
    void group2all(TAccumulator &pair_accumulator, const Tgroup &group, const std::vector<int> &index) {
        if (index.size() == 1) {
            group2all(pair_accumulator, group, index[0]);
        } else {
            for (auto &other_group : spc.groups) {
                if (&other_group != &group) {
                    group2group(pair_accumulator, group, other_group, index);
                }
            }
        }
    }

    /**
     * @brief Cross pairing of particles among a union of groups. No internal pairs within any group are considered.
     *
     * If the distance between any two groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam T
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group_index  list of groups
     */
    template <typename TAccumulator, typename T>
    void groups2self(TAccumulator &pair_accumulator, const T &group_index) {
        for (auto group1_ndx_it = group_index.begin(); group1_ndx_it < group_index.end(); ++group1_ndx_it) {
            //no such move exists that the internal energy has to be recalculated
            //groupInternal(pair_accumulator, spc.groups[*group1_ndx_it]);
            for (auto group2_ndx_it = std::next(group1_ndx_it); group2_ndx_it < group_index.end(); group2_ndx_it++) {
                group2group(pair_accumulator, spc.groups[*group1_ndx_it], spc.groups[*group2_ndx_it]);
            }
        }
    }

    /**
     * @brief Cross pairing of particles between a union of groups and its complement in space.
     *
     * If the distance between any two groups is greater or equal to the group cutoff distance, the particle pairing
     * between them is skipped.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam T
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param group_index  list of groups
     */
    template <typename TAccumulator, typename T>
    void groups2all(TAccumulator &pair_accumulator, const T &group_index) {
        groups2self(pair_accumulator, group_index);
        auto index_complement = indexComplement(spc.groups.size(), group_index);
        for (auto group1_ndx : group_index) {
            for (auto group2_ndx : index_complement) {
                group2group(pair_accumulator, spc.groups[group1_ndx], spc.groups[group2_ndx]);
            }
        }
    }

    /**
     * @brief Cross pairing between all particles in the space.
     *
     * If the distance between particles' groups is greater or equal to the group cutoff distance, no calculation is
     * performed.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @param pair_accumulator  accumulator of interacting pairs of particles
     */
    template <typename TAccumulator>
    void all(TAccumulator &pair_accumulator) {
        for (auto group_it = spc.groups.begin(); group_it < spc.groups.end(); ++group_it) {
            groupInternal(pair_accumulator, *group_it);
            for (auto other_group_it = std::next(group_it); other_group_it < spc.groups.end(); other_group_it++) {
                group2group(pair_accumulator, *group_it, *other_group_it);
            }
        }
    }

    /**
     * @brief Cross pairing between all particles in the space.
     *
     * If the distance between particles' groups is greater or equal to the group cutoff distance, no calculation is
     * performed.
     *
     * @tparam TAccumulator  an accumulator with '+=' operator overloaded to add a pair of particles as references
     *                       {T&, T&}
     * @tparam TCondition  a function returning bool and having a group as an argument
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param condition  a group filter if internal energy of the group shall be added
     */
    template <typename TAccumulator, typename TCondition>
    void all(TAccumulator &pair_accumulator, TCondition condition) {
        for (auto group_it = spc.groups.begin(); group_it < spc.groups.end(); ++group_it) {
            if (condition(*group_it)) {
                 groupInternal(pair_accumulator, *group_it);
            }
            for (auto other_group_it = std::next(group_it); other_group_it < spc.groups.end(); other_group_it++) {
                 group2group(pair_accumulator, *group_it, *other_group_it);
            }
        }
    }
};

/**
 * @brief Computes pair quantity difference for a systen perturbation. Such quantity can be energy using nonponded
 * pair potential
 * .
 * @tparam TPolicy  a pairing policy
 */
template <typename TPolicy>
class GroupPairing {
    const Space &spc;
    TPolicy pairing;

  protected:
    /**
     * @brief Computes pair quantity difference if only a single group has changed.
     *
     * @tparam TAccumulator
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param change
     */
    template <typename TAccumulator>
    void accumulateGroup(TAccumulator &pair_accumulator, const Change &change) {
        const auto &change_data = change.groups.at(0);
        const auto &group = spc.groups.at(change_data.index);
        if (change_data.atoms.size() == 1) {
            // faster algorithm if only a single particle moves
            pairing.group2all(pair_accumulator, group, change_data.atoms[0]);
            if (change_data.internal) {
                pairing.groupInternal(pair_accumulator, group, change_data.atoms[0]);
            }
        } else {
            const bool change_all = change_data.atoms.empty(); // all particles or only their subset?
            if (change_all) {
                pairing.group2all(pair_accumulator, group);
                if (change_data.internal) {
                    pairing.groupInternal(pair_accumulator, group);
                }
            } else {
                pairing.group2all(pair_accumulator, group, change_data.atoms);
                if (change_data.internal) {
                    pairing.groupInternal(pair_accumulator, group, change_data.atoms);
                }
            }
        }
    }

    /**
     * @brief Computes pair quantity difference if the number of particles has changed.
     *
     * Particles have to be explicitly enumerated in the atom indices of the changed group. Implicit addition of atoms
     * with a group is not supported yet. Note that we do not have to care about missing (removed) particles at all.
     * They are taken into account in the original (old) space where they are present.
     *
     * @param pair_accumulator  accumulator of interacting pairs of particles
     * @param change
     */
    template <typename TAccumulator>
    void accumulateSpeciation(TAccumulator &pair_accumulator, const Change &change) {
        assert(change.dN);
        const auto &moved = change.touchedGroupIndex(); // index of moved groups
        const auto &fixed = indexComplement(int(spc.groups.size()), moved) | ranges::to<std::vector>; // index of static groups
        auto filter_active = [](int size) { return ranges::views::filter([size](const auto i) { return i < size; }); };

        // loop over all changed groups
        for (auto change_group1_it = change.groups.begin(); change_group1_it < change.groups.end(); ++change_group1_it) {
            auto &group1 = spc.groups.at(change_group1_it->index);
            // filter only active particles
            const std::vector<int> index1 = change_group1_it->atoms | filter_active(group1.size()) | ranges::to<std::vector>;
            if (!index1.empty()) {
                // particles added into the group: compute (changed group) <-> (static group)
                pairing.group2groups(pair_accumulator, group1, fixed, index1);
            }
            // loop over successor changed groups (hence avoid double counting group1×group2 and group2×group1)
            for (auto change_group2_it = std::next(change_group1_it); change_group2_it < change.groups.end(); ++change_group2_it) {
                auto &group2 = spc.groups.at(change_group2_it->index);
                const std::vector<int> index2 = change_group2_it->atoms | filter_active(group2.size()) | ranges::to<std::vector>;
                if (!index1.empty() || !index2.empty()) {
                    // particles added into one or other group: compute (changed group) <-> (changed group)
                    pairing.group2group(pair_accumulator, group1, group2, index1, index2);
                }
            }
            if (!index1.empty() && !molecules.at(group1.id).rigid) {
                // compute internal energy in the changed group
                if (change_group1_it->all) {
                    pairing.groupInternal(pair_accumulator, group1);
                } else {
                    pairing.groupInternal(pair_accumulator, group1, index1);
                };
            }
        }
    }

  public:
    /**
     * @brief Computes pair quantity difference from changed particles.
     *
     * The internal energy contribution, i.e., the contribution from the intra group interactions, is added
     * only if a single group is changed or if all changed.
     *
     * @param change
     * @param pair_accumulator  accumulator of interacting pairs of particles
     */
    template <typename TAccumulator>
    void accumulate(TAccumulator &pair_accumulator, const Change &change) {
        assert(std::is_sorted(change.groups.begin(), change.groups.end()));
        if (change.all) {
            pairing.all(pair_accumulator);
        } else if (change.dV) {
            // sum all interaction energies except the internal energies of incompressible molecules
            pairing.all(pair_accumulator, [](auto &group) { return group.atomic || group.compressible; });
        } else if (!change.dN) {
            if (change.groups.size() == 1) {
                // if only a single group changes use faster algorithm and optionally add the internal energy
                accumulateGroup(pair_accumulator, change);
            } else {
                // if multiple groups move, no internal energies are computed
                const auto &moved = change.touchedGroupIndex(); // index of moved groups
                pairing.groups2all(pair_accumulator, moved);
            }
        } else { // change.dN
            accumulateSpeciation(pair_accumulator, change);
        }
    }

    GroupPairing(Space &spc) : spc(spc), pairing(spc) {}

    void from_json(const json &j) {
        pairing.from_json(j);
    }

    void to_json(json &j) const {
        pairing.to_json(j);
    }

    // FIXME a temporal fix for non-refactorized NonbondedCached
    template <typename Accumulator>
    void group2group(Accumulator &pair_accumulator, const Space::Tgroup &group1, const Space::Tgroup &group2) {
        pairing.group2group(std::forward<Accumulator &>(pair_accumulator), std::forward<const Space::Tgroup &>(group1),
                            std::forward<const Space::Tgroup &>(group2));
    }
};

/**
 * @brief Computes change in the non-bonded energy, assuming pairwise additive energy terms.
 *
 * @tparam TPairEnergy  a functor to compute non-bonded energy between two particles
 * @tparam TPairingPolicy  pairing policy to effectively sum up the pairwise additive non-bonded energy
 */
template <typename TPairEnergy, typename TPairingPolicy>
class Nonbonded : public Energybase {
  protected:
    typedef BasicEnergyAccumulator<TPairEnergy> TAccumulator;
    const Space &spc;              //!< space to operate on
    TPairEnergy pair_energy; //!< a functor to compute non-bonded energy between two particles, see PairEnergy
    TPairingPolicy pairing;  //!< pairing policy to effectively sum up the pairwise additive non-bonded energy

  public:
    Nonbonded(const json &j, Space &spc, BasePointerVector<Energybase> &pot)
        : spc(spc), pair_energy(spc, pot), pairing(spc) {
        name = "nonbonded";
        from_json(j);
    }

    void from_json(const json &j) {
        pair_energy.from_json(j);
        pairing.from_json(j);
    }

    void to_json(json &j) const override {
        pair_energy.to_json(j);
        pairing.to_json(j);
    }

     double energy(Change &change) override {
        TAccumulator energy_accumulator(pair_energy, 0.0);
        pairing.accumulate(energy_accumulator, change);
        return static_cast<double>(energy_accumulator);
    }

    /**
     * @brief Calculates the force on all particles.
     *
     * @todo A stub. Change to reflect only active particle, see Space::activeParticles().
     */
    void force(std::vector<Point> &forces) override {
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


/**
 * @brief Computes non-bonded energy contribution from changed particles. Cache group2group energy once calculated,
 * until a new trial configuration is provided. Not for general use as only partially implemented!
 *
 * Original implementation, only refurbished. Generally suboptimal as only PairingPolicy::group2group method
 * may be called.
 * No internal energy is ever computed. Cannot deal with particle count changes. And other unmentioned constrains.
 *
 * @tparam TPairEnergy  a functor to compute non-bonded energy between two particles
 * @tparam TPairingPolicy  pairing policy to effectively sum up the pairwise additive non-bonded energy
 */
template <typename TPairEnergy, typename TPairingPolicy>
class NonbondedCached : public Nonbonded<TPairEnergy, TPairingPolicy> {
    typedef Nonbonded<TPairEnergy, TPairingPolicy> Base;
    typedef BasicEnergyAccumulator<TPairEnergy> TAccumulator;
    Eigen::MatrixXf energy_cache;
    using Base::spc;

    template <typename TGroup>
    double g2g(const TGroup &g1, const TGroup &g2) {
        int i = &g1 - spc.groups.data();
        int j = &g2 - spc.groups.data();
        if (j < i) {
            std::swap(i, j);
        }
        if (Energybase::key == Energybase::TRIAL_MONTE_CARLO_STATE) { // if this is from the trial system
            TAccumulator energy_accumulator(Base::pair_energy);
            Base::pairing.group2group(energy_accumulator, g1, g2);
            energy_cache(i, j) = static_cast<double>(energy_accumulator);  // update the cache
        }
        return energy_cache(i, j); // return (cached) value
    }

    template <typename TGroup>
    double g2g(const TGroup &g1, const TGroup &g2, [[maybe_unused]] const std::vector<int> &index) {
        // index not implemented
        return g2g(g1, g2);
    }

  public:
    NonbondedCached(const json &j, Space &spc, BasePointerVector<Energybase> &pot) : Base(j, spc, pot) {
        Base::name += "EM";
        init();
    }

    /**
     * @brief Cache pair interactions in matrix.
     */
    void init() override {
        const auto groups_size = spc.groups.size();
        energy_cache.resize(groups_size, groups_size);
        energy_cache.setZero();
        TAccumulator u(Base::pair_energy);
        for (auto i = 0; i < groups_size - 1; ++i) {
            for (auto j = i + 1; j < groups_size; ++j) {
                u = 0.0;
                Base::pairing.group2group(u, spc.groups[i], spc.groups[j]);
                energy_cache(i, j) = static_cast<double>(u);
            }
        }
    }

    double energy(Change &change) override {
        // Only g2g may be called there to compute (and cache) energy!
        double energy_sum = 0.0;
        if (change) {
            if (change.all || change.dV) {
                for (auto i = spc.groups.begin(); i < spc.groups.end(); ++i) {
                    for (auto j = std::next(i); j < Base::spc.groups.end(); ++j) {
                        energy_sum += g2g(*i, *j);
                    }
                }
            } else {
                if (change.groups.size() == 1) { // if exactly ONE molecule is changed
                    auto &d = change.groups[0];
                    auto &g1 = spc.groups.at(d.index);
                    for (auto g2_it = spc.groups.begin(); g2_it < spc.groups.end(); ++g2_it) {
                        if (&g1 != &(*g2_it)) {
                            energy_sum += g2g(g1, *g2_it, d.atoms);
                        }
                    }
                } else {                                     // many molecules are changed
                    auto moved = change.touchedGroupIndex(); // index of moved groups
                    // moved<->moved
                    if (change.moved2moved) {
                        for (auto i = moved.begin(); i < moved.end(); ++i) {
                            for (auto j = std::next(i); j < moved.end(); ++j) {
                                energy_sum += g2g(spc.groups[*i], spc.groups[*j]);
                            }
                        }
                    }
                    // moved<->static
#if true
                    // classic version
                    auto fixed = indexComplement(spc.groups.size(), moved); // index of static groups
                    for (auto i : moved) {
                        for (auto j : fixed) {
                            energy_sum += g2g(spc.groups[i], spc.groups[j]);
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
                            energy_sum += g2g(spc.groups[moved[i]], spc.groups[fixed[j]]);
                        }
                    }
#endif
                }
            }
            // more todo!
        }
        return energy_sum;
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
            energy_cache.triangularView<Eigen::StrictlyUpper>() =
                (other->energy_cache).template triangularView<Eigen::StrictlyUpper>();
        } else {
            for (auto &d : change.groups) {
                for (int i = 0; i < d.index; i++) {
                    energy_cache(i, d.index) = other->energy_cache(i, d.index);
                }
                for (size_t i = d.index + 1; i < spc.groups.size(); i++) {
                    energy_cache(d.index, i) = other->energy_cache(d.index, i);
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

/**
 * @brief Oscillating energy on a single particle
 *
 * This is 2D version of the oscillating potential used
 * to illustrate parallel tempering in the book
 * "Understanding Molecular Simulation" by D. Frenkel.
 */
class Example2D : public Energybase {
  private:
    bool use_2d = true;        // Set to false to apply energy only along x (as by the book)
    double scale_energy = 1.0; // effective temperature
    const Point &particle;     // reference to 1st particle in the system
    void to_json(json &j) const override;

  public:
    Example2D(const json &j, Space &spc);
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
