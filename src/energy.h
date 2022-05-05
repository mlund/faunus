#pragma once

#include "bonds.h"
#include "externalpotential.h" // Energybase implemented here
#include "potentials_base.h"
#include "sasa.h"
#include "space.h"
#include "aux/iteratorsupport.h"
#include "aux/pairmatrix.h"
#include "smart_montecarlo.h"
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>
#include <numeric>
#include <algorithm>
#include <concepts>

struct freesasa_parameters_fwd; // workaround for freesasa unnamed struct that cannot be forward declared

#if defined(__cpp_lib_parallel_algorithm) && __has_include(<tbb/tbb.h>)
#include <execution>
#endif

#if defined(__cpp_lib_parallel_algorithm) &&                                                                           \
    __has_include(<tbb/tbb.h>) && ((defined(__clang__) && __clang_major__ >= 10) || (defined(__GNUC__) && __GNUC__ >= 10))
#define HAS_PARALLEL_TRANSFORM_REDUCE
#endif

namespace Faunus {

namespace ReactionCoordinate {
class ReactionCoordinateBase;
}

namespace Potential {
class PairPotentialBase;
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
 * @brief Check if particles are outside the simulation container
 *
 * If any particles is ouside, infinite energy is returned; zero otherwirse.
 * This is not needed for cuboidal geometry as particles are always wrapped using PBC.
 */
class ContainerOverlap : public Energybase {
  private:
    const Space& spc;
    bool groupIsOutsideContainer(const Change::GroupChange& group_change) const;
    double energyOfAllGroups() const;

  public:
    explicit ContainerOverlap(const Space& spc);
    double energy(const Change& change) override;
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
    using Tcomplex = std::complex<double>;
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
    enum Policies { PBC, PBCEigen, IPBC, IPBCEigen, METALSLIT, INVALID }; //!< Possible k-space updating schemes
    Policies policy = PBC;                                     //!< Policy for updating k-space
    explicit EwaldData(const json& j);                         //!< Initialize from json
};

NLOHMANN_JSON_SERIALIZE_ENUM(EwaldData::Policies, {
                                                      {EwaldData::INVALID, nullptr},
                                                      {EwaldData::PBC, "PBC"},
                                                      {EwaldData::PBCEigen, "PBCEigen"},
                                                      {EwaldData::IPBC, "IPBC"},
                                                      {EwaldData::METALSLIT, "METALSLIT"},
                                                      {EwaldData::IPBCEigen, "IPBCEigen"},
                                                  })

void to_json(json& j, const EwaldData& d);

/**
 * @brief Base class for Ewald k-space updates policies
 */
class EwaldPolicyBase {
  public:
    std::string cite; //!< Optional reference, preferably DOI, to further information
    virtual ~EwaldPolicyBase() = default;
    virtual void updateBox(EwaldData&, const Point&) const = 0; //!< Prepare k-vectors according to given box vector
    virtual void updateComplex(EwaldData& d,
                               const Space::GroupVector& groups) const = 0; //!< Update all k vectors
    virtual void
    updateComplex(EwaldData& d, const Change& change, const Space::GroupVector& groups,
                  const Space::GroupVector& oldgroups) const = 0; //!< Update subset of k vectors. Require `old` pointer
    virtual double selfEnergy(const EwaldData& d, Change& change,
                              Space::GroupVector& groups) = 0; //!< Self energy contribution due to a change
    virtual double surfaceEnergy(const EwaldData& d, const Change& change,
                                 const Space::GroupVector& groups) = 0; //!< Surface energy contribution due to a change
    virtual double reciprocalEnergy(const EwaldData& d) = 0;            //!< Total reciprocal energy

    /**
     * @brief Represent charges and positions using an Eigen facade (Map)
     *
     * Requires that all groups are fully active, i.e. does not work for GCMC.
     *
     * @param groups Vector of groups to represent
     * @return tuple with positions, charges
     */
    static auto mapGroupsToEigen(Space::GroupVector& groups) {
        auto is_partially_inactive = [](const Group& group) { return group.size() != group.capacity(); };
        if (ranges::cpp20::any_of(groups, is_partially_inactive)) {
            throw std::runtime_error("Eigen optimized Ewald not available with inactive groups");
        }
        auto first_particle = groups.front().begin();
        auto last_particle = groups.back().end();
        auto pos = asEigenMatrix(first_particle, last_particle,
                                 &Particle::pos); // N x 3
        auto charge = asEigenVector(first_particle, last_particle,
                                    &Particle::charge); // N x 1
        return std::make_tuple(pos, charge);
    }

    static auto mapGroupsToEigen(const Space::GroupVector& groups) {
        auto is_partially_inactive = [](const Group& group) { return group.size() != group.capacity(); };
        if (ranges::cpp20::any_of(groups, is_partially_inactive)) {
            throw std::runtime_error("Eigen optimized Ewald not available with inactive groups");
        }
        auto first_particle = groups.front().begin();
        auto last_particle = groups.back().end();
        auto pos = asEigenMatrix(first_particle, last_particle,
                                 &Particle::pos); // N x 3
        auto charge = asEigenVector(first_particle, last_particle,
                                    &Particle::charge); // N x 1
        return std::make_tuple(pos, charge);
    }

    static std::unique_ptr<EwaldPolicyBase> makePolicy(EwaldData::Policies); //!< Policy factory
};

/**
 * @brief Ion-Ion Ewald using periodic boundary conditions (PBC)
 */
struct PolicyIonIon : public EwaldPolicyBase {
    PolicyIonIon();
    void updateBox(EwaldData& d, const Point& box) const override;
    void updateComplex(EwaldData& data, const Space::GroupVector& groups) const override;
    void updateComplex(EwaldData& d, const Change& change, const Space::GroupVector& groups,
                       const Space::GroupVector& oldgroups) const override;
    double selfEnergy(const EwaldData& d, Change& change, Space::GroupVector& groups) override;
    double surfaceEnergy(const EwaldData& data, const Change& change, const Space::GroupVector& groups) override;
    double reciprocalEnergy(const EwaldData& d) override;
};

/**
 * Add documentation here...
 * @todo register image particles when updating k-vectors
 */
struct PolicyIonIonMetalSlit : public PolicyIonIon {
    PolicyIonIonMetalSlit();
    void updateBox(EwaldData& d, const Point& box) const override;
    void updateComplex(EwaldData& data, const Space::GroupVector& groups) const override;
    void updateComplex(EwaldData& d, const Change& change, const Space::GroupVector& groups,
                       const Space::GroupVector& oldgroups) const override;
    double selfEnergy(const EwaldData& d, Change& change, Space::GroupVector& groups) override;
    double surfaceEnergy(const EwaldData& data, const Change& change, const Space::GroupVector& groups) override;
    double reciprocalEnergy(const EwaldData& d) override;
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
    void updateComplex(EwaldData&, const Space::GroupVector&) const override;
    double reciprocalEnergy(const EwaldData &) override;
};

/**
 * @brief Ion-Ion Ewald with isotropic periodic boundary conditions (IPBC)
 */
struct PolicyIonIonIPBC : public PolicyIonIon {
    using PolicyIonIon::updateComplex;
    PolicyIonIonIPBC();
    void updateBox(EwaldData &, const Point &) const override;
    void updateComplex(EwaldData&, const Space::GroupVector&) const override;
    void updateComplex(EwaldData&, const Change&, const Space::GroupVector&, const Space::GroupVector&) const override;
};

/**
 * @brief Ion-Ion Ewald with isotropic periodic boundary conditions (IPBC) using Eigen operations
 * @warning Incomplete and under construction
 */
struct PolicyIonIonIPBCEigen : public PolicyIonIonIPBC {
    using PolicyIonIonIPBC::updateComplex;
    void updateComplex(EwaldData&, const Space::GroupVector&) const override;
};

/**
 * @brief Ewald summation reciprocal energy
 * @todo energy() currently has the responsibility to update k-vectors.
 *       This is error prone and should be handled *before* this step.
 */
class Ewald : public Energybase {
  protected:
    const Space& spc;
    EwaldData data;
    std::shared_ptr<EwaldPolicyBase> policy; //!< Policy for updating k-space
    const Space::GroupVector* old_groups = nullptr;

  public:
    Ewald(const Space& spc, const EwaldData& data);
    Ewald(const json& j, const Space& spc);
    void init() override;
    void setOldGroups(const Space::GroupVector& old_groups); //!< Optimization if old groups are available (optional)
    void updateState(const Change& change) override;
    double energy(const Change& change) override;
    void sync(Energybase* energybase, const Change& change) override;
    void to_json(json& j) const override;
    void force(std::vector<Point>& forces) override; // update forces on all particles
};

/**
 * @todo Add documentation
 *
 * The only thing we need to do here, is to add the real <-> imaginary interaction.
 * The reciprocal part is fully captured by the base class
 */
class MetallicEwald : public Ewald {
  public:
    MetallicEwald(const json& j, const Space& spc);

    double energy(const Change& change) override;
};

/**
 * @brief Pressure term for NPT ensemble
 */
class Isobaric : public Energybase {
  private:
    const Space& spc;
    double pressure = 0.0;                                     //!< Applied pressure
    static const std::map<std::string, double> pressure_units; //!< Possible ways pressure can be given
  public:
    Isobaric(const json& j, const Space& spc);
    double energy(const Change& change) override;
    void to_json(json& j) const override;
};

/**
 * @brief Constrain system using reaction coordinates
 *
 * If outside specified `range`, infinity energy is returned, causing rejection.
 */
class Constrain : public Energybase {
  private:
    std::string type;
    std::unique_ptr<ReactionCoordinate::ReactionCoordinateBase> coordinate;

  public:
    Constrain(const json& j, Space& space);
    double energy(const Change& change) override;
    void to_json(json& j) const override;
};

/**
 * The keys of the `intra` map are group index and the values
 * is a vector of `BondData`. For bonds between groups, fill
 * in `inter` which is evaluated for every update of call to
 * `energy`.
 *
 * @todo Optimize.
 */
class Bonded : public Energybase {
  private:
    using BondVector = BasePointerVector<Potential::BondData>;
    const Space& spc;
    BondVector external_bonds;                                      //!< inter-molecular bonds
    std::map<int, BondVector> internal_bonds;                       //!< intra-molecular bonds; key is group index
    void updateGroupBonds(const Space::GroupType& group);           //!< Update/set bonds internally in group
    double sumBondEnergy(const BondVector& bonds) const;            //!< sum energy in vector of BondData
    double internalGroupEnergy(const Change::GroupChange& changed); //!< Energy from internal bonds
    double sumEnergy(const BondVector& bonds, const ranges::cpp20::range auto& particle_indices) const;
    void updateInternalBonds(); //!< finds and adds all intra-molecular bonds of active molecules

  public:
    Bonded(const Space& spc, const BondVector& external_bonds);
    Bonded(const json& j, const Space& spc);
    void to_json(json& j) const override;
    double energy(const Change& change) override;    //!< brute force -- refine this!
    void force(std::vector<Point>& forces) override; //!< Calculates the forces on all particles
};

/**
 * @brief Sum energy in vector of BondData for matching particle indices
 * @param bonds List of bonds
 * @param particle_indices Particle indices to calculate the energy for
 *
 * To speed up the bond search, the given indices must be ordered which allows
 * for binary search which on large systems provides superior performance compared
 * to simplistic search which scales as number_of_bonds x number_of_moved_particles
 */
double Bonded::sumEnergy(const Bonded::BondVector& bonds, const ranges::cpp20::range auto& particle_indices) const {
    assert(std::is_sorted(particle_indices.begin(), particle_indices.end()));

    auto index_is_included = [&](auto index) {
        return std::binary_search(particle_indices.begin(), particle_indices.end(), index);
    };
    auto affected_bonds = bonds | ranges::cpp20::views::filter([&](const auto& bond) {
                              return ranges::cpp20::any_of(bond->indices, index_is_included);
                          });
    auto bond_energy = [dist = spc.geometry.getDistanceFunc()](const auto& bond) { return bond->energyFunc(dist); };
    return std::transform_reduce(affected_bonds.begin(), affected_bonds.end(), 0.0, std::plus<>(), bond_energy);
}

/**
 * @brief Provides a complementary set of ints with respect to the iota set of a given size.
 * @remark It is used as a helper function for pair interactions.
 *
 * @param size the iota superset contains all integers in the range [0, size)
 * @param range an original set of integers (must be sorted)
 * @return a set of ints complementary to the original set
 */
auto indexComplement(std::integral auto size, const ranges::cpp20::range auto& range) {
    namespace rv = ranges::cpp20::views;
    return rv::iota(0, static_cast<int>(size)) |
           rv::filter([&](auto i) { return !std::binary_search(range.begin(), range.end(), i); });
}

/**
 * @brief Provides a fast inlineable interface for non-bonded pair potential energy computation.
 *
 * @tparam TPairPotential  a pair potential to compute with
 * @tparam allow_anisotropic_pair_potential  pass also a distance vector to the pair potential, slower
 */
template <Potential::RequirePairPotential TPairPotential, bool allow_anisotropic_pair_potential = true>
class PairEnergy {
    const Space::GeometryType& geometry;       //!< geometry to operate with
    TPairPotential pair_potential;             //!< pair potential function/functor
    Space& spc;                                //!< space to init ParticleSelfEnergy with addPairPotentialSelfEnergy
    BasePointerVector<Energybase>& potentials; //!< registered non-bonded potentials, see addPairPotentialSelfEnergy
  public:
    /**
     * @param spc
     * @param potentials  registered non-bonded potentials
     */
    PairEnergy(Space& spc, BasePointerVector<Energybase>& potentials)
        : geometry(spc.geometry)
        , spc(spc)
        , potentials(potentials) {}

    /**
     * @brief Computes pair potential energy.
     *
     * @param a  particle
     * @param b  particle
     * @return pair potential energy between particles a and b
     */
    template <typename T> inline double potential(const T& a, const T& b) const {
        assert(&a != &b); // a and b cannot be the same particle
        if constexpr (allow_anisotropic_pair_potential) {
            const Point r = geometry.vdist(a.pos, b.pos);
            return pair_potential(a, b, r.squaredNorm(), r);
        } else {
            return pair_potential(a, b, geometry.sqdist(a.pos, b.pos), {0, 0, 0});
        }
    }

    // just a temporary placement until PairForce class template will be implemented
    template <typename ParticleType> inline Point force(const ParticleType& a, const ParticleType& b) const {
        assert(&a != &b);                                       // a and b cannot be the same particle
        const Point b_towards_a = geometry.vdist(a.pos, b.pos); // vector b -> a = a - b
        return pair_potential.force(a, b, b_towards_a.squaredNorm(), b_towards_a);
    }

    /**
     * @brief A functor alias for potential().
     * @see potential()
     */
    template <typename... Args> inline auto operator()(Args&&... args) {
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

    void from_json(const json& j) {
        Potential::from_json(j, pair_potential);
        if (!pair_potential.isotropic && !allow_anisotropic_pair_potential) {
            throw std::logic_error("Only isotropic pair potentials are allowed.");
        }
        addPairPotentialSelfEnergy();
    }

    void to_json(json& j) const { pair_potential.to_json(j); }
};

template <typename T>
concept RequirePairEnergy = requires(T instance) {
    instance.potential(Particle(), Particle());
    instance.force(Particle(), Particle());
    instance.addPairPotentialSelfEnergy();
};

/**
 * @brief Interface for energy accumulators
 *
 * The energy accumulator is used to add up energies between two particles.
 * This can be done instantly (see `InstantEnergyAccumulator`) or delaying
 * the evaluation until the energy is needed (`DelayedEnergyAccumulator`).
 * The latter may be used with parallelism.
 *
 * @todo See https://www.youtube.com/watch?v=3LsRYnRDSRA for a bizarre example
 *       where a custom `struct Tpair { const Particle &first, second; };`
 *       outperforms `std::pair` due to missed compiler optimization.
 */
class EnergyAccumulatorBase {
  protected:
    double value = 0.0;                                               //!< accumulated energy
    using ParticleRef = const std::reference_wrapper<const Particle>; //!< Particle reference
    using ParticlePair = std::pair<ParticleRef, ParticleRef>;         //!< References to two particles

  public:
    enum class Scheme { SERIAL, OPENMP, PARALLEL, INVALID };
    Scheme scheme = Scheme::SERIAL;

    EnergyAccumulatorBase(double value);
    virtual ~EnergyAccumulatorBase() = default;
    virtual void reserve(size_t number_of_particles);
    virtual void clear();
    virtual void from_json(const json &j);
    virtual void to_json(json &j) const;

    virtual explicit operator double();
    virtual EnergyAccumulatorBase& operator=(double new_value) = 0;
    virtual EnergyAccumulatorBase& operator+=(double new_value) = 0;
    virtual EnergyAccumulatorBase& operator+=(ParticlePair&& pair) = 0;

    template <typename TOtherAccumulator> inline EnergyAccumulatorBase& operator+=(TOtherAccumulator& acc) {
        value += static_cast<double>(acc);
        return *this;
    }
};

NLOHMANN_JSON_SERIALIZE_ENUM(EnergyAccumulatorBase::Scheme, {{EnergyAccumulatorBase::Scheme::INVALID, nullptr},
                                                             {EnergyAccumulatorBase::Scheme::SERIAL, "serial"},
                                                             {EnergyAccumulatorBase::Scheme::OPENMP, "openmp"},
                                                             {EnergyAccumulatorBase::Scheme::PARALLEL, "parallel"}})

template <class T>
concept RequireEnergyAccumulator = std::is_base_of_v<EnergyAccumulatorBase, T>;

/**
 * @brief A basic accumulator which immediately computes and adds energy of a pair of particles upon addition using
 * the PairEnergy templated class.
 *
 * Generally this is the original way how the pairwise nonbonded energy has been computed in Faunus. Due to compiler
 * optimization, templated class method 'PairEnergy.potential' may be inlined to significantly improve performance.
 *
 * @tparam PairEnergy  pair energy implementing a potential(a, b) method for particles a and b
 */
template <RequirePairEnergy PairEnergy> class InstantEnergyAccumulator : public EnergyAccumulatorBase {
  private:
    const PairEnergy& pair_energy; //!< recipe to compute non-bonded energy between two particles, see PairEnergy

  public:
    InstantEnergyAccumulator(const PairEnergy& pair_energy, const double value = 0.0)
        : EnergyAccumulatorBase(value), pair_energy(pair_energy) {}

    inline InstantEnergyAccumulator& operator=(const double new_value) override {
        value = new_value;
        return *this;
    }

    inline InstantEnergyAccumulator& operator+=(const double new_value) override {
        value += new_value;
        return *this;
    }

    inline InstantEnergyAccumulator& operator+=(ParticlePair&& pair) override {
        // keep this short to get inlined
        value += pair_energy.potential(pair.first.get(), pair.second.get());
        return *this;
    }

    void from_json(const json &j) override {
        EnergyAccumulatorBase::from_json(j);
        if (scheme != Scheme::SERIAL) {
            faunus_logger->warn("unsupported summation scheme; falling back to 'serial'");
        }
    }
};

/**
 * Stores a vector of particle pairs and postpones the energy evaluation until
 * `operator double()` is called. Looping over the vector can be done in serial (as a fallback);
 * using OpenMP; or using C++17 parallel algorithms if available.
 */
template <RequirePairEnergy PairEnergy> class DelayedEnergyAccumulator : public EnergyAccumulatorBase {
  private:
    std::vector<ParticlePair> particle_pairs;
    const PairEnergy& pair_energy; //!< recipe to compute non-bonded energy between two particles, see PairEnergy
    const size_t max_particles_in_buffer = 10000; //!< this can be modified to suit memory requirements

  public:
    explicit DelayedEnergyAccumulator(const PairEnergy& pair_energy, const double value = 0.0)
        : EnergyAccumulatorBase(value), pair_energy(pair_energy) {}

    /** Reserve memory for (N-1)*N/2 interaction pairs */
    void reserve(size_t number_of_particles) override {
        try {
            number_of_particles = std::min(number_of_particles, max_particles_in_buffer);
            const auto number_of_pairs = (number_of_particles - 1U) * number_of_particles / 2U;
            faunus_logger->debug(fmt::format("reserving memory for {} energy pairs ({} MB)", number_of_pairs,
                                             number_of_pairs * sizeof(ParticlePair) / (1024U * 1024U)));
            particle_pairs.reserve(number_of_pairs);
        } catch (std::exception& e) {
            throw std::runtime_error(
                fmt::format("cannot allocate memory for energy pairs: {}. Use another summation policy.", e.what()));
        }
    }

    void clear() override {
        value = 0.0;
        particle_pairs.clear();
    }

    DelayedEnergyAccumulator& operator=(const double new_value) override {
        clear();
        value = new_value;
        return *this;
    }

    inline DelayedEnergyAccumulator& operator+=(const double new_value) override {
        value += new_value;
        return *this;
    }

    inline DelayedEnergyAccumulator& operator+=(ParticlePair&& pair) override {
        assert(particle_pairs.capacity() > 0);
        if (particle_pairs.size() == particle_pairs.capacity()) {
            operator double(); // sum stored pairs and reset buffer
        }
        particle_pairs.emplace_back(pair);
        return *this;
    }

    explicit operator double() override {
        switch (scheme) {
        case Scheme::OPENMP:
            value += accumulateOpenMP();
            break;
        case Scheme::PARALLEL:
            value += accumulateParallel();
            break;
        default:
            value += accumulateSerial();
        }
        particle_pairs.clear();
        return value;
    }

  private:
    double accumulateSerial() const {
        double sum = 0.0;
        for (const auto [particle1, particle2] : particle_pairs) {
            sum += pair_energy.potential(particle1.get(), particle2.get());
        }
        return sum;
    }

    double accumulateParallel() const {
#if defined(HAS_PARALLEL_TRANSFORM_REDUCE)
        return std::transform_reduce(
            std::execution::par, particle_pairs.cbegin(), particle_pairs.cend(), 0.0, std::plus<>(),
            [&](const auto& pair) { return pair_energy.potential(pair.first.get(), pair.second.get()); });
#else
        return accumulateSerial(); // fallback
#endif
    }

    double accumulateOpenMP() const {
        double sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
        for (const auto& pair : particle_pairs) {
            sum += pair_energy.potential(pair.first.get(), pair.second.get());
        }
        return sum;
    }
};

template <RequirePairEnergy TPairEnergy>
std::unique_ptr<EnergyAccumulatorBase> createEnergyAccumulator(const json& j, const TPairEnergy& pair_energy,
                                                               double initial_value) {
    std::unique_ptr<EnergyAccumulatorBase> accumulator;
    if (j.value("summation_policy", EnergyAccumulatorBase::Scheme::SERIAL) != EnergyAccumulatorBase::Scheme::SERIAL) {
        accumulator = std::make_unique<DelayedEnergyAccumulator<TPairEnergy>>(pair_energy, initial_value);
        faunus_logger->debug("activated delayed energy summation");
    } else {
        accumulator = std::make_unique<InstantEnergyAccumulator<TPairEnergy>>(pair_energy, initial_value);
        faunus_logger->debug("activated instant energy summation");
    }
    accumulator->from_json(j);
    return accumulator;
}

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
    PairMatrix<double> cutoff_squared;    //!< matrix with group-to-group cutoff distances squared in angstrom squared
    Space::GeometryType& geometry;        //!< geometry to compute the inter group distance with
    friend void from_json(const json&, GroupCutoff&);
    friend void to_json(json&, const GroupCutoff&);
    void setSingleCutoff(double cutoff);

  public:
    /**
     * @brief Determines if two groups are separated beyond the cutoff distance.
     * @return true if the group-to-group distance is beyond the cutoff distance, false otherwise
     */
    inline bool cut(const Group& group1, const Group& group2) {
        if (group1.isAtomic() || group2.isAtomic()) {
            return false; // atomic groups have ill-defined mass centers
        }
        return geometry.sqdist(group1.mass_center, group2.mass_center) >= cutoff_squared(group1.id, group2.id);
    }

    double getCutoff(size_t id1, size_t id2) const;

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
    explicit GroupCutoff(Space::GeometryType& geometry);
};

void from_json(const json&, GroupCutoff &);
void to_json(json&, const GroupCutoff &);

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
 * @see InstantEnergyAccumulator, GroupCutoff
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
    explicit GroupPairingPolicy(Space& spc)
        : spc(spc)
        , cut(spc.geometry) {}

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
    template <RequireEnergyAccumulator TAccumulator, typename T>
    inline void particle2particle(TAccumulator& pair_accumulator, const T& a, const T& b) const {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void groupInternal(TAccumulator& pair_accumulator, const TGroup& group) {
        const auto &moldata = group.traits();
        if (!moldata.rigid) {
            const int group_size = group.size();
            for (int i = 0; i < group_size - 1; ++i) {
                for (int j = i + 1; j < group_size; ++j) {
                    // This compound condition is faster than an outer atomic condition;
                    // tested on bulk example in GCC 9.2.
                    if (group.isAtomic() || !moldata.isPairExcluded(i, j)) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void groupInternal(TAccumulator& pair_accumulator, const TGroup& group, const std::size_t index) {
        const auto &moldata = group.traits();
        if (!moldata.rigid) {
            if (group.isAtomic()) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup, typename TIndex>
    void groupInternal(TAccumulator& pair_accumulator, const TGroup& group, const TIndex& index) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void group2group(TAccumulator& pair_accumulator, const TGroup& group1, const TGroup& group2) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void group2group(TAccumulator& pair_accumulator, const TGroup& group1, const TGroup& group2,
                     const std::vector<std::size_t>& index1) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void group2group(TAccumulator& pair_accumulator, const TGroup& group1, const TGroup& group2,
                     const std::vector<std::size_t>& index1, const std::vector<std::size_t>& index2) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup, typename TGroups>
    void group2groups(TAccumulator& pair_accumulator, const TGroup& group, const TGroups& groups) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup, typename TGroups>
    void group2groups(TAccumulator& pair_accumulator, const TGroup& group, const TGroups& group_index,
                      const std::vector<std::size_t>& index) {
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
    template <RequireEnergyAccumulator TAccumulator, typename Tgroup>
    void group2all(TAccumulator& pair_accumulator, const Tgroup& group) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TGroup>
    void group2all(TAccumulator& pair_accumulator, const TGroup& group, const int index) {
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
    template <RequireEnergyAccumulator TAccumulator, typename Tgroup>
    void group2all(TAccumulator& pair_accumulator, const Tgroup& group, const std::vector<std::size_t>& index) {
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
    template <RequireEnergyAccumulator TAccumulator, typename T>
    void groups2self(TAccumulator& pair_accumulator, const T& group_index) {
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
    template <RequireEnergyAccumulator TAccumulator, typename T>
    void groups2all(TAccumulator& pair_accumulator, const T& group_index) {
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
    template <RequireEnergyAccumulator TAccumulator> void all(TAccumulator& pair_accumulator) {
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
    template <RequireEnergyAccumulator TAccumulator, typename TCondition>
    void all(TAccumulator& pair_accumulator, TCondition condition) {
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
    template <RequireEnergyAccumulator TAccumulator>
    void accumulateGroup(TAccumulator& pair_accumulator, const Change& change) {
        const auto &change_data = change.groups.at(0);
        const auto& group = spc.groups.at(change_data.group_index);
        if (change_data.relative_atom_indices.size() == 1) {
            // faster algorithm if only a single particle moves
            pairing.group2all(pair_accumulator, group, change_data.relative_atom_indices[0]);
            if (change_data.internal) {
                pairing.groupInternal(pair_accumulator, group, change_data.relative_atom_indices[0]);
            }
        } else {
            const bool change_all = change_data.relative_atom_indices.empty(); // all particles or only their subset?
            if (change_all) {
                pairing.group2all(pair_accumulator, group);
                if (change_data.internal) {
                    pairing.groupInternal(pair_accumulator, group);
                }
            } else {
                pairing.group2all(pair_accumulator, group, change_data.relative_atom_indices);
                if (change_data.internal) {
                    pairing.groupInternal(pair_accumulator, group, change_data.relative_atom_indices);
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
    template <RequireEnergyAccumulator TAccumulator>
    void accumulateSpeciation(TAccumulator& pair_accumulator, const Change& change) {
        assert(change.matter_change);
        const auto &moved = change.touchedGroupIndex(); // index of moved groups
        const auto fixed =
            indexComplement(spc.groups.size(), moved) | ranges::to<std::vector>; // index of static groups
        auto filter_active = [](int size) { return ranges::views::filter([size](const auto i) { return i < size; }); };

        // loop over all changed groups
        for (auto change_group1_it = change.groups.begin(); change_group1_it < change.groups.end(); ++change_group1_it) {
            const auto& group1 = spc.groups.at(change_group1_it->group_index);
            // filter only active particles
            const auto index1 =
                change_group1_it->relative_atom_indices | filter_active(group1.size()) | ranges::to<std::vector>;
            if (!index1.empty()) {
                // particles added into the group: compute (changed group) <-> (static group)
                pairing.group2groups(pair_accumulator, group1, fixed, index1);
            }
            // loop over successor changed groups (hence avoid double counting group1×group2 and group2×group1)
            for (auto change_group2_it = std::next(change_group1_it); change_group2_it < change.groups.end(); ++change_group2_it) {
                const auto& group2 = spc.groups.at(change_group2_it->group_index);
                const auto index2 =
                    change_group2_it->relative_atom_indices | filter_active(group2.size()) | ranges::to<std::vector>;
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
    template <RequireEnergyAccumulator TAccumulator>
    void accumulate(TAccumulator& pair_accumulator, const Change& change) {
        assert(std::is_sorted(change.groups.begin(), change.groups.end()));
        if (change.everything) {
            pairing.all(pair_accumulator);
        } else if (change.volume_change) {
            // sum all interaction energies except the internal energies of incompressible molecules
            pairing.all(pair_accumulator, [](auto& group) { return group.isAtomic() || group.traits().compressible; });
        } else if (!change.matter_change) {
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
    void group2group(Accumulator& pair_accumulator, const Space::GroupType& group1, const Space::GroupType& group2) {
        pairing.group2group(std::forward<Accumulator&>(pair_accumulator), std::forward<const Space::GroupType&>(group1),
                            std::forward<const Space::GroupType&>(group2));
    }
};

class NonbondedBase : public Energybase {
  public:
    virtual double particleParticleEnergy(const Particle &particle1, const Particle &particle2) = 0;
    virtual double groupGroupEnergy(const Group& group1, const Group& group2) = 0;
};

/**
 * @brief Computes change in the non-bonded energy, assuming pairwise additive energy terms.
 *
 * @tparam TPairEnergy  a functor to compute non-bonded energy between two particles
 * @tparam TPairingPolicy  pairing policy to effectively sum up the pairwise additive non-bonded energy
 */
template <RequirePairEnergy TPairEnergy, typename TPairingPolicy> class Nonbonded : public NonbondedBase {
  protected:
    const Space& spc;        //!< space to operate on
    TPairEnergy pair_energy; //!< a functor to compute non-bonded energy between two particles, see PairEnergy
    TPairingPolicy pairing;  //!< pairing policy to effectively sum up the pairwise additive non-bonded energy
    std::shared_ptr<EnergyAccumulatorBase>
        energy_accumulator; //!< energy accumulator used for storing and summing pair-wise energies

  public:
    Nonbonded(const json& j, Space& spc, BasePointerVector<Energybase>& pot)
        : spc(spc), pair_energy(spc, pot), pairing(spc) {
        name = "nonbonded";
        from_json(j);
        energy_accumulator = createEnergyAccumulator(j, pair_energy, 0.0);
        energy_accumulator->reserve(spc.numParticles()); // attempt to reduce memory fragmentation
    }

    double particleParticleEnergy(const Particle& particle1, const Particle& particle2) override {
        return pair_energy(particle1, particle2);
    }

    double groupGroupEnergy(const Group &group1, const Group& group2) override {
        InstantEnergyAccumulator<TPairEnergy> accumulator(pair_energy);
        pairing.group2group(accumulator, group1, group2);
        return static_cast<double>(accumulator);
    }

    void from_json(const json &j) {
        pair_energy.from_json(j);
        pairing.from_json(j);
    }

    void to_json(json &j) const override {
        pair_energy.to_json(j);
        pairing.to_json(j);
        energy_accumulator->to_json(j);
    }

    double energy(const Change& change) override {
        energy_accumulator->clear();
        // down-cast to avoid slow, virtual function calls:
        if (auto ptr = std::dynamic_pointer_cast<InstantEnergyAccumulator<TPairEnergy>>(energy_accumulator)) {
            pairing.accumulate(*ptr, change);
        } else if (auto ptr = std::dynamic_pointer_cast<DelayedEnergyAccumulator<TPairEnergy>>(energy_accumulator)) {
            pairing.accumulate(*ptr, change);
        } else {
            pairing.accumulate(*energy_accumulator, change);
        }
        return static_cast<double>(*energy_accumulator);
    }

    /**
     * @brief Calculates the force on all particles.
     *
     * @todo A stub. Change to reflect only active particle, see Space::activeParticles().
     */
    void force(std::vector<Point> &forces) override {
        // just a temporary hack; perhaps better to allow PairForce instead of the PairEnergy template
        assert(forces.size() == spc.particles.size() && "the forces size must match the particle size");
        for (size_t i = 0; i < spc.particles.size() - 1; ++i) {
            for (size_t j = i + 1; j < spc.particles.size(); ++j) {
                const Point f = pair_energy.force(spc.particles[i], spc.particles[j]);
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
template <RequirePairEnergy TPairEnergy, typename TPairingPolicy>
class NonbondedCached : public Nonbonded<TPairEnergy, TPairingPolicy> {
    using Base = Nonbonded<TPairEnergy, TPairingPolicy>;
    using TAccumulator = InstantEnergyAccumulator<TPairEnergy>;
    Eigen::MatrixXf energy_cache;
    using Base::spc;

    template <typename TGroup>
    double g2g(const TGroup &g1, const TGroup &g2) {
        int i = &g1 - spc.groups.data();
        int j = &g2 - spc.groups.data();
        if (j < i) {
            std::swap(i, j);
        }
        if (Energybase::state == Energybase::MonteCarloState::TRIAL) { // if this is from the trial system
            TAccumulator energy_accumulator(Base::pair_energy);
            Base::pairing.group2group(energy_accumulator, g1, g2);
            energy_cache(i, j) = static_cast<double>(energy_accumulator);  // update the cache
        }
        return energy_cache(i, j); // return (cached) value
    }

    template <typename TGroup>
    double g2g(const TGroup& g1, const TGroup& g2, [[maybe_unused]] const std::vector<std::size_t>& index) {
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
                Base::pairing.group2group(u, spc.groups.at(i), spc.groups.at(j));
                energy_cache(i, j) = static_cast<double>(u);
            }
        }
    }

    double energy(const Change& change) override {
        // Only g2g may be called there to compute (and cache) energy!
        double energy_sum = 0.0;
        if (change) {
            if (change.everything || change.volume_change) {
                for (auto i = spc.groups.begin(); i < spc.groups.end(); ++i) {
                    for (auto j = std::next(i); j < Base::spc.groups.end(); ++j) {
                        energy_sum += g2g(*i, *j);
                    }
                }
            } else {
                if (change.groups.size() == 1) { // if exactly ONE molecule is changed
                    auto &d = change.groups[0];
                    auto& g1 = spc.groups.at(d.group_index);
                    for (auto g2_it = spc.groups.begin(); g2_it < spc.groups.end(); ++g2_it) {
                        if (&g1 != &(*g2_it)) {
                            energy_sum += g2g(g1, *g2_it, d.relative_atom_indices);
                        }
                    }
                } else {                                     // many molecules are changed
                    auto moved = change.touchedGroupIndex(); // index of moved groups
                    // moved<->moved
                    if (change.moved_to_moved_interactions) {
                        for (auto i = moved.begin(); i < moved.end(); ++i) {
                            for (auto j = std::next(i); j < moved.end(); ++j) {
                                energy_sum += g2g(spc.groups[*i], spc.groups[*j]);
                            }
                        }
                    }
                    // moved<->static
#if true
                    // classic version
                    const auto fixed = indexComplement(spc.groups.size(), moved) | ranges::to_vector; // static groups
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
    void sync(Energybase* base_ptr, const Change& change) override {
        auto other = dynamic_cast<decltype(this)>(base_ptr);
        assert(other);
        if (change.everything || change.volume_change) {
            energy_cache.triangularView<Eigen::StrictlyUpper>() =
                (other->energy_cache).template triangularView<Eigen::StrictlyUpper>();
        } else {
            for (const auto& d : change.groups) {
                for (int i = 0; i < d.group_index; i++) {
                    energy_cache(i, d.group_index) = other->energy_cache(i, d.group_index);
                }
                for (size_t i = d.group_index + 1; i < spc.groups.size(); i++) {
                    energy_cache(d.group_index, i) = other->energy_cache(d.group_index, i);
                }
            }
        }
    }
};

#ifdef ENABLE_FREESASA
/**
 * @brief Interface to the FreeSASA C-library. Experimental and unoptimized.
 * https://freesasa.github.io/
 *
 * @todo - Implement partial evaluation refelcting `change` object
 *       - Average volume currently mixes accepted/rejected states
 */
class FreeSASAEnergy : public Energybase {
  private:
    std::vector<double> positions; //!< Flattened position buffer for all particles
    std::vector<double> radii;     //!< Radii buffer for all particles
    std::vector<double> sasa;      //!< Target buffer for calculated surface areas

    const Space& spc;
    double cosolute_molarity = 0.;                       //!< co-solute concentration (mol/l)
    std::unique_ptr<freesasa_parameters_fwd> parameters; //!< Parameters for freesasa
    Average<double> mean_surface_area;

    void to_json(json &j) const override;
    void sync(Energybase* energybase_ptr, const Change& change) override;
    void updateSASA(const Change& change);
    void init() override;

    /**
     * @brief Copies radii from Space to internal buffer
     * @param begin Iterator to first particle
     * @param end Iterator to beyond last particle
     * @param change Change object (currently unused)
     */
    template <typename Tfirst, typename Tend>
    void updateRadii(Tfirst begin, Tend end, [[maybe_unused]] const Change& change) {
        const auto number_of_particles = std::distance(begin, end);
        radii.clear();
        radii.reserve(number_of_particles);
        std::transform(begin, end, std::back_inserter(radii),
                       [](const Particle& particle) { return particle.traits().sigma * 0.5; });
    }

    /**
     * @brief Copies positions from Space to internal (flattened) buffer
     * @param begin Iterator to first particle
     * @param end Iterator to beyond last particle
     * @param change Change object (currently unused)
     */
    template <typename Tfirst, typename Tend>
    void updatePositions(Tfirst begin, Tend end, [[maybe_unused]] const Change& change) {
        const auto number_of_particles = std::distance(begin, end);
        positions.clear();
        positions.reserve(3 * number_of_particles);
        for (const auto& particle : spc.activeParticles()) {
            const auto* xyz = particle.pos.data();
            positions.insert(positions.end(), xyz, xyz + 3);
        }
    }

  public:
    /**
     * @param spc
     * @param cosolute_molarity in particles per angstrom cubed
     * @param probe_radius in angstrom
     */
    FreeSASAEnergy(const Space& spc, double cosolute_molarity, double probe_radius);
    FreeSASAEnergy(const json& j, const Space& spc);
    double energy(const Change& change) override;
    const std::vector<double>& getAreas() const { return sasa; }
}; //!< SASA energy from transfer free energies
#endif

/**
 * @brief class for calculating SASA energies calculating SASA of each particle every step
 *
 * This calculates SASA for all particles, i.e. the `Change` object is ignored and a full
 * (and slow) system update is performed on each call to `energy()`. This class is used as
 * a reference and base class for more clever implementations.
 */
class SASAEnergyReference : public Energybase {
  private:
    void to_json(json& j) const override;
    void sync(Energybase* energybase_ptr, const Change& change) override;

  protected:
    void init() override;
    using index_type = size_t;
    const Space& spc;                     //!< Space to operate on
    std::vector<double> areas;            //!< Target buffer for calculated surface areas
    double cosolute_molarity = 0.0;       //!< co-solute concentration (mol/l)
    std::unique_ptr<SASA::SASABase> sasa; //!< performs neighbour searching and subsequent sasa calculation

    /** absolute index of particle in space particle vector */
    inline auto indexOf(const Particle& particle) const {
        return static_cast<index_type>(std::addressof(particle) - std::addressof(spc.particles.at(0)));
    }

  public:
    SASAEnergyReference(const Space& spc, double cosolute_molarity, double probe_radius, int slices_per_atom = 25,
                        bool dense_container = true);
    SASAEnergyReference(const json& j, const Space& spc);
    const std::vector<double>& getAreas() const;
    double energy(const Change& change) override;
};

/**
 * @brief class for calculating SASA energies calculating SASA of particles based on change object every step
 */
class SASAEnergy : public SASAEnergyReference {
  private:
    std::vector<std::vector<index_type>>
        current_neighbours; //!< holds cached neighbour indices for each particle in ParticleVector
    std::vector<index_type> changed_indices; //!< paritcle indices whose SASA changed based on change object

    void sync(Energybase* energybase_ptr, const Change& change) override;
    void init() override;

    void updateChangedIndices(const Change& change);
    void insertChangedNeighboursOf(index_type index, std::set<index_type>& target_indices) const;

  public:
    SASAEnergy(const Space& spc, double cosolute_molarity, double probe_radius, int slices_per_atom = 25,
               bool dense_container = true);
    SASAEnergy(const json& j, const Space& spc);
    double energy(const Change& change) override;
};

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
    double energy(const Change& change) override;
};

/**
 * @brief Aggregate and sum energy terms
 */
class Hamiltonian : public Energybase, public BasePointerVector<Energybase> {
  private:
    double maximum_allowed_energy = pc::infty; //!< Maximum allowed energy change
    std::vector<double> latest_energies;       //!< Placeholder for the lastest energies for each energy term
    decltype(vec)& energy_terms;               //!< Alias for `vec`
    void addEwald(const json& j, Space& spc);  //!< Adds an instance of reciprocal space Ewald energies (if appropriate)
    void checkBondedMolecules() const;         //!< Warn if bonded molecules and no bonded energy term
    void to_json(json& j) const override;
    void force(PointVector& forces) override;
    std::unique_ptr<Energybase> createEnergy(Space& spc, const std::string& name, const json& j);

  public:
    Hamiltonian(Space& spc, const json& j);
    void init() override;
    void updateState(const Change& change) override;
    void sync(Energybase* other_hamiltonian, const Change& change) override;
    double energy(const Change& change) override;      //!< Energy due to changes
    const std::vector<double>& latestEnergies() const; //!< Energies for each term from the latest call to `energy()`
};
} // namespace Energy
} // namespace Faunus
