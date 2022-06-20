#pragma once

#include "bonds.h"
#include "externalpotential.h" // Energybase implemented here
#include "potentials_base.h"
#include "sasa.h"
#include "space.h"
#include "celllistimpl.h"
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
    enum Policies { PBC, PBCEigen, IPBC, IPBCEigen, INVALID }; //!< Possible k-space updating schemes
    Policies policy = PBC;                                     //!< Policy for updating k-space
    explicit EwaldData(const json& j);                         //!< Initialize from json
};

NLOHMANN_JSON_SERIALIZE_ENUM(EwaldData::Policies, {
                                                      {EwaldData::INVALID, nullptr},
                                                      {EwaldData::PBC, "PBC"},
                                                      {EwaldData::PBCEigen, "PBCEigen"},
                                                      {EwaldData::IPBC, "IPBC"},
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
  private:
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
