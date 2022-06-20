#pragma once
#include "energy.h"

namespace Faunus {

namespace Potential {
class PairPotentialBase;
}

namespace Energy {
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
 * @brief Determines if two groups are separated beyond the cutoff distance.
 *
 * The distance between centers of mass is considered. The cutoff distance can be specified independently for each
 * group pair to override the default value.
 *
 * @see GroupPairingPolicy
 */
class GroupCutoff {
    double default_cutoff_squared = pc::max_value;
    PairMatrix<double> cutoff_squared; //!< matrix with group-to-group cutoff distances squared in angstrom squared
    Space::GeometryType& geometry;     //!< geometry to compute the inter group distance with
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
    template <typename... Args> inline auto operator()(Args&&... args) { return cut(std::forward<Args>(args)...); }

    /**
     * @brief Sets the geometry.
     * @param geometry  geometry to compute the inter group distance with
     */
    explicit GroupCutoff(Space::GeometryType& geometry);
};

void from_json(const json&, GroupCutoff&);
void to_json(json&, const GroupCutoff&);

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

    const Space& getSpace() const { return spc; }

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

    explicit EnergyAccumulatorBase(double value);
    virtual ~EnergyAccumulatorBase() = default;
    virtual void reserve(size_t number_of_particles);
    virtual void clear();
    virtual void from_json(const json& j);
    virtual void to_json(json& j) const;

    virtual explicit operator double();
    virtual EnergyAccumulatorBase& operator=(double new_value) = 0;
    virtual EnergyAccumulatorBase& operator+=(double new_value) = 0;
    virtual EnergyAccumulatorBase& operator+=(ParticlePair&& pair) = 0;
    virtual void updateState(const Change& change); //!< If internal state needs updating with a system change
    virtual void sync(const EnergyAccumulatorBase& other, const Change& change); //!< Sync state with other accumulator

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
        : EnergyAccumulatorBase(value)
        , pair_energy(pair_energy) {}

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

    void from_json(const json& j) override {
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
  protected:
    std::vector<ParticlePair> particle_pairs;
    const PairEnergy& pair_energy; //!< recipe to compute non-bonded energy between two particles, see PairEnergy
    const size_t max_particles_in_buffer = 10000; //!< this can be modified to suit memory requirements

  public:
    explicit DelayedEnergyAccumulator(const PairEnergy& pair_energy, const double value = 0.0)
        : EnergyAccumulatorBase(value)
        , pair_energy(pair_energy) {}

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
                                                               double initial_energy) {
    std::unique_ptr<EnergyAccumulatorBase> accumulator;
    if (j.value("summation_policy", EnergyAccumulatorBase::Scheme::SERIAL) != EnergyAccumulatorBase::Scheme::SERIAL) {
        accumulator = std::make_unique<DelayedEnergyAccumulator<TPairEnergy>>(pair_energy, initial_energy);
        faunus_logger->debug("activated delayed energy summation");
    } else {
        accumulator = std::make_unique<InstantEnergyAccumulator<TPairEnergy>>(pair_energy, initial_energy);
        faunus_logger->debug("activated instant energy summation");
    }
    accumulator->from_json(j);
    return accumulator;
}

} // namespace Energy
} // namespace Faunus