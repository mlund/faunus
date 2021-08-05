#pragma once
#include "core.h"
#include "particle.h"

namespace Faunus::Geometry {
// @todo typedef re-defined to avoid geometry.h dependency; consider another fix for this.
typedef std::function<Point(const Point&, const Point&)> DistanceFunction;
} // namespace Faunus::Geometry

namespace Faunus::Potential {

/**
 * @brief Base class for bonded potentials
 *
 * This stores data on the bond type; atom indices; json keywords;
 * and potentially also the `energyFunc` functor (nullptr per default).
 *
 * The `forceFunc` functor returns a vector of forces acting on the
 * participating atoms.
 *
 * @todo Memory inefficient to store energy and force functors. Virtual functions?
 */
struct BondData {
    enum Variant { HARMONIC = 0, FENE, FENEWCA, HARMONIC_TORSION, GROMOS_TORSION, PERIODIC_DIHEDRAL, HARMONIC_DIHEDRAL, INVALID };
    std::vector<int> indices; //!< Absolute indiced of participating particles

    /** Calculates potential energy of bonded atoms(kT) */
    std::function<double(Geometry::DistanceFunction)> energyFunc = nullptr;

    using IndexAndForce = std::pair<int, Point>; //!< Force (second) on particle w. absolute index (first)

    /** Calculates forces on bonded atoms (kT/Å) */
    std::function<std::vector<IndexAndForce>(Geometry::DistanceFunction)> forceFunc = nullptr;

    virtual void from_json(const json&) = 0;
    virtual void to_json(json&) const = 0;
    virtual int numindex() const = 0;                                    //!< Required number of atom indices for bond
    virtual Variant type() const = 0;                                    //!< Returns bond type (sett `Variant` enum)
    virtual std::shared_ptr<BondData> clone() const = 0;                 //!< Make shared pointer *copy* of data
    virtual void setEnergyFunction(const ParticleVector& particles) = 0; //!< Set energy function; store particles ref.
    bool hasEnergyFunction() const;                                      //!< test if energy function has been set
    bool hasForceFunction() const;                                       //!< test if force function has been set
    void shiftIndices(const int offset);                                 //!< Add offset to particle indices
    BondData() = default;
    explicit BondData(const std::vector<int>& indices);
    virtual ~BondData() = default;
};

NLOHMANN_JSON_SERIALIZE_ENUM(BondData::Variant, {{BondData::Variant::INVALID, nullptr},
                                                 {BondData::Variant::HARMONIC, "harmonic"},
                                                 {BondData::Variant::FENE, "fene"},
                                                 {BondData::Variant::FENEWCA, "fene+wca"},
                                                 {BondData::Variant::HARMONIC_TORSION, "harmonic_torsion"},
                                                 {BondData::Variant::GROMOS_TORSION, "gromos_torsion"},
                                                 {BondData::Variant::PERIODIC_DIHEDRAL, "periodic_dihedral"},
                                                 {BondData::Variant::HARMONIC_DIHEDRAL, "harmonic_dihedral"}})

struct StretchData : public BondData {
    int numindex() const override;
    StretchData() = default;
    StretchData(const std::vector<int>& indices);
};

struct TorsionData : public BondData {
    int numindex() const override;
    TorsionData() = default;
    TorsionData(const std::vector<int>& indices);
};

/**
 * @brief Harmonic Bond
 *
 * U(r) = k/2 * (r - r_eq)^2
 */
struct HarmonicBond : public StretchData {
    double half_force_constant = 0.0;
    double equilibrium_distance = 0.0;
    Variant type() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json& j) override;
    void to_json(json& j) const override;
    void setEnergyFunction(const ParticleVector& particles) override; //!< Set energy and force functors
    HarmonicBond() = default;
    HarmonicBond(double k, double req, const std::vector<int>& indices);
};

/**
 * @brief FENE bond
 *
 * U(r) = -k/2 * r_max^2 * ln(1 - r^2 / r_max^2) if r < r_max, ∞ otherwise
 */
struct FENEBond : public StretchData {
    double half_force_constant = 0.0;
    double max_squared_distance = 0.0;
    Variant type() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json& j) override;
    void to_json(json& j) const override;
    void setEnergyFunction(const ParticleVector& particles) override;
    FENEBond() = default;
    FENEBond(double k, double rmax, const std::vector<int>& indices);
};

/**
 * @brief FENE+WCA bond
 */
struct FENEWCABond : public StretchData {
    double half_force_constant = 0.0;
    double max_distance_squared = 0.0;
    double epsilon = 0.0;
    double sigma_squared = 0.0;
    std::array<double, 4> k = {{0.0, 0.0, 0.0, 0.0}};
    Variant type() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json& j) override;
    void to_json(json& j) const override;
    void setEnergyFunction(const ParticleVector& calculateDistance) override;
    FENEWCABond() = default;
    FENEWCABond(double k, double rmax, double epsilon, double sigma, const std::vector<int>& indices);
};

/**
 * @brief Harmonic torsion
 *
 * U(a) = k/2 * (a - a_eq)^2
 */
struct HarmonicTorsion : public TorsionData {
    double half_force_constant = 0.0;
    double equilibrium_angle = 0.0;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json& j) override;
    void to_json(json& j) const override;
    Variant type() const override;
    void setEnergyFunction(const ParticleVector& particles) override;
    HarmonicTorsion() = default;
    HarmonicTorsion(double k, double aeq, const std::vector<int>& indices);
};

/**
 * @brief Gromos torsion
 *
 * U(a) = k/2 * (cos(a) - cos(a_eq))^2
 */
struct GromosTorsion : public TorsionData {
    double half_force_constant = 0.0;
    double cosine_equilibrium_angle = 0.0;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json& j) override;
    void to_json(json& j) const override;
    Variant type() const override;
    void setEnergyFunction(const ParticleVector& calculateDistance) override;
    GromosTorsion() = default;
    GromosTorsion(double k, double cos_aeq, const std::vector<int>& indices);
};

/**
 * @brief Periodic dihedral
 *
 * U(a) = k * (1 + cos(n * a - phi))
 */
struct PeriodicDihedral : public BondData {
    double force_constant = 0.0;
    double phase_angle = 0.0;
    double periodicity = 1.0;
    int numindex() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json& j) override;
    void to_json(json& j) const override;
    Variant type() const override;
    void setEnergyFunction(const ParticleVector& particles) override;
    PeriodicDihedral() = default;
    PeriodicDihedral(double k, double phi, double n, const std::vector<int>& indices);
};

/**
 * @brief Harmonic dihedral
 *
 * U(a) = k/2 * (phi - phi_eq)^2)
 */
struct HarmonicDihedral : public BondData {
    double half_force_constant = 0.0;
    double equilibrium_dihedral = 0.0;
    int numindex() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json& j) override;
    void to_json(json& j) const override;
    Variant type() const override;
    void setEnergyFunction(const ParticleVector& particles) override;
    HarmonicDihedral() = default;
    HarmonicDihedral(double k, double deq, const std::vector<int>& indices);
};

void from_json(const json& j, std::shared_ptr<BondData>& bond);
void to_json(json& j, const std::shared_ptr<const BondData>& bond);
void to_json(json& j, const BondData& bond);

} // namespace Faunus::Potential