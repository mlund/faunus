#pragma once

#include "core.h"
//#include "auxiliary.h"
#include "particle.h"

namespace Faunus {

namespace Geometry {
// this typedef is re-defined here so as not to have geometry.h as a dependency.
// consider another fix for this.
typedef std::function<Point(const Point&, const Point&)> DistanceFunction;
} // namespace Geometry

namespace Potential {
/**
 * @brief Base class for bonded potentials
 *
 * This stores data on the bond type; atom indices; json keywords;
 * and potentially also the `energyFunc` functor (nullptr per default).
 *
 * The `forceFunc` functor returns a vector of forces acting on the
 * participating atoms
 */
struct BondData {
    enum Variant { HARMONIC = 0, FENE, FENEWCA, HARMONIC_TORSION, GROMOS_TORSION, PERIODIC_DIHEDRAL, INVALID };
    std::vector<int> index; //!< Index of participating atoms

    /** Calculates potential energy of bonded atoms(kT) */
    std::function<double(Geometry::DistanceFunction)> energyFunc = nullptr;

    using ParticleForce = std::pair<int, Point>; //!< Force (second) on particle w. absolute index (first)

    /** Calculates forces on bonded atoms (kT/Å) */
    std::function<std::vector<ParticleForce>(Geometry::DistanceFunction)> forceFunc = nullptr;

    virtual void from_json(const json&) = 0;
    virtual void to_json(json&) const = 0;
    virtual int numindex() const = 0;                                    //!< Required number of atom indices for bond
    virtual Variant type() const = 0;                                    //!< Returns bond type (sett `Variant` enum)
    virtual std::shared_ptr<BondData> clone() const = 0;                 //!< Make shared pointer *copy* of data
    virtual void setEnergyFunction(const ParticleVector& particles) = 0; //!< Set energy function; store particles ref.
    bool hasEnergyFunction() const;                                      //!< test if energy function has been set
    bool hasForceFunction() const;                                       //!< test if force function has been set
    void shift(const int offset);                                        //!< Add offset to index
    BondData() = default;
    explicit BondData(const std::vector<int>& index);
    virtual ~BondData() = default;
};

NLOHMANN_JSON_SERIALIZE_ENUM(BondData::Variant, {{BondData::Variant::INVALID, nullptr},
                                                 {BondData::Variant::HARMONIC, "harmonic"},
                                                 {BondData::Variant::FENE, "fene"},
                                                 {BondData::Variant::FENEWCA, "fene+wca"},
                                                 {BondData::Variant::HARMONIC_TORSION, "harmonic_torsion"},
                                                 {BondData::Variant::GROMOS_TORSION, "gromos_torsion"},
                                                 {BondData::Variant::PERIODIC_DIHEDRAL, "periodic_dihedral"}})

struct StretchData : public BondData {
    int numindex() const override;
    StretchData() = default;
    StretchData(const std::vector<int>& index);
};

struct TorsionData : public BondData {
    int numindex() const override;
    TorsionData() = default;
    TorsionData(const std::vector<int>& index);
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
    HarmonicBond(double k, double req, const std::vector<int>& index);
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
    FENEBond(double k, double rmax, const std::vector<int>& index);
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
    FENEWCABond(double k, double rmax, double epsilon, double sigma, const std::vector<int>& index);
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
    HarmonicTorsion(double k, double aeq, const std::vector<int>& index);
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
    GromosTorsion(double k, double cos_aeq, const std::vector<int>& index);
};

/**
 * @brief Periodic dihedral
 *
 * U(a) = k * (1 + cos(n * a - phi))
 */
struct PeriodicDihedral : public BondData {
    double force_constant = 0.0;
    double dihedral_angle = 0.0;
    double periodicity = 1.0;
    int numindex() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json& j) override;
    void to_json(json& j) const override;
    Variant type() const override;
    void setEnergyFunction(const ParticleVector& particles) override;
    PeriodicDihedral() = default;
    PeriodicDihedral(double k, double phi, double n, const std::vector<int>& index);
};

/*
 * Serialize to/from json
 */

void from_json(const json& j, std::shared_ptr<BondData>& bond);
void to_json(json& j, const std::shared_ptr<const BondData>& bond);
void to_json(json& j, const BondData& bond);

[[deprecated("Use bonds.find<TClass>() method instead.")]] inline auto
filterBonds(const std::vector<std::shared_ptr<BondData>>& bonds, BondData::Variant bondtype) {
    std::vector<std::shared_ptr<BondData>> filt;
    filt.reserve(bonds.size());
    std::copy_if(bonds.begin(), bonds.end(), std::back_inserter(filt),
                 [bondtype](const auto& d) { return d->type() == bondtype; });
    return filt;
} //!< Filter bond container for matching bond type and return _reference_ to original

} // namespace Potential
} // namespace Faunus
