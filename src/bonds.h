#pragma once

#include "core.h"
#include "auxiliary.h"
#include "particle.h"

namespace Faunus {

namespace Geometry {
// this typedef is re-defined here so as not to have geometry.h as a dependency.
// consider another fix for this.
typedef std::function<Point(const Point &, const Point &)> DistanceFunction;
} // namespace Geometry

namespace Potential {
/**
 * @brief Base class for bonded potentials
 *
 * This stores data on the bond type; atom indices; json keywords;
 * and potentially also the energy function (nullptr per default).
 */
struct BondData {
    enum Variant { HARMONIC = 0, FENE, FENEWCA, HARMONIC_TORSION, GROMOS_TORSION, PERIODIC_DIHEDRAL, NONE };
    std::vector<int> index;
    std::function<double(Geometry::DistanceFunction)> energy = nullptr; //!< potential energy (kT)

    virtual void from_json(const json &) = 0;
    virtual void to_json(json &) const = 0;
    virtual int numindex() const = 0;                    //!< Required number of atom indices for bond
    virtual Variant type() const = 0;                    //!< Returns bond type (sett `Variant` enum)
    virtual std::string name() const = 0;                //!< Name/key of bond type used in for json I/O
    virtual std::shared_ptr<BondData> clone() const = 0; //!< Make shared pointer *copy* of data
    bool hasEnergyFunction() const;                      //!< test if energy function has been set
    void shift(int offset);                              //!< Shift indices
    BondData() = default;
    BondData(const std::vector<int> &index);
    virtual ~BondData() = default;
};

struct StretchData : public BondData {
    int numindex() const override { return 2; }
    StretchData() = default;
    StretchData(const std::vector<int> &index) : BondData(index) {};
};

struct TorsionData : public BondData {
    int numindex() const override { return 3; }
    TorsionData() = default;
    TorsionData(const std::vector<int> &index) : BondData(index) {};
};

/**
 * @brief Harmonic Bond
 *
 * U(r) = k/2 * (r - r_eq)^2
 */
struct HarmonicBond : public StretchData {
    double k_half = 0, req = 0;
    Variant type() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
    HarmonicBond() = default;
    HarmonicBond(double k, double req, const std::vector<int> &index);
};

/**
 * @brief FENE bond
 *
 * U(r) = -k/2 * r_max^2 * ln(1 - r^2 / r_max^2) if r < r_max, ∞ otherwise
 */
struct FENEBond : public StretchData {
    double k_half = 0, rmax_squared = 0;
    Variant type() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
    FENEBond() = default;
    FENEBond(double k, double rmax, const std::vector<int> &index);
};

/**
 * @brief FENE+WCA bond
 */
struct FENEWCABond : public StretchData {
    double k_half = 0, rmax_squared = 0, epsilon = 0, sigma_squared = 0;
    std::array<double, 4> k = {{0, 0, 0, 0}};
    Variant type() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
    FENEWCABond() = default;
    FENEWCABond(double k, double rmax, double epsilon, double sigma, const std::vector<int> &index);
};

/**
 * @brief Harmonic torsion
 *
 * U(a) = k/2 * (a - a_eq)^2
 */
struct HarmonicTorsion : public TorsionData {
    double k_half = 0, aeq = 0;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Variant type() const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
    HarmonicTorsion() = default;
    HarmonicTorsion(double k, double aeq, const std::vector<int> &index);
};

/**
 * @brief Gromos torsion
 *
 * U(a) = k/2 * (cos(a) - cos(a_eq))^2
 */
struct GromosTorsion : public TorsionData {
    double k_half = 0, cos_aeq = 0;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Variant type() const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
    GromosTorsion() = default;
    GromosTorsion(double k, double cos_aeq, const std::vector<int> &index);
};

/**
 * @brief Periodic dihedral
 *
 * U(a) = k * (1 + cos(n * a - phi))
 */
struct PeriodicDihedral : public BondData {
    double k = 0, phi = 0, n = 1;
    int numindex() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Variant type() const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
    PeriodicDihedral() = default;
    PeriodicDihedral(double k, double phi, double n, const std::vector<int> &index);
};

/*
 * Serialize to/from json
 */

void from_json(const json &j, std::shared_ptr<BondData> &b);
void to_json(json &j, const std::shared_ptr<const BondData> &b);
void to_json(Faunus::json &j, const BondData &b);

void setBondEnergyFunction(std::shared_ptr<BondData> &b,
                           const ParticleVector &p); //!< Set the bond energy function of `BondData` which
                                                     //!< require a reference to the particle vector

[[deprecated("Use bonds.find<TClass>() method instead.")]]
inline auto filterBonds(const std::vector<std::shared_ptr<BondData>> &bonds, BondData::Variant bondtype) {
    std::vector<std::shared_ptr<BondData>> filt;
    filt.reserve(bonds.size());
    std::copy_if(bonds.begin(), bonds.end(), std::back_inserter(filt),
                 [bondtype](const auto &d) { return d->type() == bondtype; });
    return filt;
} //!< Filter bond container for matching bond type and return _reference_ to original

} // namespace Potential
} // namespace Faunus
