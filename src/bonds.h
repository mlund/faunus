#pragma once

#include "core.h"
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
    bool exclude = false;           //!< True if exclusion of non-bonded interaction should be attempted
    bool keepelectrostatics = true; //!< If `exclude==true`, try to keep electrostatic interactions
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

/**
 * @brief Harmonic Bond
 */
struct HarmonicBond : public BondData {
    double k_half = 0, req = 0;
    int numindex() const override;
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
 */
struct FENEBond : public BondData {
    std::array<double, 4> k = {{0, 0, 0, 0}};
    int numindex() const override;
    Variant type() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
}; // end of FENE

/**
 * @brief FENE+WCA bond
 */
struct FENEWCABond : public BondData {
    std::array<double, 4> k = {{0, 0, 0, 0}};
    int numindex() const override;
    Variant type() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
}; // end of FENE+WCA

struct HarmonicTorsion : public BondData {
    double k_half = 0, aeq = 0;
    int numindex() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Variant type() const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
    HarmonicTorsion() = default;
    HarmonicTorsion(double k, double aeq, const std::vector<int> &index);
}; // end of HarmonicTorsion

struct GromosTorsion : public BondData {
    double k_half = 0, cos_aeq = 0;
    int numindex() const override;
    std::shared_ptr<BondData> clone() const override;
    void from_json(const json &j) override;
    void to_json(json &j) const override;
    Variant type() const override;
    std::string name() const override;
    void setEnergyFunction(const ParticleVector &p);
    GromosTorsion() = default;
    GromosTorsion(double k, double cos_aeq, const std::vector<int> &index);
}; // end of GromosTorsion

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
}; // end of PeriodicDihedral

/*
 * Serialize to/from json
 */

void to_json(json &j, const std::shared_ptr<BondData> &b);
void from_json(const json &j, std::shared_ptr<BondData> &b);

void setBondEnergyFunction(std::shared_ptr<BondData> &b,
                           const ParticleVector &p); //!< Set the bond energy function of `BondData` which
                                                     //!< require a reference to the particle vector

inline auto filterBonds(const std::vector<std::shared_ptr<BondData>> &bonds, BondData::Variant bondtype) {
    std::vector<std::shared_ptr<BondData>> filt;
    filt.reserve(bonds.size());
    std::copy_if(bonds.begin(), bonds.end(), std::back_inserter(filt),
                 [bondtype](const auto &d) { return d->type() == bondtype; });
    return filt;
} //!< Filter bond container for matching bond type and return _reference_ to original

} // namespace Potential
} // namespace Faunus
