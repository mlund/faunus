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
    virtual ~BondData();
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

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] BondData") {
    std::shared_ptr<BondData> b;

    // exact match required
    CHECK_THROWS(b = R"({ "harmoNIC": {"index":[2,3], "k":0.5, "req":2.1}} )"_json;);

    // test harmonic
    SUBCASE("HarmonicBond") {
        json j = R"({ "harmonic": {"index":[2,3], "k":0.5, "req":2.1}} )"_json;
        b = j;
        CHECK(j == json(b));
        CHECK_THROWS(b = R"({"harmonic": { "index":[2], "k":0.5, "req":2.1}} )"_json);
        CHECK_THROWS(b = R"({"harmonic": { "index":[2,3], "req":2.1}} )"_json);
        CHECK_THROWS(b = R"({"harmonic": { "index":[2,3], "k":2.1}} )"_json);
    }

    // test fene
    SUBCASE("FENEBond") {
        json j = R"({"fene": { "index":[2,3], "k":1, "rmax":2.1 }} )"_json;
        b = j;
        CHECK(j == json(b));
        CHECK_THROWS(b = R"({"fene": { "index":[2,3,4], "k":1, "rmax":2.1}} )"_json);
        CHECK_THROWS(b = R"({"fene": { "index":[2,3], "rmax":2.1}} )"_json);
        CHECK_THROWS(b = R"({"fene": { "index":[2,3], "k":1}} )"_json);
    }

    // test fene+wca
    SUBCASE("FENEWCABond") {
        json j = R"({"fene+wca": { "index":[2,3], "k":1, "rmax":2.1, "eps":2.48, "sigma":2}} )"_json;
        b = j;
        CHECK(j == json(b));
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3,4], "k":1, "rmax":2.1, "eps":2.48, "sigma":2}} )"_json);
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3], "rmax":2.1, "eps":2.48, "sigma":2}} )"_json);
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3], "k":1, "eps":2.48, "sigma":2}} )"_json);
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3], "k":1, "rmax":2.1, "eps":2.48}} )"_json);
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3], "k":1, "rmax":2.1, "sigma":2}} )"_json);
    }

    // test harmonic
    SUBCASE("HarmonicTorsion") {
        json j = R"({ "harmonic_torsion": {"index":[0,1,2], "k":0.5, "aeq":60}} )"_json;
        b = j;
        CHECK(j == json(b));
        CHECK_THROWS(b = R"({"harmonic_torsion": { "index":[2], "k":0.5, "aeq":2.1}} )"_json);
        CHECK_THROWS(b = R"({"harmonic_torsion": { "index":[0,1,2], "aeq":2.1}} )"_json);
        CHECK_THROWS(b = R"({"harmonic_torsion": { "index":[0,1,3], "k":2.1}} )"_json);
    }

    // test bond filter
    SUBCASE("filterBonds()") {
        std::vector<std::shared_ptr<BondData>> bonds = {
            R"({"fene":      {"index":[2,3], "k":1, "rmax":2.1, "eps":2.48}} )"_json,
            R"({"harmonic" : {"index":[2,3], "k":0.5, "req":2.1} } )"_json};
        auto filt = filterBonds(bonds, BondData::HARMONIC);
        CHECK(filt.size() == 1);
        CHECK(filt[0]->type() == BondData::HARMONIC);
        CHECK(filt[0] == bonds[1]); // filt should contain references to bonds
    }
}
#endif
} // namespace Potential
} // namespace Faunus
