#include <doctest/doctest.h>
#include "bonds.h"
#include "geometry.h"
#include "auxiliary.h"

namespace Faunus::Potential {

void to_json(Faunus::json& j, const std::shared_ptr<const BondData>& bond) {
    to_json(j, *bond);
}

void to_json(Faunus::json& j, const BondData& bond) {
    json val;
    bond.to_json(val);
    val["index"] = bond.index;
    j = {{bond.name(), val}};
}

void from_json(const json& j, std::shared_ptr<BondData>& bond) {
    try {
        const auto& [key, j_params] = jsonSingleItem(j);
        if (key == HarmonicBond().name()) {
            bond = std::make_shared<HarmonicBond>();
        } else if (key == FENEBond().name()) {
            bond = std::make_shared<FENEBond>();
        } else if (key == FENEWCABond().name()) {
            bond = std::make_shared<FENEWCABond>();
        } else if (key == HarmonicTorsion().name()) {
            bond = std::make_shared<HarmonicTorsion>();
        } else if (key == GromosTorsion().name()) {
            bond = std::make_shared<GromosTorsion>();
        } else if (key == PeriodicDihedral().name()) {
            bond = std::make_shared<PeriodicDihedral>();
        } else {
            throw ConfigurationError("'{}': unknown bond", key);
        }

        try {
            bond->from_json(j_params);
            bond->index = j_params.at("index").get<decltype(bond->index)>();
            if (bond->index.size() != bond->numindex()) {
                usageTip.pick(key);
                throw ConfigurationError("exactly {} indices required", bond->numindex());
            }
        } catch (std::exception& e) {
            throw ConfigurationError("'{}': {}", key, e.what());
        }
    } catch (std::exception& e) {
        throw ConfigurationError("bond: {}", e.what()).attachJson(j);
    }
}

void setBondEnergyFunction(std::shared_ptr<BondData> &b, const ParticleVector &p) {
    if (b->type() == BondData::HARMONIC)
        std::dynamic_pointer_cast<HarmonicBond>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::FENE)
        std::dynamic_pointer_cast<FENEBond>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::FENEWCA)
        std::dynamic_pointer_cast<FENEWCABond>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::HARMONIC_TORSION)
        std::dynamic_pointer_cast<HarmonicTorsion>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::GROMOS_TORSION)
        std::dynamic_pointer_cast<GromosTorsion>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::PERIODIC_DIHEDRAL)
        std::dynamic_pointer_cast<PeriodicDihedral>(b)->setEnergyFunction(p);
    else {
        assert(false); // we should never reach here
    }
}

void BondData::shift(int offset) {
    for (auto &i : index)
        i += offset;
}

bool BondData::hasEnergyFunction() const { return energyFunc != nullptr; }

BondData::BondData(const std::vector<int> &index) : index(index) {}

HarmonicBond::HarmonicBond(double k, double req, const std::vector<int> &index)
    : StretchData(index), k_half(k / 2), req(req) {}

void HarmonicBond::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2; // k
    req = j.at("req").get<double>() * 1.0_angstrom;                               // req
}

void HarmonicBond::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / 1.0_kJmol * 1.0_angstrom * 1.0_angstrom}, {"req", req / 1.0_angstrom}};
}

std::string HarmonicBond::name() const { return "harmonic"; }

std::shared_ptr<BondData> HarmonicBond::clone() const { return std::make_shared<HarmonicBond>(*this); }

BondData::Variant HarmonicBond::type() const { return BondData::HARMONIC; }

/**
 * @param particle Particle vector to all particles in the system
 *
 * This sets both the `energyFunc` and `forceFunc` function objects
 * for calculating the potential energy and the forces on the
 * participating atoms
 */
void HarmonicBond::setEnergyFunction(const ParticleVector &particle) {
    energyFunc = [&](Geometry::DistanceFunction dist) { // potential energy functor
        double d = req - dist(particle[index[0]].pos, particle[index[1]].pos).norm();
        return k_half * d * d; // kT/Å
    };
    forceFunc = [&](Geometry::DistanceFunction dist) -> std::vector<Point> { // force functor
        auto rvec = dist(particle[index[0]].pos, particle[index[1]].pos); // vector between points
        double r = rvec.norm();                                           // distance between particles
        auto force = 2.0 * k_half * (req - r) * rvec / r;                 // force on first particle
        return {force, -force};                                           // force on both particles (kT/Å)
    };
}

FENEBond::FENEBond(double k, double rmax, const std::vector<int> &index)
    : StretchData(index), k_half(k / 2), rmax_squared(rmax * rmax) {}

std::shared_ptr<BondData> FENEBond::clone() const { return std::make_shared<FENEBond>(*this); }

BondData::Variant FENEBond::type() const { return BondData::FENE; }

void FENEBond::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2; // k
    rmax_squared = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2); // rmax
}

void FENEBond::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / (1.0_kJmol / std::pow(1.0_angstrom, 2))}, {"rmax", std::sqrt(rmax_squared) / 1.0_angstrom}};
}

std::string FENEBond::name() const { return "fene"; }

void FENEBond::setEnergyFunction(const ParticleVector &p) {
    energyFunc = [&](Geometry::DistanceFunction dist) {
        double r_squared = dist(p[index[0]].pos, p[index[1]].pos).squaredNorm();
        return (r_squared >= rmax_squared) ? pc::infty
                                           : -k_half * rmax_squared * std::log(1 - r_squared / rmax_squared);
    };
}

FENEWCABond::FENEWCABond(double k, double rmax, double epsilon, double sigma, const std::vector<int> &index)
    : StretchData(index), k_half(k / 2), rmax_squared(rmax * rmax), epsilon(epsilon), sigma_squared(sigma * sigma) {}

std::shared_ptr<BondData> FENEWCABond::clone() const { return std::make_shared<FENEWCABond>(*this); }

BondData::Variant FENEWCABond::type() const { return BondData::FENEWCA; }

void FENEWCABond::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2;
    rmax_squared = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2);
    epsilon = j.at("eps").get<double>() * 1.0_kJmol;
    sigma_squared = std::pow(j.at("sigma").get<double>() * 1.0_angstrom, 2);
}

void FENEWCABond::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / (1.0_kJmol / std::pow(1.0_angstrom, 2))},
         {"rmax", std::sqrt(rmax_squared) / 1.0_angstrom},
         {"eps", epsilon / 1.0_kJmol},
         {"sigma", std::sqrt(sigma_squared) / 1.0_angstrom}};
}

std::string FENEWCABond::name() const { return "fene+wca"; }
void FENEWCABond::setEnergyFunction(const ParticleVector &p) {
    energyFunc = [&](Geometry::DistanceFunction dist) {
        double r_squared = dist(p[index[0]].pos, p[index[1]].pos).squaredNorm();
        double wca = 0;
        double x = sigma_squared;
        if (r_squared <= x * 1.2599210498948732) {
            x = x / r_squared;
            x = x * x * x;
            wca = epsilon * (x * x - x + 0.25);
        }
        return (r_squared > rmax_squared) ? pc::infty
                                          : -k_half * rmax_squared * std::log(1 - r_squared / rmax_squared) + wca;
    };
}

void HarmonicTorsion::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_rad, 2) / 2;
    aeq = j.at("aeq").get<double>() * 1.0_deg;
}

void HarmonicTorsion::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / (1.0_kJmol / std::pow(1.0_rad, 2))}, {"aeq", aeq / 1.0_deg}};
    _roundjson(j, 6);
}

HarmonicTorsion::HarmonicTorsion(double k, double aeq, const std::vector<int> &index)
    : TorsionData(index), k_half(k / 2), aeq(aeq) {}

BondData::Variant HarmonicTorsion::type() const { return BondData::HARMONIC_TORSION; }

std::string HarmonicTorsion::name() const { return "harmonic_torsion"; }

std::shared_ptr<BondData> HarmonicTorsion::clone() const { return std::make_shared<HarmonicTorsion>(*this); }

void HarmonicTorsion::setEnergyFunction(const ParticleVector &p) {
    energyFunc = [&](Geometry::DistanceFunction dist) {
        Point ray1 = dist(p[index[0]].pos, p[index[1]].pos);
        Point ray2 = dist(p[index[2]].pos, p[index[1]].pos);
        double angle = std::acos(ray1.dot(ray2) / ray1.norm() / ray2.norm());
        return k_half * (angle - aeq) * (angle - aeq);
    };
}

void GromosTorsion::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / 2;   // k
    cos_aeq = cos(j.at("aeq").get<double>() * 1.0_deg); // cos(angle)
}

void GromosTorsion::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / 1.0_kJmol}, {"aeq", std::acos(cos_aeq) / 1.0_deg}};
}

GromosTorsion::GromosTorsion(double k, double cos_aeq, const std::vector<int> &index)
    : TorsionData(index), k_half(k / 2), cos_aeq(cos_aeq) {}

BondData::Variant GromosTorsion::type() const { return BondData::GROMOS_TORSION; }

std::string GromosTorsion::name() const { return "gromos_torsion"; }

std::shared_ptr<BondData> GromosTorsion::clone() const { return std::make_shared<GromosTorsion>(*this); }

void GromosTorsion::setEnergyFunction(const ParticleVector &p) {
    energyFunc = [&](Geometry::DistanceFunction dist) {
        Point ray1 = dist(p[index[0]].pos, p[index[1]].pos);
        Point ray2 = dist(p[index[2]].pos, p[index[1]].pos);
        double dcos = cos_aeq - ray1.dot(ray2) / (ray1.norm() * ray2.norm());
        return k_half * dcos * dcos;
    };
}

int PeriodicDihedral::numindex() const { return 4; }

std::shared_ptr<BondData> PeriodicDihedral::clone() const { return std::make_shared<PeriodicDihedral>(*this); }

void PeriodicDihedral::from_json(const Faunus::json &j) {
    k = j.at("k").get<double>() * 1.0_kJmol;   // k
    n = j.at("n").get<double>();               // multiplicity/periodicity n
    phi = j.at("phi").get<double>() * 1.0_deg; // angle
}

void PeriodicDihedral::to_json(Faunus::json &j) const {
    j = {{"k", k / 1.0_kJmol}, {"n", n}, {"phi", phi / 1.0_deg}};
}

PeriodicDihedral::PeriodicDihedral(double k, double phi, double n, const std::vector<int> &index)
    : BondData(index), k(k), phi(phi), n(n) {}

BondData::Variant PeriodicDihedral::type() const { return BondData::PERIODIC_DIHEDRAL; }

std::string PeriodicDihedral::name() const { return "periodic_dihedral"; }

void PeriodicDihedral::setEnergyFunction(const ParticleVector &p) {
    energyFunc = [&](Geometry::DistanceFunction dist) {
        Point vec1 = dist(p[index[1]].pos, p[index[0]].pos);
        Point vec2 = dist(p[index[2]].pos, p[index[1]].pos);
        Point vec3 = dist(p[index[3]].pos, p[index[2]].pos);
        Point norm1 = vec1.cross(vec2);
        Point norm2 = vec2.cross(vec3);
        // atan2( [v1×v2]×[v2×v3]⋅[v2/|v2|], [v1×v2]⋅[v2×v3] )
        double angle = atan2((norm1.cross(norm2)).dot(vec2) / vec2.norm(), norm1.dot(norm2));
        return k * (1 + cos(n * angle - phi));
    };
}

StretchData::StretchData(const std::vector<int> &index) : BondData(index) {}
int StretchData::numindex() const { return 2; }

TorsionData::TorsionData(const std::vector<int> &index) : BondData(index) {}
int TorsionData::numindex() const { return 3; }

TEST_SUITE_BEGIN("Bonds");

TEST_CASE("[Faunus] BondData") {
    using doctest::Approx;

    ParticleVector p_4a(2, Particle());
    p_4a[0].pos = {2.0, 1.0, -2.0};
    p_4a[1].pos = {2.0, 1.0, 2.0};
    ParticleVector p_60deg_4a(3, Particle());
    p_60deg_4a[0].pos = {1.0, 1.0 + std::sqrt(3), 4.0};
    p_60deg_4a[1].pos = {1.0, 1.0, 1.0};
    p_60deg_4a[2].pos = {1.0, 5.0, 1.0};

    Geometry::DistanceFunction distance = [](const Point &a, const Point &b) -> Point { return b - a; };
    Geometry::DistanceFunction distance_3a = [](const Point &, const Point &) -> Point { return {0, 3, 0}; };
    Geometry::DistanceFunction distance_5a = [](const Point &, const Point &) -> Point { return {0, 3, 4}; };

    typedef std::shared_ptr<BondData> BondDataPtr;
    BondDataPtr bond_ptr;

    SUBCASE("HarmonicBond") {
        SUBCASE("HarmonicBond Energy") {
            HarmonicBond bond(100.0, 5.0, {0, 1});
            bond.setEnergyFunction(p_4a);
            CHECK_EQ(bond.energyFunc(distance_5a), Approx(0));
            CHECK_EQ(bond.energyFunc(distance_3a), Approx(200));
            CHECK_EQ(bond.energyFunc(distance), Approx(50));
        }
        SUBCASE("HarmonicBond Force") {
            HarmonicBond bond(100.0, 4, {0, 1});
            bond.setEnergyFunction(p_4a);
            auto forces = bond.forceFunc(distance_3a);
            CHECK(forces.size() == 2);
            CHECK(forces[0].x() == Approx(0));
            CHECK(forces[0].y() == Approx(100));
            CHECK(forces[0].y() == Approx(-forces[1].y()));
            CHECK(forces[0].z() == Approx(0));
        }
        SUBCASE("HarmonicBond JSON") {
            json j = R"({"harmonic": {"index":[1,2], "k":10.0, "req":2.0}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<HarmonicBond>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energyFunc(distance), Approx(10.0_kJmol / 2 * 4));
        }
        SUBCASE("HarmonicBond JSON Invalid") {
            CHECK_NOTHROW((R"({"harmonic": {"index":[0,9], "k":0.5, "req":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"harmoNIC": {"index":[2,3], "k":0.5, "req":2.1}})"_json)
                             .get<BondDataPtr>()); // exact match required
            CHECK_THROWS(
                (R"({"harmonic": {"index":[2], "k":0.5, "req":2.1}})"_json).get<BondDataPtr>());   // 2 atom indices
            CHECK_THROWS((R"({"harmonic": {"index":[2,3], "req":2.1}})"_json).get<BondDataPtr>()); // k missing
            CHECK_THROWS((R"({"harmonic": {"index":[2,3], "k":2.1}})"_json).get<BondDataPtr>());   // req missing
        }
    }

    SUBCASE("FENEBond") {
        SUBCASE("FENEBond Energy") {
            FENEBond bond(100.0, 5.0, {0, 1});
            bond.setEnergyFunction(p_4a);
            CHECK_EQ(bond.energyFunc(distance_5a), pc::infty);
            CHECK_EQ(bond.energyFunc(distance_3a), Approx(557.86));
            CHECK_EQ(bond.energyFunc(distance), Approx(1277.06));
        }
        SUBCASE("FENEBond JSON") {
            json j = R"({"fene": {"index":[1,2], "k":8, "rmax":6.0 }})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<FENEBond>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energyFunc(distance), Approx(84.641_kJmol));
        }
        SUBCASE("FENEBond JSON Invalid") {
            CHECK_NOTHROW((R"({"fene": {"index":[0,9], "k":0.5, "rmax":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"FENE": {"index":[0,9], "k":0.5, "rmax":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene": {"index":[2,3,4], "k":1, "rmax":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene": {"index":[2,3], "rmax":2.1}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene": {"index":[2,3], "k":1}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("FENEWCABond") {
        SUBCASE("FENEWCABond Energy") {
            FENEWCABond bond(100.0, 5.0, 20.0, 3.2, {0, 1});
            bond.setEnergyFunction(p_4a);
            CHECK_EQ(bond.energyFunc(distance_5a), pc::infty);
            CHECK_EQ(bond.energyFunc(distance_3a), Approx(557.86 + 18.931));
            CHECK_EQ(bond.energyFunc(distance), Approx(1277.06));
        }
        SUBCASE("FENEWCABond JSON") {
            json j = R"({"fene+wca": {"index":[1,2], "k":8, "rmax":6.0, "eps":3.5, "sigma":4.5}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<FENEWCABond>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energyFunc(distance), Approx(92.805_kJmol));
        }
        SUBCASE("FENEWCABond JSON Invalid") {
            CHECK_NOTHROW(
                (R"({"fene+wca": {"index":[0,9], "k":1, "rmax":2.1, "eps":2.48, "sigma":2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3,4], "k":1, "rmax":2.1, "eps":2.48, "sigma":2}})"_json)
                             .get<BondDataPtr>());
            CHECK_THROWS(
                (R"({"fene+wca": {"index":[2,3], "rmax":2.1, "eps":2.48, "sigma":2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3], "k":1, "eps":2.48, "sigma":2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3], "k":1, "rmax":2.1, "eps":2.48}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"fene+wca": {"index":[2,3], "k":1, "rmax":2.1, "sigma":2}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("HarmonicTorsion") {
        SUBCASE("HarmonicTorsion Energy") {
            HarmonicTorsion bond(100.0, 45.0_deg, {0, 1, 2});
            bond.setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond.energyFunc(distance), Approx(100.0 / 2 * std::pow(15.0_deg, 2)));
        }
        SUBCASE("HarmonicTorsion JSON") {
            json j = R"({"harmonic_torsion": {"index":[0,1,2], "k":0.5, "aeq":65}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<HarmonicTorsion>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energyFunc(distance), Approx(0.5_kJmol / 2 * std::pow(5.0_deg, 2)));
        }
        SUBCASE("HarmonicTorsion JSON Invalid") {
            CHECK_NOTHROW((R"({"harmonic_torsion": {"index":[0,1,9], "k":0.5, "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"harmonic_torsion": {"index":[2], "k":0.5, "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"harmonic_torsion": {"index":[0,1,2], "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"harmonic_torsion": {"index":[0,1,3], "k":0.5}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("GromosTorsion") {
        SUBCASE("GromosTorsion Energy") {
            GromosTorsion bond(100.0, cos(45.0_deg), {0, 1, 2});
            bond.setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond.energyFunc(distance), Approx(100.0 / 2 * std::pow(cos(60.0_deg) - cos(45.0_deg), 2)));
        }
        SUBCASE("GromosTorsion JSON") {
            json j = R"({"gromos_torsion": {"index":[0,1,2], "k":0.5, "aeq":65}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<GromosTorsion>(bond_ptr)->setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond_ptr->energyFunc(distance),
                     Approx(0.5_kJmol / 2 * std::pow(cos(60.0_deg) - cos(65.0_deg), 2)));
        }
        SUBCASE("GromosTorsion JSON Invalid") {
            CHECK_NOTHROW((R"({"gromos_torsion": {"index":[0,1,9], "k":0.5, "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"gromos_torsion": {"index":[2], "k":0.5, "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"gromos_torsion": {"index":[0,1,2], "aeq":30.5}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"gromos_torsion": {"index":[0,1,3], "k":0.5}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("PeriodicDihedral") {
        ParticleVector p_45deg(4, Particle());
        p_45deg[1].pos = {0.0, 0.0, 0.0};
        p_45deg[2].pos = {0.0, 0.0, 2.0};
        p_45deg[0].pos = {5.0, 0.0, 0.0};
        p_45deg[3].pos = {10.0, 10.0, 2.0};

        ParticleVector p_90deg(p_45deg);
        p_90deg[3].pos[0] *= 0;
        ParticleVector p_60deg(p_45deg);
        p_60deg[3].pos[1] *= std::sqrt(3);
        ParticleVector p_120deg(p_60deg);
        p_120deg[3].pos[0] *= -1;

        SUBCASE("PeriodicDihedral Energy") {
            PeriodicDihedral bond(100.0, 0.0_deg, 3, {0, 1, 2, 3});
            bond.setEnergyFunction(p_120deg);
            CHECK_EQ(bond.energyFunc(distance), Approx(200.0));
            bond.setEnergyFunction(p_60deg);
            CHECK_EQ(bond.energyFunc(distance), Approx(0.0));
            bond.setEnergyFunction(p_90deg);
            CHECK_EQ(bond.energyFunc(distance), Approx(100.0));
        }
        SUBCASE("PeriodicDihedral JSON") {
            json j = R"({"periodic_dihedral": {"index":[0,1,2,3], "k":10, "phi":0.0, "n": 3}})"_json;
            bond_ptr = j;
            CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<PeriodicDihedral>(bond_ptr)->setEnergyFunction(p_90deg);
            CHECK_EQ(bond_ptr->energyFunc(distance), Approx(10.0_kJmol));
        }
        SUBCASE("PeriodicDihedral JSON Invalid") {
            CHECK_NOTHROW(
                (R"({"periodic_dihedral": {"index":[0,1,2,9], "k":0.5, "phi":2.1, "n": 2}})"_json).get<BondDataPtr>());
            CHECK_THROWS(
                (R"({"periodic_dihedral": {"index":[0,1,2], "k":0.5, "phi":2.1, "n": 2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"periodic_dihedral": {"index":[0,1,2,3], "phi":2.1, "n": 2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"periodic_dihedral": {"index":[0,1,2,3], "k":0.5, "n": 2}})"_json).get<BondDataPtr>());
            CHECK_THROWS((R"({"periodic_dihedral": {"index":[0,1,2,3], "k":0.5, "phi":2.1}})"_json).get<BondDataPtr>());
        }
    }

    SUBCASE("Find") {
        BasePointerVector<BondData> bonds;
        bonds.emplace_back<FENEBond>(1.0, 2.1, std::vector<int>{2, 3});
        bonds.emplace_back<HarmonicBond>(1.0, 2.1, std::vector<int>{2, 3});
        auto harmonic_bonds = bonds.find<HarmonicBond>();
        CHECK(harmonic_bonds.size() == 1);
        CHECK(harmonic_bonds.front()->type() == BondData::HARMONIC);
        CHECK(harmonic_bonds.front() == bonds.back()); // harmonic_bonds should contain references to bonds
    }
}
TEST_SUITE_END();

} // namespace Faunus::Potential
