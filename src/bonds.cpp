#include <doctest/doctest.h>
#include "bonds.h"
#include "geometry.h"
#include "auxiliary.h"

namespace Faunus::Potential {

void to_json(Faunus::json& j, const std::shared_ptr<const BondData>& bond) { to_json(j, *bond); }

void to_json(Faunus::json& j, const BondData& bond) {
    json val;
    bond.to_json(val);
    val["index"] = bond.indices;
    j = {{bond.type(), val}};
}

void from_json(const json& j, std::shared_ptr<BondData>& bond) {
    try {
        const auto& [bondtype, parameters] = jsonSingleItem(j);
        const BondData::Variant variant = json(bondtype);
        switch (variant) {
        case BondData::HARMONIC:
            bond = std::make_shared<HarmonicBond>();
            break;
        case BondData::FENE:
            bond = std::make_shared<FENEBond>();
            break;
        case BondData::FENEWCA:
            bond = std::make_shared<FENEWCABond>();
            break;
        case BondData::HARMONIC_TORSION:
            bond = std::make_shared<HarmonicTorsion>();
            break;
        case BondData::GROMOS_TORSION:
            bond = std::make_shared<GromosTorsion>();
            break;
        case BondData::PERIODIC_DIHEDRAL:
            bond = std::make_shared<PeriodicDihedral>();
            break;
        default:
            throw ConfigurationError("invalid bondtype '{}'", bondtype);
        }
        try {
            bond->from_json(parameters);
            bond->indices = parameters.at("index").get<decltype(bond->indices)>();
            if (bond->indices.size() != bond->numindex()) {
                usageTip.pick(bondtype);
                throw ConfigurationError("exactly {} indices required", bond->numindex());
            }
        } catch (std::exception& e) {
            throw ConfigurationError("{} -> {}", bondtype, e.what());
        }
    } catch (std::exception& e) {
        throw ConfigurationError("bond making: {}", e.what()).attachJson(j);
    }
}

void BondData::shiftIndices(const int offset) {
    for (auto& i : indices) {
        i += offset;
    }
}

bool BondData::hasEnergyFunction() const { return energyFunc != nullptr; }

bool BondData::hasForceFunction() const { return forceFunc != nullptr; }

BondData::BondData(const std::vector<int>& indices) : indices(indices) {}

HarmonicBond::HarmonicBond(double k, double req, const std::vector<int>& indices)
    : StretchData(indices), half_force_constant(k / 2.0), equilibrium_distance(req) {}

void HarmonicBond::from_json(const Faunus::json& j) {
    half_force_constant = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2.0; // k
    equilibrium_distance = j.at("req").get<double>() * 1.0_angstrom;                             // req
}

void HarmonicBond::to_json(Faunus::json& j) const {
    j = {{"k", 2.0 * half_force_constant / 1.0_kJmol * 1.0_angstrom * 1.0_angstrom},
         {"req", equilibrium_distance / 1.0_angstrom}};
}

std::shared_ptr<BondData> HarmonicBond::clone() const { return std::make_shared<HarmonicBond>(*this); }

BondData::Variant HarmonicBond::type() const { return BondData::HARMONIC; }

/**
 * @param particles Particle vector to all particles in the system
 *
 * This sets both the `energyFunc` and `forceFunc` function objects
 * for calculating the potential energy and the forces on the
 * participating atoms
 */
void HarmonicBond::setEnergyFunction(const ParticleVector& particles) {
    energyFunc = [&](Geometry::DistanceFunction calculateDistance) { // potential energy functor
        const auto& particle1 = particles[indices[0]];
        const auto& particle2 = particles[indices[1]];
        const auto distance = equilibrium_distance - calculateDistance(particle1.pos, particle2.pos).norm();
        return half_force_constant * distance * distance; // kT/Å
    };
    forceFunc = [&](Geometry::DistanceFunction calculateDistance) { // force functor
        const auto& particle1 = particles[indices[0]];
        const auto& particle2 = particles[indices[1]];
        const auto distance_vec = calculateDistance(particle1.pos, particle2.pos);
        const auto distance = distance_vec.norm(); // distance between particles
        auto force = 2.0 * half_force_constant * (equilibrium_distance - distance) * distance_vec / distance;
        return std::vector<ParticleForce>({{indices[0], force}, {indices[1], -force}}); // force on both particles
    };
}

FENEBond::FENEBond(double k, double rmax, const std::vector<int>& indices)
    : StretchData(indices), half_force_constant(k / 2), max_squared_distance(rmax * rmax) {}

std::shared_ptr<BondData> FENEBond::clone() const { return std::make_shared<FENEBond>(*this); }

BondData::Variant FENEBond::type() const { return BondData::FENE; }

void FENEBond::from_json(const Faunus::json& j) {
    half_force_constant = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2; // k
    max_squared_distance = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2);             // rmax
}

void FENEBond::to_json(Faunus::json& j) const {
    j = {{"k", 2.0 * half_force_constant / (1.0_kJmol / std::pow(1.0_angstrom, 2))},
         {"rmax", std::sqrt(max_squared_distance) / 1.0_angstrom}};
}

void FENEBond::setEnergyFunction(const ParticleVector& particles) {
    energyFunc = [&](Geometry::DistanceFunction calcDistance) {
        const auto squared_distance = calcDistance(particles[indices[0]].pos, particles[indices[1]].pos).squaredNorm();
        if (squared_distance >= max_squared_distance) {
            return pc::infty;
        }
        return -half_force_constant * max_squared_distance * std::log(1.0 - squared_distance / max_squared_distance);
    };
}

FENEWCABond::FENEWCABond(double k, double rmax, double epsilon, double sigma, const std::vector<int>& indices)
    : StretchData(indices), half_force_constant(k / 2.0), max_distance_squared(rmax * rmax), epsilon(epsilon),
      sigma_squared(sigma * sigma) {}

std::shared_ptr<BondData> FENEWCABond::clone() const { return std::make_shared<FENEWCABond>(*this); }

BondData::Variant FENEWCABond::type() const { return BondData::FENEWCA; }

void FENEWCABond::from_json(const Faunus::json& j) {
    half_force_constant = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2.0;
    max_distance_squared = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2);
    epsilon = j.at("eps").get<double>() * 1.0_kJmol;
    sigma_squared = std::pow(j.at("sigma").get<double>() * 1.0_angstrom, 2);
}

void FENEWCABond::to_json(Faunus::json& j) const {
    j = {{"k", 2.0 * half_force_constant / (1.0_kJmol / std::pow(1.0_angstrom, 2))},
         {"rmax", std::sqrt(max_distance_squared) / 1.0_angstrom},
         {"eps", epsilon / 1.0_kJmol},
         {"sigma", std::sqrt(sigma_squared) / 1.0_angstrom}};
}

void FENEWCABond::setEnergyFunction(const ParticleVector& particles) {
    energyFunc = [&](Geometry::DistanceFunction calculateDistance) {
        double wca = 0.0;
        constexpr auto two_to_the_power_of_two_sixths = 1.2599210498948732; // 2^((1/6)^2)
        const auto squared_distance =
            calculateDistance(particles[indices[0]].pos, particles[indices[1]].pos).squaredNorm();
        if (squared_distance <= sigma_squared * two_to_the_power_of_two_sixths) {
            double sigma6 = sigma_squared / squared_distance;
            sigma6 = sigma6 * sigma6 * sigma6;
            wca = epsilon * (sigma6 * sigma6 - sigma6 + 0.25);
        }
        if (squared_distance > max_distance_squared) {
            return pc::infty;
        }
        return -half_force_constant * max_distance_squared * std::log(1.0 - squared_distance / max_distance_squared) +
               wca;
    };
}

void HarmonicTorsion::from_json(const Faunus::json& j) {
    half_force_constant = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_rad, 2) / 2.0;
    equilibrium_angle = j.at("aeq").get<double>() * 1.0_deg;
}

void HarmonicTorsion::to_json(Faunus::json& j) const {
    j = {{"k", 2 * half_force_constant / (1.0_kJmol / std::pow(1.0_rad, 2))}, {"aeq", equilibrium_angle / 1.0_deg}};
    _roundjson(j, 6);
}

HarmonicTorsion::HarmonicTorsion(double k, double aeq, const std::vector<int>& indices)
    : TorsionData(indices), half_force_constant(k / 2.0), equilibrium_angle(aeq) {}

BondData::Variant HarmonicTorsion::type() const { return BondData::HARMONIC_TORSION; }

std::shared_ptr<BondData> HarmonicTorsion::clone() const { return std::make_shared<HarmonicTorsion>(*this); }

void HarmonicTorsion::setEnergyFunction(const ParticleVector& particles) {
    energyFunc = [&](Geometry::DistanceFunction calculateDistance) {
        const auto ray1 = calculateDistance(particles[indices[0]].pos, particles[indices[1]].pos);
        const auto ray2 = calculateDistance(particles[indices[2]].pos, particles[indices[1]].pos);
        const auto angle = std::acos(ray1.dot(ray2) / (ray1.norm() * ray2.norm()));
        return half_force_constant * (angle - equilibrium_angle) * (angle - equilibrium_angle);
    };
}

void GromosTorsion::from_json(const Faunus::json& j) {
    half_force_constant = j.at("k").get<double>() * 1.0_kJmol / 2.0;     // k
    cosine_equilibrium_angle = cos(j.at("aeq").get<double>() * 1.0_deg); // cos(angle)
}

void GromosTorsion::to_json(Faunus::json& j) const {
    j = {{"k", 2.0 * half_force_constant / 1.0_kJmol}, {"aeq", std::acos(cosine_equilibrium_angle) / 1.0_deg}};
}

GromosTorsion::GromosTorsion(double k, double cos_aeq, const std::vector<int>& indices)
    : TorsionData(indices), half_force_constant(k / 2.0), cosine_equilibrium_angle(cos_aeq) {}

BondData::Variant GromosTorsion::type() const { return BondData::GROMOS_TORSION; }

std::shared_ptr<BondData> GromosTorsion::clone() const { return std::make_shared<GromosTorsion>(*this); }

void GromosTorsion::setEnergyFunction(const ParticleVector& particles) {
    energyFunc = [&](Geometry::DistanceFunction calculateDistance) {
        auto ray1 = calculateDistance(particles[indices[0]].pos, particles[indices[1]].pos);
        auto ray2 = calculateDistance(particles[indices[2]].pos, particles[indices[1]].pos);
        const auto x = cosine_equilibrium_angle - ray1.dot(ray2) / (ray1.norm() * ray2.norm());
        return half_force_constant * x * x;
    };
}

int PeriodicDihedral::numindex() const { return 4; }

std::shared_ptr<BondData> PeriodicDihedral::clone() const { return std::make_shared<PeriodicDihedral>(*this); }

void PeriodicDihedral::from_json(const Faunus::json& j) {
    force_constant = j.at("k").get<double>() * 1.0_kJmol; // k
    periodicity = j.at("n").get<double>();                // multiplicity/periodicity n
    dihedral_angle = j.at("phi").get<double>() * 1.0_deg; // angle
}

void PeriodicDihedral::to_json(Faunus::json& j) const {
    j = {{"k", force_constant / 1.0_kJmol}, {"n", periodicity}, {"phi", dihedral_angle / 1.0_deg}};
}

PeriodicDihedral::PeriodicDihedral(double k, double phi, double n, const std::vector<int>& indices)
    : BondData(indices), force_constant(k), dihedral_angle(phi), periodicity(n) {}

BondData::Variant PeriodicDihedral::type() const { return BondData::PERIODIC_DIHEDRAL; }

void PeriodicDihedral::setEnergyFunction(const ParticleVector& particles) {
    energyFunc = [&](Geometry::DistanceFunction dist) {
        auto v1 = dist(particles[indices[1]].pos, particles[indices[0]].pos);
        auto v2 = dist(particles[indices[2]].pos, particles[indices[1]].pos);
        auto v3 = dist(particles[indices[3]].pos, particles[indices[2]].pos);
        auto norm1 = v1.cross(v2);
        auto norm2 = v2.cross(v3);
        // atan2( [v1×v2]×[v2×v3]⋅[v2/|v2|], [v1×v2]⋅[v2×v3] )
        const auto angle = std::atan2((norm1.cross(norm2)).dot(v2) / v2.norm(), norm1.dot(norm2));
        return force_constant * (1.0 + std::cos(periodicity * angle - dihedral_angle));
    };
}

StretchData::StretchData(const std::vector<int>& indices) : BondData(indices) {}
int StretchData::numindex() const { return 2; }

TorsionData::TorsionData(const std::vector<int>& indices) : BondData(indices) {}
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

    Geometry::DistanceFunction distance = [](const Point& a, const Point& b) -> Point { return b - a; };
    Geometry::DistanceFunction distance_3a = [](const Point&, const Point&) -> Point { return {0, 3, 0}; };
    Geometry::DistanceFunction distance_5a = [](const Point&, const Point&) -> Point { return {0, 3, 4}; };

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
            CHECK(forces[0].second.x() == Approx(0));
            CHECK(forces[0].second.y() == Approx(100));
            CHECK(forces[0].second.y() == Approx(-forces[1].second.y()));
            CHECK(forces[0].second.z() == Approx(0));
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
