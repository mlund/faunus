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
        case BondData::HARMONIC_DIHEDRAL:
            bond = std::make_shared<HarmonicDihedral>();
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

HarmonicBond::HarmonicBond(double force_constant, double equilibrium_distance, const std::vector<int>& indices)
    : StretchData(indices), half_force_constant(force_constant / 2.0), equilibrium_distance(equilibrium_distance) {}

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
    energyFunc = [&](Geometry::DistanceFunction calculateDistance) -> double {
        const auto& particle1 = particles[indices[0]];
        const auto& particle2 = particles[indices[1]];
        const auto distance = equilibrium_distance - calculateDistance(particle1.pos, particle2.pos).norm();
        return half_force_constant * distance * distance; // kT/Å
    };
    forceFunc = [&](Geometry::DistanceFunction calculateDistance) -> std::vector<IndexAndForce> {
        const auto& particle1 = particles[indices[0]];
        const auto& particle2 = particles[indices[1]];
        const auto distance_vec = calculateDistance(particle1.pos, particle2.pos);
        const auto distance = distance_vec.norm(); // distance between particles
        auto force = 2.0 * half_force_constant * (equilibrium_distance - distance) * distance_vec / distance;
        return {{indices[0], force}, {indices[1], -force}};
    };
}

FENEBond::FENEBond(double force_constant, double max_distance, const std::vector<int>& indices)
    : StretchData(indices), half_force_constant(0.5 * force_constant),
      max_squared_distance(max_distance * max_distance) {}

std::shared_ptr<BondData> FENEBond::clone() const { return std::make_shared<FENEBond>(*this); }

BondData::Variant FENEBond::type() const { return BondData::FENE; }

void FENEBond::from_json(const Faunus::json& j) {
    half_force_constant = 0.5 * j.at("k").get<double>() * 1.0_kJmol / (1.0_angstrom * 1.0_angstrom);
    max_squared_distance = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2);
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
    forceFunc = [&](Geometry::DistanceFunction distance) -> std::vector<IndexAndForce> {
        const Point ba = distance(particles[indices[0]].pos, particles[indices[1]].pos); // b->a
        const auto squared_distance = ba.squaredNorm();
        if (squared_distance >= max_squared_distance) {
            throw std::runtime_error("Fene potential: Force undefined for distances greater than rmax.");
        };
        const auto force_magnitude =
            -2.0 * half_force_constant * ba.norm() / (1.0 - squared_distance / max_squared_distance);
        Point force0 = force_magnitude * ba.normalized();
        Point force1 = -force0;
        return {{indices[0], force0}, {indices[1], force1}};
    };
}

FENEWCABond::FENEWCABond(double force_constant, double max_distance, double epsilon, double sigma,
                         const std::vector<int>& indices)
    : StretchData(indices), half_force_constant(force_constant / 2.0),
      max_distance_squared(max_distance * max_distance), epsilon(epsilon), sigma_squared(sigma * sigma) {}

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
        constexpr auto two_to_the_power_of_two_sixths = 1.01944064370214482816981563263103378007648819; // 2^((1/6)^2)
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
    forceFunc = [&](Geometry::DistanceFunction calculateDistance) -> std::vector<IndexAndForce> {
        const Point ba = calculateDistance(particles[indices[0]].pos, particles[indices[1]].pos); // b->a
        const auto distance = ba.norm();
        const auto squared_distance = distance * distance;
        if (squared_distance >= max_distance_squared) {
            throw std::runtime_error("Fene+WCA potential: Force undefined for distances greater than rmax.");
        }
        double wca_force = 0.0;
        constexpr auto two_to_the_power_of_two_sixths = 1.01944064370214482816981563263103378007648819; // 2^((1/6)^2)
        if (squared_distance <= sigma_squared * two_to_the_power_of_two_sixths) {
            double sigma6 = sigma_squared / squared_distance;
            sigma6 = sigma6 * sigma6 * sigma6;
            wca_force = -24.0 * epsilon * (2.0 * sigma6 * sigma6 - sigma6) / distance;
        }
        const auto force_magnitude =
            -(2.0 * half_force_constant * distance / (1.0 - squared_distance / max_distance_squared) + wca_force);
        Point force0 = force_magnitude * ba.normalized();
        Point force1 = -force0;
        return {{indices[0], force0}, {indices[1], force1}};
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

HarmonicTorsion::HarmonicTorsion(double force_constant, double equilibrium_angle, const std::vector<int>& indices)
    : TorsionData(indices), half_force_constant(force_constant / 2.0), equilibrium_angle(equilibrium_angle) {}

BondData::Variant HarmonicTorsion::type() const { return BondData::HARMONIC_TORSION; }

std::shared_ptr<BondData> HarmonicTorsion::clone() const { return std::make_shared<HarmonicTorsion>(*this); }

void HarmonicTorsion::setEnergyFunction(const ParticleVector& particles) {
    energyFunc = [&](Geometry::DistanceFunction calculateDistance) -> double {
        const auto vec1 = calculateDistance(particles[indices[0]].pos, particles[indices[1]].pos).normalized();
        const auto vec2 = calculateDistance(particles[indices[2]].pos, particles[indices[1]].pos).normalized();
        const auto delta_angle = equilibrium_angle - std::acos(vec1.dot(vec2));
        return half_force_constant * delta_angle * delta_angle;
    };
    forceFunc = [&](Geometry::DistanceFunction calculateDistance) -> std::vector<IndexAndForce> {
        if constexpr (false) {
            // see https://github.com/OpenMD/OpenMD/blob/master/src/primitives/Bend.cpp
            Point vec1 = calculateDistance(particles[indices[0]].pos, particles[indices[1]].pos);
            Point vec2 = calculateDistance(particles[indices[2]].pos, particles[indices[1]].pos);
            const auto inverse_norm1 = 1.0 / vec1.norm();
            const auto inverse_norm2 = 1.0 / vec2.norm();
            vec1 *= inverse_norm1;
            vec2 *= inverse_norm2;
            const auto cosine_angle = vec1.dot(vec2);
            const auto inverse_sine_angle = 1.0 / std::sqrt(std::fabs(1.0 - cosine_angle * cosine_angle));
            const auto angle = std::acos(cosine_angle);
            const auto prefactor = 2.0 * half_force_constant * (angle - equilibrium_angle) * inverse_sine_angle;
            Point force0 = prefactor * inverse_norm1 * (vec2 - vec1 * cosine_angle);
            Point force2 = prefactor * inverse_norm2 * (vec1 - vec2 * cosine_angle);
            Point force1 = -(force0 + force2); // no net force
            return {{indices[0], force0}, {indices[1], force1}, {indices[2], force2}};
        } else { // @todo which is faster?
            const Point ba = calculateDistance(particles[indices[0]].pos, particles[indices[1]].pos); // b->a
            const Point bc = calculateDistance(particles[indices[2]].pos, particles[indices[1]].pos); // b->c
            const auto inverse_norm_ba = 1.0 / ba.norm();
            const auto inverse_norm_bc = 1.0 / bc.norm();
            const auto angle = std::acos(ba.dot(bc) * inverse_norm_ba * inverse_norm_bc);
            const auto force_magnitude = -2.0 * half_force_constant * (angle - equilibrium_angle);
            const Point plane_abc = ba.cross(bc).eval();
            Point force0 = (force_magnitude * inverse_norm_ba * ba.cross(plane_abc).normalized());
            Point force2 = (force_magnitude * inverse_norm_bc * -bc.cross(plane_abc).normalized());
            Point force1 = -(force0 + force2); // Newton's third law
            return {{indices[0], force0}, {indices[1], force1}, {indices[2], force2}};
        }
    };
}

void GromosTorsion::from_json(const Faunus::json& j) {
    half_force_constant = 0.5 * j.at("k").get<double>() * 1.0_kJmol;
    cosine_equilibrium_angle = std::cos(j.at("aeq").get<double>() * 1.0_deg);
}

void GromosTorsion::to_json(Faunus::json& j) const {
    j = {{"k", 2.0 * half_force_constant / 1.0_kJmol}, {"aeq", std::acos(cosine_equilibrium_angle) / 1.0_deg}};
}

GromosTorsion::GromosTorsion(double force_constant, double cosine_equilibrium_angle, const std::vector<int>& indices)
    : TorsionData(indices), half_force_constant(0.5 * force_constant),
      cosine_equilibrium_angle(cosine_equilibrium_angle) {}

BondData::Variant GromosTorsion::type() const { return BondData::GROMOS_TORSION; }

std::shared_ptr<BondData> GromosTorsion::clone() const { return std::make_shared<GromosTorsion>(*this); }

void GromosTorsion::setEnergyFunction(const ParticleVector& particles) {
    energyFunc = [&](Geometry::DistanceFunction calculateDistance) {
        auto vec1 = calculateDistance(particles[indices[0]].pos, particles[indices[1]].pos).normalized();
        auto vec2 = calculateDistance(particles[indices[2]].pos, particles[indices[1]].pos).normalized();
        const auto cosine_angle_displacement = cosine_equilibrium_angle - vec1.dot(vec2);
        return half_force_constant * cosine_angle_displacement * cosine_angle_displacement;
    };
    forceFunc = [&](Geometry::DistanceFunction distance) -> std::vector<IndexAndForce> {
        const Point ba = distance(particles[indices[0]].pos, particles[indices[1]].pos); // b->a
        const Point bc = distance(particles[indices[2]].pos, particles[indices[1]].pos); // b->c
        const auto inverse_norm_ba = 1.0 / ba.norm();
        const auto inverse_norm_bc = 1.0 / bc.norm();
        const auto cosine_angle = ba.dot(bc) * inverse_norm_ba * inverse_norm_bc;
        const auto angle = std::acos(cosine_angle);
        const auto force_magnitude =
            -2.0 * half_force_constant * std::sin(angle) * (cosine_angle - cosine_equilibrium_angle);
        const Point plane_abc = ba.cross(bc).eval();
        Point force0 = force_magnitude * inverse_norm_ba * ba.cross(plane_abc).normalized();
        Point force2 = force_magnitude * inverse_norm_bc * -bc.cross(plane_abc).normalized();
        Point force1 = -(force0 + force2); // Newton's third law
        return {{indices[0], force0}, {indices[1], force1}, {indices[2], force2}};
    };
}

int PeriodicDihedral::numindex() const { return 4; }

std::shared_ptr<BondData> PeriodicDihedral::clone() const { return std::make_shared<PeriodicDihedral>(*this); }

void PeriodicDihedral::from_json(const Faunus::json& j) {
    force_constant = j.at("k").get<double>() * 1.0_kJmol;
    periodicity = j.at("n").get<double>();
    phase_angle = j.at("phi").get<double>() * 1.0_deg;
}

void PeriodicDihedral::to_json(Faunus::json& j) const {
    j = {{"k", force_constant / 1.0_kJmol}, {"n", periodicity}, {"phi", phase_angle / 1.0_deg}};
}

PeriodicDihedral::PeriodicDihedral(double force_constant, double phase_angle, double periodicity,
                                   const std::vector<int>& indices)
    : BondData(indices), force_constant(force_constant), phase_angle(phase_angle), periodicity(periodicity) {}

BondData::Variant PeriodicDihedral::type() const { return BondData::PERIODIC_DIHEDRAL; }

void PeriodicDihedral::setEnergyFunction(const ParticleVector& particles) {
    // Torsion on the form a(0) - b(1) - c(2) - d(3)
    energyFunc = [&](Geometry::DistanceFunction distance) {
        auto ab = distance(particles[indices[1]].pos, particles[indices[0]].pos); // a->b
        auto bc = distance(particles[indices[2]].pos, particles[indices[1]].pos); // b->c
        auto cd = distance(particles[indices[3]].pos, particles[indices[2]].pos); // c->d
        auto normal_abc = ab.cross(bc).eval();                                    // ab x bc
        auto normal_bcd = bc.cross(cd).eval();                                    // bc x cd
        // atan2( [ab×bc]×[bc×cd]⋅[bc/|bc|], [ab×bc]⋅[bc×cd] )
        const auto dihedral_angle =
            std::atan2((normal_abc.cross(normal_bcd)).dot(bc) / bc.norm(), normal_abc.dot(normal_bcd));
        return force_constant * (1.0 + std::cos(periodicity * dihedral_angle - phase_angle));
    };
    forceFunc = [&](Geometry::DistanceFunction distance) -> std::vector<IndexAndForce> {
        auto ab = distance(particles[indices[1]].pos, particles[indices[0]].pos); // a->b
        auto bc = distance(particles[indices[2]].pos, particles[indices[1]].pos); // b->c
        auto cd = distance(particles[indices[3]].pos, particles[indices[2]].pos); // c->d
        auto normal_abc = ab.cross(bc).eval();
        auto normal_bcd = bc.cross(cd).eval();
        const auto dihedral_angle =
            std::atan2((normal_abc.cross(normal_bcd)).dot(bc) / bc.norm(), normal_abc.dot(normal_bcd));

        // Calculation of the energy derivative with respect to the dihedral angle.
        const auto force_magnitude =
            periodicity * force_constant * std::sin(periodicity * dihedral_angle - phase_angle);

        // Calculation of the dihedral angle derivative with respect to the position vector.
        const auto inverse_norm_ab = 1.0 / ab.norm();
        const auto inverse_norm_bc = 1.0 / bc.norm();
        const auto inverse_norm_cd = 1.0 / cd.norm();
        const auto angle_abc = std::acos(-ab.dot(bc) * inverse_norm_ab * inverse_norm_bc);
        const auto angle_bcd = std::acos(-bc.dot(cd) * inverse_norm_bc * inverse_norm_cd);
        const auto theta_a_derivative = inverse_norm_ab / std::sin(angle_abc);
        const auto theta_d_derivative = inverse_norm_cd / std::sin(angle_bcd);

        // Calculation of directional vectors on particle a and d.
        const Point ortho_normalized_abc = -normal_abc.normalized();   // normalized vector orthogonal to the plane abc.
        const Point ortho_normalized_bcd = cd.cross(-bc).normalized(); // normalized vector orthogonal to the plane bcd.

        // Calculation of forces on particle a and d
        Point force_a = force_magnitude * ortho_normalized_abc * theta_a_derivative;
        Point force_d = force_magnitude * ortho_normalized_bcd * theta_d_derivative;

        // Calculation of force and associated vectors for atom c.
        const Point bc_midpoint = 0.5 * bc;
        const Point torque_c = -(bc_midpoint.cross(force_d) + 0.5 * cd.cross(force_d) - 0.5 * ab.cross(force_a));
        Point force_c = torque_c.cross(bc_midpoint) / bc_midpoint.squaredNorm();
        Point force_b = -(force_a + force_c + force_d); // Newton's third law for force on atom b.
        return {{indices[0], force_a}, {indices[1], force_b}, {indices[2], force_c}, {indices[3], force_d}};
    };
}

void HarmonicDihedral::from_json(const Faunus::json& j) {
    half_force_constant = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_rad, 2) / 2.0;
    equilibrium_dihedral = j.at("deq").get<double>() * 1.0_deg;
}

void HarmonicDihedral::to_json(Faunus::json& j) const {
    j = {{"k", 2 * half_force_constant / (1.0_kJmol / std::pow(1.0_rad, 2))}, {"deq", equilibrium_dihedral / 1.0_deg}};
    _roundjson(j, 6);
}

int HarmonicDihedral::numindex() const { return 4; }

HarmonicDihedral::HarmonicDihedral(double force_constant, double equilibrium_dihedral, const std::vector<int>& indices)
    : BondData(indices), half_force_constant(force_constant / 2.0), equilibrium_dihedral(equilibrium_dihedral) {}

BondData::Variant HarmonicDihedral::type() const { return BondData::HARMONIC_DIHEDRAL; }

std::shared_ptr<BondData> HarmonicDihedral::clone() const { return std::make_shared<HarmonicDihedral>(*this); }

void HarmonicDihedral::setEnergyFunction(const ParticleVector& particles) {
    // Torsion on the form a(0) - b(1) - c(2) - d(3)
    energyFunc = [&](Geometry::DistanceFunction distance) {
        auto ab = distance(particles[indices[1]].pos, particles[indices[0]].pos); // a->b
        auto bc = distance(particles[indices[2]].pos, particles[indices[1]].pos); // b->c
        auto cd = distance(particles[indices[3]].pos, particles[indices[2]].pos); // c->d
        auto normal_abc = ab.cross(bc).eval();                                    // ab x bc
        auto normal_bcd = bc.cross(cd).eval();                                    // bc x cd
        // atan2( [ab×bc]×[bc×cd]⋅[bc/|bc|], [ab×bc]⋅[bc×cd] )
        const auto dihedral_angle =
            std::atan2((normal_abc.cross(normal_bcd)).dot(bc) / bc.norm(), normal_abc.dot(normal_bcd));
        const auto delta_dihedral = equilibrium_dihedral - dihedral_angle;
        return half_force_constant * delta_dihedral * delta_dihedral;
    };
    forceFunc = [&](Geometry::DistanceFunction distance) -> std::vector<IndexAndForce> {
        auto ab = distance(particles[indices[1]].pos, particles[indices[0]].pos); // a->b
        auto bc = distance(particles[indices[2]].pos, particles[indices[1]].pos); // b->c
        auto cd = distance(particles[indices[3]].pos, particles[indices[2]].pos); // c->d
        auto normal_abc = ab.cross(bc).eval();
        auto normal_bcd = bc.cross(cd).eval();
        const auto dihedral_angle =
            std::atan2((normal_abc.cross(normal_bcd)).dot(bc) / bc.norm(), normal_abc.dot(normal_bcd));

        // Calculation of the energy derivative with respect to the dihedral angle.
        const auto force_magnitude = -2.0 * half_force_constant * (dihedral_angle - equilibrium_dihedral);

        // Calculation of the dihedral angle derivative with respect to the position vector.
        const auto inverse_norm_ab = 1.0 / ab.norm();
        const auto inverse_norm_bc = 1.0 / bc.norm();
        const auto inverse_norm_cd = 1.0 / cd.norm();
        const auto angle_abc = std::acos(-ab.dot(bc) * inverse_norm_ab * inverse_norm_bc);
        const auto angle_bcd = std::acos(-bc.dot(cd) * inverse_norm_bc * inverse_norm_cd);
        const auto theta_a_derivative = inverse_norm_ab / std::sin(angle_abc);
        const auto theta_d_derivative = inverse_norm_cd / std::sin(angle_bcd);

        // Calculation of directional vectors on particle a and d.
        const Point ortho_normalized_abc = -normal_abc.normalized();   // normalized vector orthogonal to the plane abc.
        const Point ortho_normalized_bcd = cd.cross(-bc).normalized(); // normalized vector orthogonal to the plane bcd.

        // Calculation of forces on particle a and d
        Point force_a = force_magnitude * ortho_normalized_abc * theta_a_derivative;
        Point force_d = force_magnitude * ortho_normalized_bcd * theta_d_derivative;

        // Calculation of force and associated vectors for atom c.
        const Point bc_midpoint = 0.5 * bc;
        const Point torque_c = -(bc_midpoint.cross(force_d) + 0.5 * cd.cross(force_d) - 0.5 * ab.cross(force_a));
        Point force_c = torque_c.cross(bc_midpoint) / bc_midpoint.squaredNorm();
        Point force_b = -(force_a + force_c + force_d); // Newton's third law for force on atom b.
        return {{indices[0], force_a}, {indices[1], force_b}, {indices[2], force_c}, {indices[3], force_d}};
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
    ParticleVector p_90deg_4a(3, Particle());
    p_90deg_4a[0].pos = {0.0, 1.0, 0.0};
    p_90deg_4a[1].pos = {0.0, 0.0, 0.0};
    p_90deg_4a[2].pos = {1.0, 0.0, 0.0};

    Geometry::DistanceFunction distance = [](auto& a, auto& b) -> Point { return a - b; };
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
        SUBCASE("FENEBond Force") {
            FENEBond bond(100.0, 5.0, {0, 1});
            bond.setEnergyFunction(p_4a);
            auto forces = bond.forceFunc(distance_3a);
            CHECK(forces.size() == 2);
            CHECK(forces[0].first == 0);
            CHECK(forces[1].first == 1);
            CHECK(forces[0].second.x() == Approx(0.0));
            CHECK(forces[0].second.y() == Approx(-468.75));
            CHECK(forces[0].second.z() == Approx(0.0));
            CHECK(forces[1].second.x() == Approx(0.0));
            CHECK(forces[1].second.y() == Approx(468.75));
            CHECK(forces[1].second.z() == Approx(0.0));
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
        SUBCASE("FENEWCABond Force") {
            FENEWCABond bond(100, 5.0, 20.0, 3.2, {0, 1});
            bond.setEnergyFunction(p_4a);
            auto forces = bond.forceFunc(distance_3a);
            CHECK(forces.size() == 2);
            CHECK(forces[0].first == 0);
            CHECK(forces[1].first == 1);
            CHECK(forces[0].second.x() == Approx(0.0));
            CHECK(forces[0].second.y() == Approx(-10.1974323155));
            CHECK(forces[0].second.z() == Approx(0.0));
            CHECK(forces[1].second.x() == Approx(0.0));
            CHECK(forces[1].second.y() == Approx(10.1974323155));
            CHECK(forces[1].second.z() == Approx(0.0));
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
        SUBCASE("HarmonicTorsion Force") {
            HarmonicTorsion bond(1, 45.0_deg, {0, 1, 2});
            bond.setEnergyFunction(p_90deg_4a);
            auto forces = bond.forceFunc(distance);
            CHECK(forces.size() == 3);
            CHECK(forces[0].first == 0);
            CHECK(forces[1].first == 1);
            CHECK(forces[2].first == 2);
            CHECK(forces[0].second.x() == Approx(0.78539816));
            CHECK(forces[0].second.y() == Approx(0.0));
            CHECK(forces[0].second.z() == Approx(0.0));
            CHECK(forces[1].second.x() == Approx(-0.78539816));
            CHECK(forces[1].second.y() == Approx(-0.78539816));
            CHECK(forces[1].second.z() == Approx(0.0));
            CHECK(forces[2].second.x() == Approx(0.0));
            CHECK(forces[2].second.y() == Approx(0.78539816));
            CHECK(forces[2].second.z() == Approx(0.0));
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
            GromosTorsion bond(100.0, std::cos(45.0_deg), {0, 1, 2});
            bond.setEnergyFunction(p_60deg_4a);
            CHECK_EQ(bond.energyFunc(distance), Approx(100.0 / 2 * std::pow(cos(60.0_deg) - cos(45.0_deg), 2)));
        }
        SUBCASE("GromosTorsion Force") {
            GromosTorsion bond(100.0, std::cos(45.0_deg), {0, 1, 2});
            bond.setEnergyFunction(p_90deg_4a);
            auto forces = bond.forceFunc(distance);
                CHECK(forces.size() == 3);
                CHECK(forces[0].first == 0);
                CHECK(forces[1].first == 1);
                CHECK(forces[2].first == 2);
                CHECK(forces[0].second.x() == Approx(-70.7106781187));
                CHECK(forces[0].second.y() == Approx(0.0));
                CHECK(forces[0].second.z() == Approx(0.0));
                CHECK(forces[1].second.x() == Approx(70.7106781187));
                CHECK(forces[1].second.y() == Approx(70.7106781187));
                CHECK(forces[1].second.z() == Approx(0.0));
                CHECK(forces[2].second.x() == Approx(0.0));
                CHECK(forces[2].second.y() == Approx(-70.7106781187));
                CHECK(forces[2].second.z() == Approx(0.0));
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
        ParticleVector p_45deg(4);
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
        SUBCASE("PeriodicDihedral Forces") {
            PeriodicDihedral bond(100.0, 0.0_deg, 3, {0, 1, 2, 3});
            bond.setEnergyFunction(p_90deg);
            auto forces = bond.forceFunc(distance);
            CHECK(forces.size() == 4);
            CHECK(forces[0].first == 0);
            CHECK(forces[1].first == 1);
            CHECK(forces[2].first == 2);
            CHECK(forces[3].first == 3);
            CHECK(forces[0].second.x() == Approx(0));
            CHECK(forces[0].second.y() == Approx(60));
            CHECK(forces[0].second.z() == Approx(0));
            CHECK(forces[1].second.x() == Approx(0));
            CHECK(forces[1].second.y() == Approx(-60));
            CHECK(forces[1].second.z() == Approx(0));
            CHECK(forces[2].second.x() == Approx(-30));
            CHECK(forces[2].second.y() == Approx(0));
            CHECK(forces[2].second.z() == Approx(0));
            CHECK(forces[3].second.x() == Approx(30));
            CHECK(forces[3].second.y() == Approx(0));
            CHECK(forces[3].second.z() == Approx(0));
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

        SUBCASE("HarmonicDihedral") {
        ParticleVector p_45deg(4);
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

            SUBCASE("HarmonicDihedral Energy") {
            HarmonicDihedral bond(100.0, 90.0_deg, {0, 1, 2, 3});
            bond.setEnergyFunction(p_120deg);
                CHECK_EQ(bond.energyFunc(distance), Approx(13.7077838904));
            bond.setEnergyFunction(p_60deg);
                CHECK_EQ(bond.energyFunc(distance), Approx(13.7077838904));
            bond.setEnergyFunction(p_90deg);
                CHECK_EQ(bond.energyFunc(distance), Approx(0.0));
        }
            SUBCASE("HarmonicDihedral Forces") {
            HarmonicDihedral bond(100.0, 90.0_deg, {0, 1, 2, 3});
            bond.setEnergyFunction(p_120deg);
            auto forces = bond.forceFunc(distance);
                CHECK(forces.size() == 4);
                CHECK(forces[0].first == 0);
                CHECK(forces[1].first == 1);
                CHECK(forces[2].first == 2);
                CHECK(forces[3].first == 3);
                CHECK(forces[0].second.x() == Approx(0));
                CHECK(forces[0].second.y() == Approx(10.471975512));
                CHECK(forces[0].second.z() == Approx(0));
                CHECK(forces[1].second.x() == Approx(0));
                CHECK(forces[1].second.y() == Approx(-10.471975512));
                CHECK(forces[1].second.z() == Approx(0));
                CHECK(forces[2].second.x() == Approx(-2.2672492053));
                CHECK(forces[2].second.y() == Approx(-1.308996939));
                CHECK(forces[2].second.z() == Approx(0));
                CHECK(forces[3].second.x() == Approx(2.2672492053));
                CHECK(forces[3].second.y() == Approx(1.308996939));
                CHECK(forces[3].second.z() == Approx(0));
        }
            SUBCASE("HarmonicDihedral JSON") {
            json j = R"({"harmonic_dihedral": {"index":[0,1,2,3], "k":100.0, "deq":90}})"_json;
            bond_ptr = j;
                CHECK_EQ(json(bond_ptr), j);
            std::dynamic_pointer_cast<HarmonicDihedral>(bond_ptr)->setEnergyFunction(p_120deg);
                CHECK_EQ(bond_ptr->energyFunc(distance), Approx(100.0_kJmol / 2.0 * std::pow(30.0_deg, 2)));
        }
            SUBCASE("HarmonicDihedral JSON Invalid") {
                CHECK_NOTHROW(
                (R"({"harmonic_dihedral": {"index":[0,1,2,3], "k":0.5, "deq":90}})"_json).get<BondDataPtr>());
                CHECK_THROWS(
                (R"({"harmonic_dihedral": {"index":[0,1,2], "k":0.5, "deq":90}})"_json).get<BondDataPtr>());
                CHECK_THROWS((R"({"harmonic_dihedral": {"index":[0,1,2,3], "deq":90}})"_json).get<BondDataPtr>());
                CHECK_THROWS((R"({"harmonic_dihedral": {"index":[0,1,2,3], "k":0.5}})"_json).get<BondDataPtr>());
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
