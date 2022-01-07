#include "regions.h"
#include "space.h"
#include <range/v3/algorithm/any_of.hpp>
#include <iostream>

namespace Faunus {
namespace Region {

RegionBase::RegionBase(RegionType type) : type(type) {}

std::optional<double> RegionBase::volume() const { return std::nullopt; }

bool RegionBase::inside(const Particle& particle) const { return isInside(particle.pos); }

bool RegionBase::inside(const Group& group) const {
    if (use_group_mass_center) {
        if (auto mass_center = group.massCenter()) {
            return isInside(mass_center.value());
        }
    }
    return ranges::cpp20::any_of(group, [&](const Particle& particle) { return inside(particle); });
}

/**
 * Expects an object where KEY is an arbitrary, user-defined name and
 * the VALUE is another object defining rhe region type and other
 * required properties. Example:
 *
 * ~~~ yaml
 * mysubspace:
 *     policy: within_molecule_type
 *     molecule: water
 *     threshold: 7
 *     com: true
 * ~~~
 */
std::unique_ptr<RegionBase> createRegion(const Space& spc, const json& j) {
    const auto type = j.at("policy").get<RegionType>();
    switch (type) {
    case RegionType::WITHIN_MOLID:
        return std::make_unique<WithinMoleculeType>(spc, j);
    case RegionType::WITHIN_PARTICLE:
        return std::make_unique<SphereAroundParticle>(spc, j);
    case RegionType::WITHIN_ELLIPSOID:
        return std::make_unique<MovingEllipsoid>(spc, j);
    default:
        throw ConfigurationError("unknown region type");
    }
}
void to_json(json& j, const RegionBase& region) {
    region.to_json(j);
    j["policy"] = region.type;
}

void WithinMoleculeType::to_json(json& j) const {
    j = {{"threshold", sqrt(threshold_squared)}, {"com", use_mass_center}, {"molecule", molecules.at(molid).name}};
}

/**
 * @param spc Space to operate on
 * @param molecule_name Molecule type to target
 * @param threshold Radial distance threshold to particles or mass center in target groups
 * @param use_mass_center Use mass-center distance for molecular groups
 * @throw if `use_mass_center==true` and the given molecule type is not molecular, i.e. has ill-defined mass center
 */
WithinMoleculeType::WithinMoleculeType(const Space& spc, std::string_view molecule_name, double threshold,
                                       bool use_mass_center)
    : RegionBase(RegionType::WITHIN_MOLID), spc(spc), molid(findMoleculeByName(molecule_name).id()),
      use_mass_center(use_mass_center), threshold_squared(threshold * threshold) {
    if (use_mass_center && !Faunus::molecules.at(molid).isMolecular()) {
        throw ConfigurationError("center of mass threshold ill-defined for `{}`", molecule_name);
    }
}

WithinMoleculeType::WithinMoleculeType(const Space& spc, const json& j)
    : WithinMoleculeType(spc, j.at("molname").get<std::string>(), j.at("threshold").get<double>(),
                         j.value("com", false)) {}

bool WithinMoleculeType::isInside(const Point& position) const {
    using ranges::cpp20::any_of;
    auto is_inside_group = [&](const Group& group) {
        if (use_mass_center) {
            return within_threshold(position, *group.massCenter());
        }
        return any_of(group.positions(),
                      [&](const auto& position_in_group) { return within_threshold(position, position_in_group); });
    };
    auto groups = spc.findMolecules(molid, Space::Selection::ACTIVE);
    return any_of(groups, is_inside_group);
}

std::optional<double> WithinMoleculeType::volume() const {
    if (use_mass_center) {
        return 4.0 * pc::pi / 3.0 * std::pow(threshold_squared, 1.5);
    }
    return std::nullopt;
}

bool WithinMoleculeType::within_threshold(const Point& position1, const Point& position2) const {
    return spc.geometry.sqdist(position1, position2) < threshold_squared;
}

/**
 * @param spc Space to operate on
 * @param index Particle index used to define center of region
 * @param radius Radius of spherical region
 */
SphereAroundParticle::SphereAroundParticle(const Space& spc, ParticleVector::size_type index, double radius)
    : RegionBase(RegionType::WITHIN_PARTICLE), spc(spc), particle_index(index), radius_squared(radius * radius) {}

SphereAroundParticle::SphereAroundParticle(const Space& spc, const json& j)
    : SphereAroundParticle(spc, j.at("index").get<double>(), j.at("radius").get<double>()) {}

bool SphereAroundParticle::isInside(const Point& position) const {
    return spc.geometry.sqdist(position, spc.particles.at(particle_index).pos) < radius_squared;
}
std::optional<double> SphereAroundParticle::volume() const {
    return 4.0 * pc::pi / 3.0 * std::pow(radius_squared, 1.5);
}
void SphereAroundParticle::to_json(json& j) const {
    j = {{"index", particle_index}, {"threshold", std::sqrt(radius_squared)}};
}

bool MovingEllipsoid::isInside(const Point& position) const {
    Point r21_half = 0.5 * spc.geometry.vdist(reference_position_2, reference_position_1); // half ref2 -> ref1
    const auto r21_half_len = r21_half.norm();
    if (parallel_radius_squared < r21_half_len) {
        faunus_logger->error("Parallel radius ({} Å) smaller than half distance between reference atoms ({} Å)",
                             parallel_radius, r21_half_len);
    }
    auto midpoint = spc.geometry.vdist(reference_position_2, r21_half); // between ref1 & ref2
    Point midpoint_pos = spc.geometry.vdist(position, midpoint);        // midpoint -> pos
    const auto midpoint_pos_len = midpoint_pos.norm();
    const auto cos_theta = midpoint_pos.dot(r21_half) / midpoint_pos_len / r21_half_len;
    const auto theta = std::acos(cos_theta);
    const auto x = cos_theta * midpoint_pos_len;
    const auto y = std::sin(theta) * midpoint_pos_len;
    const auto coord =
        x * x / parallel_radius_squared + y * y / (perpendicular_radius * perpendicular_radius); // normalized coord
    return coord < 1.0; // (>1.0 → outside, <1.0 → inside)
}

/**
 * @param spc Space to operate on
 * @param particle_index1 Index of first particle defining ellipsoid
 * @param particle_index2 Index of second particle defining ellipsoid
 * @param parallel_radius Ellipsoidal radius along axis connecting reference atoms
 * @param perpendicular_radius Ellipsoidal radius perpendicular to axis connecting reference atoms
 */
MovingEllipsoid::MovingEllipsoid(const Space& spc, ParticleVector::size_type particle_index1,
                                 ParticleVector::size_type particle_index2, double parallel_radius,
                                 double perpendicular_radius)
    : RegionBase(RegionType::WITHIN_ELLIPSOID), spc(spc), particle_index_1(particle_index1),
      particle_index_2(particle_index2), parallel_radius(parallel_radius), perpendicular_radius(perpendicular_radius),
      parallel_radius_squared(parallel_radius * parallel_radius),
      reference_position_1(spc.particles.at(particle_index_1).pos),
      reference_position_2(spc.particles.at(particle_index_2).pos) {}

MovingEllipsoid::MovingEllipsoid(const Space& spc, const json& j)
    : MovingEllipsoid(spc, j.at("index1").get<int>(), j.at("index2").get<int>(), j.at("parallel_radius").get<double>(),
                      j.at("perpendicular_radius").get<double>()) {}

void MovingEllipsoid::to_json(json& j) const {
    j["index1"] = particle_index_1;
    j["index2"] = particle_index_1;
    j["perpendicular_radius"] = perpendicular_radius;
    j["parallel_radius"] = parallel_radius;
}

TEST_CASE("[Faunus] Region::MovingEllipsoid") {
    Space spc;
    spc.geometry = Geometry::Chameleon(Geometry::Sphere(100), Geometry::Variant::SPHERE);
    spc.particles.resize(2);
    spc.particles.at(0).pos = {-4.0, 0.0, 0.0}; // first reference
    spc.particles.at(1).pos = {4.0, 0.0, 0.0};  // second reference
    MovingEllipsoid region(spc, 0, 1, 5, 5);

    // add checks for a number of point to verify if inside or outside
    CHECK(region.isInside({0.0, 0.0, 0.0}));
}

} // namespace Region
} // namespace Faunus
