#include "regions.h"
#include "space.h"
#include <range/v3/algorithm/any_of.hpp>
#include <iostream>

namespace Faunus {
namespace Region {

RegionBase::RegionBase(RegionType type) : type(type) {}

std::optional<double> RegionBase::volume() const { return std::nullopt; }

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
        return std::make_unique<VidarsRegion>(spc, j);
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

bool VidarsRegion::isInside(const Point& position) const {
    const auto ref1_pos = spc.particles.at(particle_index1).pos;
    const auto ref2_pos = spc.particles.at(particle_index2).pos;

    Point cylAxis = spc.geometry.vdist(ref2_pos, ref1_pos) * 0.5; // half vector between reference atoms
    if (parallel_radius < cylAxis.norm()) {                       // checking so that a is larger than length of cylAxis
        throw std::runtime_error(
            "specified radius of ellipsoid along the axis connecting reference atoms (rx) must be larger or equal "
            "to half the distance between reference atoms. Specified radius is " +
            std::to_string(parallel_radius) + " Å whereas half the distance between reference atoms is " +
            std::to_string(cylAxis.norm()) + "Å");
    }
    Point origin = ref2_pos - cylAxis;                // coordinates of middle point between reference atoms: new origo
    auto molV = spc.geometry.vdist(position, origin); // vector between selected molecule and center of geometry
    auto cosTheta = molV.dot(cylAxis) / molV.norm() / cylAxis.norm(); // cosinus of angle between coordinate
    // vector of selected molecule and axis
    // connecting reference atoms
    auto theta = acos(cosTheta);     // angle between coordinate vector of sel. molecule and axis connecting ref. atoms
    auto x = cosTheta * molV.norm(); // x coordinate of selected molecule with respect to center of geometry
    // (in plane including vectors molV and cylAxis)
    auto y = sin(theta) * molV.norm(); // y coordinate of selected molecule with respect to center of geometry
    // (in plane including vectors molV and cylAxis)
    auto coord =
        x * x / (parallel_radius * parallel_radius) +
        y * y / (perpendicular_radius * perpendicular_radius); // calculating normalized coordinate with respect to
    // dimensions of geometry (>1.0 → outside, <1.0 → inside)
    return coord <= 1.0;
}

VidarsRegion::VidarsRegion(const Space& spc, ParticleVector::size_type particle_index1,
                           ParticleVector::size_type particle_index2, double r_x, double r_y)
    : RegionBase(RegionType::WITHIN_ELLIPSOID), spc(spc), particle_index1(particle_index1),
      particle_index2(particle_index2), parallel_radius(r_x), perpendicular_radius(r_y) {}

VidarsRegion::VidarsRegion(const Space& spc, const json& j)
    : VidarsRegion(spc, j.at("index1").get<int>(), j.at("index2").get<int>(), j.at("r_x").get<double>(),
                   j.at("r_y").get<double>()) {}

void VidarsRegion::to_json(json& j) const {
    j["index1"] = particle_index1;
    j["index2"] = particle_index1;
    j["perpendicular_radius"] = perpendicular_radius;
    j["parallel_radius"] = parallel_radius;
}

} // namespace Region
} // namespace Faunus
