#include "regions.h"
#include "space.h"
#include <range/v3/algorithm/any_of.hpp>
#include <iostream>

namespace Faunus {
namespace Region {

RegionBase::RegionBase(RegionType type) : type(type) {}

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

} // namespace Region
} // namespace Faunus
