#pragma once
#include "core.h"
#include "molecule.h"
#include "group.h"
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>

/**
Possible layout:

regions:
    subspace1:
        type: within
        radius: 7
        molecules: [membrane]
        com: true
*/

namespace Faunus {
class Space; // forward declare Space

/**
 * Subspaces around particles, molecules etc.
 */
namespace Region {

enum class RegionType { WITHIN_MOLID, WITHIN_PARTICLE, WITHIN_ELLIPSOID, INVALID };

NLOHMANN_JSON_SERIALIZE_ENUM(RegionType, {{RegionType::INVALID, nullptr},
                                          {RegionType::WITHIN_MOLID, "around_molecule_type"},
                                          {RegionType::WITHIN_PARTICLE, "around_particle"},
                                          {RegionType::WITHIN_ELLIPSOID, "ellipsoid"}})

/**
 * @brief Base class for defining sub-spaces of a simulation
 *
 * A region is a sub-space of the system, for example a sphere, cuboid,
 * or other more complex shapes. The class can determine if a point
 * is inside the region. A region need not be static and can
 * follow molecules. A region may or may not have a well-defined volume.
 */
class RegionBase {
  private:
    virtual bool isInside(const Point& position) const = 0; //!< true if point is inside region
  public:
    const RegionType type;
    virtual std::optional<double> volume() const; //!< Volume of region if applicable
    virtual void to_json(json& j) const = 0;
    virtual ~RegionBase() = default;
    explicit RegionBase(RegionType type);

    bool use_group_mass_center = false; //!< Use group mass-center to check if inside region

    bool inside(const Particle& particle) const; //!< Determines if particle is inside region
    bool inside(const Group& group) const;       //!< Determines of groups is inside region

    /** Selects particles within the region */
    template <typename ParticleRange> auto filterInside(const ParticleRange& particles) const {
        namespace rv = ranges::cpp20::views;
        return particles | rv::transform(&Particle::pos) | rv::filter(&RegionBase::isInside);
    }
};

/*
 * @brief Factory function to generate all known regions from json
 */
std::unique_ptr<RegionBase> createRegion(const Space& spc, const json& j);

void to_json(json& j, const RegionBase& region);

/**
 * @brief Within a spherical cutoff distance from a molecule type
 *
 * Checks if within a radial distance from any particle in any groups
 * of the given molecular id (`molecule_name`).  If the group is molecular,
 * `com` can be used to check for a spherical
 * volume around the mass center. If so, `volume()` returns the
 * spherical volume, otherwise `nullopt`
 */
class WithinMoleculeType : public RegionBase {
  private:
    const Space& spc;                     //!< reference to space
    const MoleculeData::index_type molid; //!< molid to target
    const bool use_mass_center = false;   //!< true = with respect to center of mass
    const double threshold_squared;       //!< squared distance threshold from other particles or com
    bool within_threshold(const Point& position1, const Point& position2) const;

  public:
    WithinMoleculeType(const Space& spc, std::string_view molecule_name, double threshold, bool use_mass_center);
    WithinMoleculeType(const Space& spc, const json& j);
    bool isInside(const Point& position) const override;
    std::optional<double> volume() const override;
    void to_json(json& j) const override;
};

/**
 * Spherical region centered on a particle
 */
class SphereAroundParticle : public RegionBase {
  private:
    const Space& spc;                               //!< reference to space
    const ParticleVector::size_type particle_index; //!< Index of particle that defines the center of the region
    const double radius_squared;                    //!< squared distance threshold from other particles or com

  public:
    SphereAroundParticle(const Space& spc, ParticleVector::size_type index, double radius);
    SphereAroundParticle(const Space& spc, const json& j);
    bool isInside(const Point& position) const override;
    std::optional<double> volume() const override;
    void to_json(json& j) const override;
};

/**
 * An ellipsoid defined by two (moving) particles
 */
class MovingEllipsoid : public RegionBase {
  private:
    const Space& spc;
    const ParticleVector::size_type particle_index_1; //!< Index of first reference particle
    const ParticleVector::size_type particle_index_2; //!< Index of second reference particle
    const double parallel_radius;                     //!< ellipsoidal radius along axis connecting reference atoms
    const double perpendicular_radius; //!< ellipsoidal radius perpendicular to axis connecting reference atoms
    const double parallel_radius_squared;

    const Point& reference_position_1; //!< Reference to first reference particle position
    const Point& reference_position_2; //!< Reference to second reference particle position

  public:
    MovingEllipsoid(const Space& spc, ParticleVector::size_type particle_index1,
                    ParticleVector::size_type particle_index2, double parallel_radius, double perpendicular_radius);
    MovingEllipsoid(const Space& spc, const json& j);
    bool isInside(const Point& position) const override;
    void to_json(json& j) const override;
};

} // namespace Region
} // namespace Faunus
