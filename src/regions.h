#pragma once
#include "core.h"

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
namespace Region {

/**
 * @brief Base class for defining sub-spaces of a simulation
 *
 * A region is a sub-space of the system, for example a sphere, cuboid,
 * or other more complex shapes. The class can determine if a point
 * is inside or outside the region. A region need not be static and
 * follow molecules, for example (see `WithinGroups`).
 */
class RegionBase {
  public:
    enum RegionType { SPHERE = 0, CUBOID, WITHIN, NONE };
    const std::map<RegionType, std::string> map = {{SPHERE, "sphere"}, {WITHIN, "within"}};
    RegionType type = NONE;
    std::string name;                               //!< User defined name, may be freely changed
    virtual bool isInside(const Point &) const = 0; //!< true if point is inside region
    virtual double volume() const = 0;              //!< volume of region (-1 if ill defined)
    virtual void to_json(json &) const = 0;
    RegionBase(RegionType);
    virtual ~RegionBase() = default;
};

/**
 * @brief Serialize region to json
 */
void to_json(json &, const std::shared_ptr<RegionBase> &);

/*
 * @brief Factory function to generate all known regions from json
 */
std::shared_ptr<RegionBase> createRegion(const json &, Space &spc);

/**
 * @brief Within one or many groups given by indexes or by molecule names
 *
 * If all groups are molecular, `com` can be used to check for a spherical
 * volume around the mass center. If so, `volume()` returns the
 * spherical volume, otherwise -1.
 */
class WithinGroups : public RegionBase {
  private:
    Space &spc;                  //!< reference to space
    std::vector<size_t> indexes; //!< group indexes to check
    bool com = false;            //!< true = with respect to center of mass
    double threshold2;           //!< squared distance threshold from other particles or com

  public:
    WithinGroups(const json &, Space &);
    static const std::string type_name; //!< Static, fixed name of type
    bool isInside(const Point &) const override;
    double volume() const override;
    void to_json(json &) const override;
};

} // namespace Region
} // namespace Faunus
