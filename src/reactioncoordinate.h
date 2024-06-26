#pragma once

#include "core.h"
#include <functional>

namespace Faunus {

class Space;

namespace ReactionCoordinate {

/**
 * @brief Base class for reaction coordinates
 *
 * A reaction coordinate (RC) is an arbitrary, one-dimensional property of
 * the simulated system. This could be the x-position of a specific particle,
 * the system volume, a mass center of something else. RC's are used in
 * the penalty's energy function and can also be used to probe the system
 * during analysis.
 */
class ReactionCoordinateBase
{
  protected:
    std::function<double()> function = nullptr; //!< returns reaction coordinate
                                                //!< Default 1.0; currently unused
  public:
    explicit ReactionCoordinateBase(const json& j); //!< constructor reads resolution, min, max
    double resolution = 0.0;    //!< Resolution used when binning (histograms etc.)
    double minimum_value = 0.0; //!< Minimum allowed value
    double maximum_value = 0.0; //!< Maximum allowed value
    std::string name;           //!< Meaningful, short name. Don't use spaces or weird characters

    double operator()();                  //!< Calculates reaction coordinate
    virtual void _to_json(json& j) const; //!< json serialization
    [[nodiscard]] bool
    inRange(double coord) const; //!< Determines if coordinate is within [min,max]
    virtual ~ReactionCoordinateBase() = default;
};

void to_json(json& j, const ReactionCoordinateBase&
                          reaction_coordinate); //!< Serialize any reaction coordinate to json

std::unique_ptr<ReactionCoordinateBase>
createReactionCoordinate(const json&,
                         const Space&); //!< Factory function to create all known penalty functions

class SystemProperty : public ReactionCoordinateBase
{
  protected:
    std::string property;

  public:
    SystemProperty(const json& j, const Space& spc);
    void _to_json(json& j) const override;
};

class AtomProperty : public ReactionCoordinateBase
{
  protected:
    size_t index; // atom index
    Point dir = {0.0, 0.0, 0.0};

  public:
    std::string property;
    AtomProperty(const json& j, const Space& spc);
    void _to_json(json& j) const override;
};

/**
 * @todo Refactor so that each scheme is a derived class implementing
 *       a virtual energy function instead of the std::function object
 */
struct MoleculeProperty : public ReactionCoordinateBase
{
  private:
    size_t index; //!< Group index
    Point direction = {0.0, 0.0, 0.0};
    std::vector<size_t> indexes;
    std::string property;

    void selectAngleWithVector(const json& j, const Space& spc);
    void selectRinner(const json& j, const Space& spc);
    void selectMinimumGroupDistance(const json& j, const Space& spc);
    void selectMassCenterDistance(const json& j, const Space& spc);
    void selectDipoleAngle(const json& j, const Space& spc, Geometry::BoundaryFunction& b);
    void selectGyrationRadius(const Space& spc);
    void selectAtomAtomDistance(const json& j, const Space& spc);
    void selectMassCenterDistanceZ(const json& j, const Space& spc);
    void selectLengthOverRadiusRatio(const json& j, const Space& spc);

  public:
    MoleculeProperty(const json& j, const Space& spc);
    void _to_json(json& j) const override;
};

} // namespace ReactionCoordinate
} // namespace Faunus
