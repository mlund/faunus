#pragma once

#include "core.h"
#include "group.h"

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
class ReactionCoordinateBase {
  protected:
    std::function<double()> function = nullptr;    //!< returns reaction coordinate
    virtual double normalize(double) const; //!< Default 1.0; currently unused
  public:
    ReactionCoordinateBase(const json &);   //!< constructor reads resolution, min, max
    double resolution = 0.0;                //!< Resolution used when binning (histograms etc.)
    double minimum_value = 0.0;             //!< Minimum allowed value
    double maximum_value = 0.0;             //!< Maximum allowed value
    std::string name;                       //!< Meaningful, short name. Don't use spaces or weird characters

    double operator()(); //!< Calculates reaction coordinate
    virtual void _to_json(json &j) const;   //!< json serialization
    bool inRange(double coord) const; //!< Determines if coordinate is within [min,max]
    virtual ~ReactionCoordinateBase() = default;
};

void to_json(json &j, const ReactionCoordinateBase &r); //!< Serialize any reaction coordinate to json

std::unique_ptr<ReactionCoordinateBase>
createReactionCoordinate(const json&, Space&); //!< Factory function to create all known penalty functions

class SystemProperty : public ReactionCoordinateBase {
  protected:
    std::string property;

  public:
    SystemProperty(const json &j, Space &spc);
    void _to_json(json &j) const override;
};

class AtomProperty : public ReactionCoordinateBase {
  protected:
    size_t index; // atom index
    Point dir = {0, 0, 0};

  public:
    std::string property;
    AtomProperty(const json &j, Space &spc);
    void _to_json(json &j) const override;
};

struct MoleculeProperty : public ReactionCoordinateBase {
  protected:
    size_t index; // molecule index
    Point dir = {0, 0, 0};
    std::vector<size_t> indexes;

  public:
    std::string property;
    MoleculeProperty(const json &j, Space &spc);
    void _to_json(json &j) const override;
};

} // namespace ReactionCoordinate
} // namespace Faunus
