#pragma once

#include "core.h"
#include "group.h"

namespace Faunus {

struct Space;

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
struct ReactionCoordinateBase {
    ReactionCoordinateBase(const json &);   //!< constructor reads binwidth, min, max
    std::function<double()> f = nullptr;    //!< returns reaction coordinate
    virtual void _to_json(json &j) const;   //!< json serialization
    virtual double normalize(double) const; //!< Default 1.0; currently unused
    double binwidth = 0, min = 0, max = 0;
    std::string name; //!< Meaningful, short name. Don't use spaces or weird characters

    double operator()(); //!< Calculates reaction coordinate

    bool inRange(double coord) const; //!< Determines if coordinate is within [min,max]
    virtual ~ReactionCoordinateBase() = default;
};

void to_json(json &j, const ReactionCoordinateBase &r); //!< Serialize any reaction coordinate to json

std::shared_ptr<ReactionCoordinateBase>
createReactionCoordinate(const json &, Space &); //!< Factory function to create all known penalty functions

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] ReactionCoordinateBase") {
    using doctest::Approx;
    ReactionCoordinateBase c(R"({"range":[-1.5, 2.1], "resolution":0.2})"_json);
    CHECK(c.min == Approx(-1.5));
    CHECK(c.max == Approx(2.1));
    CHECK(c.binwidth == Approx(0.2));
    CHECK(c.inRange(-1.5) == true);
    CHECK(c.inRange(-1.51) == false);
    CHECK(c.inRange(2.11) == false);
    CHECK(c.inRange(2.1) == true);
}
#endif

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
