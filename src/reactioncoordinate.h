#pragma once

#include "core.h"
#include "group.h"

namespace Faunus {

struct Space;

namespace ReactionCoordinate {

/**
 * @brief Base class for reaction coordinates
 */
struct ReactionCoordinateBase {
    std::function<double()> f = nullptr; // returns reaction coordinate
    virtual void _to_json(json &j) const;
    virtual double normalize(double) const;
    double binwidth = 0, min = 0, max = 0;
    std::string name;

    double operator()(); //!< Calculates reaction coordinate

    bool inRange(double coord) const; //!< Determines if coordinate is within [min,max]
    virtual ~ReactionCoordinateBase() = default;
};

void to_json(json &j, const ReactionCoordinateBase &r);
void from_json(const json &j, ReactionCoordinateBase &r);

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] ReactionCoordinateBase") {
    using doctest::Approx;
    ReactionCoordinateBase c = R"({"range":[-1.5, 2.1], "resolution":0.2})"_json;
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
    AtomProperty() = default;
    AtomProperty(const json &j, Space &spc);
    void _to_json(json &j) const override;
};

struct MoleculeProperty : public AtomProperty {
  protected:
    std::vector<size_t> indexes;

  public:
    MoleculeProperty(const json &j, Space &spc);
    void _to_json(json &j) const override;
};

} // namespace ReactionCoordinate
} // namespace Faunus
