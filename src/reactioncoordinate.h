#pragma once

#include "core.h"
#include "group.h"
#include "space.h"
#include <Eigen/Dense>

namespace Faunus {

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
    SystemProperty(const json &j, Tspace &spc);
    void _to_json(json &j) const override;
};

class AtomProperty : public ReactionCoordinateBase {
  protected:
    size_t index; // atom index
    Point dir = {0, 0, 0};

  public:
    std::string property;
    AtomProperty() = default;
    AtomProperty(const json &j, Tspace &spc);
    void _to_json(json &j) const override;
};

struct MoleculeProperty : public AtomProperty {
  protected:
    std::vector<size_t> indexes;

  public:
    MoleculeProperty(const json &j, Tspace &spc);
    void _to_json(json &j) const override;
};

/**
 * @brief Reaction coordinate: molecule-molecule mass-center separation
 */
struct MassCenterSeparation : public ReactionCoordinateBase {
    Eigen::Vector3i dir = {1, 1, 1};
    std::vector<size_t> indexes;
    std::vector<std::string> type;
    MassCenterSeparation(const json &j, Tspace &spc);
    double normalize(double coord) const override; // normalize by volume element
    void _to_json(json &j) const override;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] MassCenterSeparation") {
    using doctest::Approx;
    Tspace spc;
    MassCenterSeparation c(R"({"dir":[1,1,0], "indexes":[0,8,9,18], "type":[] })"_json, spc);
    CHECK(c.dir.x() == 1);
    CHECK(c.dir.y() == 1);
    CHECK(c.dir.z() == 0);
    CHECK(c.indexes == decltype(c.indexes)({0, 8, 9, 18}));
}
#endif

} // namespace ReactionCoordinate
} // namespace Faunus
