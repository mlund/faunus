#pragma once

/**
regions:
    subspace1:
        type: sphere
        radius: 20
        origin: [0,0,5]
    subspace2:
        type: within
        radius: 7
        molecules: [Na+, Cl-]
        com: true
*/

namespace Faunus {
namespace Regions {

class RegionBase {
  private:
    virtual void _from_json(const json &) = 0;

  public:
    virtual bool isInside(const Point &) const = 0;
    virtual double volume() const = 0;
    void from_json(const json &);
};

class WithinMolecule : public RegionBase {
  private:
    double threshold;
    Point origo;
    void _from_json(const json &) override;

  public:
    bool isInside(const Point &) const override;
    double volume() const override;
};

} // namespace Regions
} // namespace Faunus
