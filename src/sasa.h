#pragma once

#ifndef FAUNUS_SASA_H
#define FAUNUS_SASA_H

#include "celllistimpl.h"
#include "particle.h"

namespace Faunus {

class Space;
namespace Geometry {
class Chameleon;
}
class Change;

namespace SASA {

using index_type = size_t;
using GeometryType = Geometry::Chameleon;

class SASABase {

  public:

    struct Neighbours {
        std::vector<index_type> indices; //!< indices of neighbouring particles in ParticleVector
        PointVector points;          //!< vectors to neighbouring particles
        index_type index;            //!< index of particle whose neighbours are in indices
    };

  protected:
    double probe_radius;            //!< radius of the probe sphere
    std::vector<double> areas;      //!< vector holding SASA area of each atom
    std::vector<double> radii;      //!< Radii buffer for all particles
    int slices_per_atom;            //!< number of slices of each sphere in SASA calculation
    const double TWOPI = 2. * M_PI; //!< 2 PI
    const Particle* first_particle; //! first particle in ParticleVector

    /**
     * @brief returns absolute index of particle in ParticleVector
     * @param particle
     */
    inline index_type indexOf(const Particle& particle) {
        return static_cast<index_type>(std::addressof(particle) - first_particle);
    }

    double calcSASAOfParticle(const Neighbours& neighbour) const;

    double exposedArcLength(std::vector<std::pair<double, double>>& arcs) const;

  public:
    void updateSASA(const std::vector<SASABase::Neighbours>& neighbours_data,
                    const std::vector<index_type>& target_indices);

    virtual void init(Space& spc) = 0;

    virtual std::vector<SASABase::Neighbours> calcNeighbourData(Space& spc,
                                                                const std::vector<index_type>& target_indices) = 0;

    virtual SASABase::Neighbours calcNeighbourDataOfParticle(Space& spc, const index_type target_index) = 0;

    virtual void update(Space& spc, const Change& change) = 0;

    const std::vector<double>& getAreas() const;

    SASABase(Space& spc, double probe_radius, int slices_per_atom);
    virtual ~SASABase() = default;
};

class SASA : public SASABase {

  public:
    void init(Space& spc) override;

    std::vector<SASABase::Neighbours> calcNeighbourData(Space& spc,
                                                        const std::vector<index_type>& target_indices) override;

    SASA::Neighbours calcNeighbourDataOfParticle(Space& spc, const index_type target_index);

    void update([[maybe_unused]] Space& spc, [[maybe_unused]] const Change& change) override;

    SASA(Space& spc, double probe_radius, int slices_per_atom);
    SASA(const json& j, Space& spc);
};

//!< TODO update function does perhaps unnecessary containsMember(Member&) checks
//!< TODO create a wrapper class for cell_list so that the Space dependence is in there and not here
template <typename CellList> class SASACellList : public SASABase {

  private:
    using CellCoord = typename CellList::Grid::CellCoord;

    std::unique_ptr<CellList> cell_list;   //!< pointer to cell list
    double cell_length;                    //!< dimension of a single cell
    std::vector<CellCoord> cell_offsets;   //!< holds offsets which define a 3x3x3 cube around central cell

  public:
    SASACellList(Space& spc, double probe_radius, int slices_per_atom);
    SASACellList(const json& j, Space& spc);
    virtual ~SASACellList() = default;

    void init(Space& spc) override;

    SASABase::Neighbours calcNeighbourDataOfParticle(Space& spc, const index_type target_index) override;

    std::vector<SASABase::Neighbours> calcNeighbourData(Space& spc,
                                                        const std::vector<index_type>& target_indices) override;

    void update(Space& spc, const Change& change) override;

  private:
    template <typename TBegin, typename TEnd> void createCellList(TBegin begin, TEnd end, GeometryType& geometry);

    void updateMatterChange(Space& spc, const Change& change);

    void updatePositionsChange(Space& spc, const Change& change);
};

using PeriodicCellList = CellList::CellListSpatial<CellList::CellListType<index_type, CellList::Grid::Grid3DPeriodic>>;
using FixedCellList = CellList::CellListSpatial<CellList::CellListType<index_type, CellList::Grid::Grid3DFixed>>;
extern template class SASACellList<PeriodicCellList>;
extern template class SASACellList<FixedCellList>;

} // namespace SASA

} // namespace Faunus

#endif // FAUNUS_SASA_H
