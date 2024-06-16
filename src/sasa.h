#pragma once

#ifndef FAUNUS_SASA_H
#define FAUNUS_SASA_H

#include "celllistimpl.h"
#include "particle.h"
#include <range/v3/numeric.hpp>
#include <numbers>

namespace Faunus {

class Space;

namespace Geometry {
class Chameleon;
}

struct Change;

namespace SASA {

using index_type = size_t;
using GeometryType = Geometry::Chameleon;

/**
 * @brief base class for calculating solvent accessible surface areas of target particles
 *        derived classes implement specific neighbour search algorithms
 *
 */
class SASABase {
  public:
    struct Neighbours {
        std::vector<index_type> indices; //!< indices of neighbouring particles in ParticleVector
        PointVector points;              //!< vectors to neighbouring particles
        index_type index;                //!< index of particle whose neighbours are in indices
    };

    bool needs_syncing = false;
    //!< flag indicating if syncing of cell_lists is needed
    //!< this is important  in case there is a particle insertion
    //!< which is rejected due to containerOverlap, because the sasa->update is not called

  protected:
    double probe_radius = 1.4;      //!< radius of the probe sphere
    std::vector<double> areas;      //!< vector holding SASA area of each atom
    std::vector<double> sasa_radii; //!< Radii buffer for all particles
    int slices_per_atom = 20;       //!< number of slices of each sphere in SASA calculation
    const double two_pi = 2.0 * std::numbers::pi;
    const Particle* first_particle; //! first particle in ParticleVector

    /**
     * @brief returns absolute index of particle in ParticleVector
     * @param particle
     */
    [[nodiscard]] inline index_type indexOf(const Particle& particle) const
    {
        return static_cast<index_type>(std::addressof(particle) - first_particle);
    }

    [[nodiscard]] double calcSASAOfParticle(const Neighbours& neighbour) const;
    double exposedArcLength(std::vector<std::pair<double, double>>& arcs) const;

  public:
    [[nodiscard]] double calcSASAOfParticle(const Space& spc, const Particle& particle) const;

    /**
     * @brief calculates total sasa of either particles or groups between given iterators
     * @param spc
     * @param begin iterator to either a first particle or group in a range
     * @param end  iterator to either end of particle or group in a range
     * @tparam TBegin
     * @tparam TEnd
     */
    template <typename TBegin, typename TEnd> double calcSASA(const Space& spc, TBegin begin, TEnd end) const
    {
        return ranges::accumulate(begin, end, 0.0, [&spc, this](auto& area, const auto& species) {
            return area + this->calcSASAOf(spc, species);
        });
    }

    template <typename TSpecies> double calcSASAOf(const Space& spc, const TSpecies& species) const;

    void updateSASA(const std::vector<SASABase::Neighbours>& neighbours_data,
                    const std::vector<index_type>& target_indices);

    virtual void init(const Space& spc) = 0;
    [[nodiscard]] virtual std::vector<SASABase::Neighbours>
    calcNeighbourData(const Space& spc, const std::vector<index_type>& target_indices) const = 0;
    [[nodiscard]] virtual SASABase::Neighbours calcNeighbourDataOfParticle(const Space& spc,
                                                                           index_type target_index) const = 0;
    virtual void update(const Space& spc, const Change& change) = 0;
    [[nodiscard]] const std::vector<double>& getAreas() const;
    SASABase(const Space& spc, double probe_radius, int slices_per_atom);
    virtual ~SASABase() = default;
};

/**
 * @brief derived class of SASABase which uses O(N^2) neighbour search
 *
 **/
class SASA : public SASABase {
  public:
    void init(const Space& spc) override;
    [[nodiscard]] std::vector<SASABase::Neighbours>
    calcNeighbourData(const Space& spc, const std::vector<index_type>& target_indices) const override;
    [[nodiscard]] SASA::Neighbours calcNeighbourDataOfParticle(const Space& spc,
                                                               index_type target_index) const override;
    void update([[maybe_unused]] const Space& spc, [[maybe_unused]] const Change& change) override;
    SASA(const Space& spc, double probe_radius, int slices_per_atom);
    SASA(const json& j, const Space& spc);
};

/**
 * @brief derived class of SASABase which uses cell lists for neighbour search
 *        currently can be used only for geometries with either 3 or 0 periodic boundaries
 *
 * @todo  update function does perhaps unnecessary containsMember(Member&) checks
 * @todo  create a wrapper class for cell_list so that the Space dependence is in there and not here
 **/
template <typename CellList> class SASACellList : public SASABase {
  private:
    using CellCoord = typename CellList::Grid::CellCoord;
    std::unique_ptr<CellList> cell_list; //!< pointer to cell list
    double cell_length;                  //!< dimension of a single cell
    std::vector<CellCoord> cell_offsets; //!< holds offsets which define a 3x3x3 cube around central cell

  public:
    SASACellList(const Space& spc, double probe_radius, int slices_per_atom);
    SASACellList(const json& j, const Space& spc);
    ~SASACellList() override = default;
    void init(const Space& spc) override;
    [[nodiscard]] SASABase::Neighbours calcNeighbourDataOfParticle(const Space& spc,
                                                                   index_type target_index) const override;
    [[nodiscard]] std::vector<SASABase::Neighbours>
    calcNeighbourData(const Space& spc, const std::vector<index_type>& target_indices) const override;
    void update(const Space& spc, const Change& change) override;

  private:
    template <typename TBegin, typename TEnd> void createCellList(TBegin begin, TEnd end, const GeometryType& geometry);
    void updateMatterChange(const Space& spc, const Change& change);
    void updatePositionsChange(const Space& spc, const Change& change);
};

using PeriodicGrid = CellList::Grid::Grid3DPeriodic;
using FixedGrid = CellList::Grid::Grid3DFixed;

template <typename TMember, typename TIndex>
using DenseContainer = CellList::Container::DenseContainer<TMember, TIndex>;
template <typename TMember, typename TIndex>
using SparseContainer = CellList::Container::SparseContainer<TMember, TIndex>;

template <class TGrid, template <typename, typename> class TContainer = DenseContainer>
using CellListType =
    CellList::CellListSpatial<CellList::CellListType<index_type, TGrid, CellList::CellListBase, TContainer>>;

using DensePeriodicCellList = CellListType<PeriodicGrid, DenseContainer>;
using DenseFixedCellList = CellListType<FixedGrid, DenseContainer>;
using SparsePeriodicCellList = CellListType<PeriodicGrid, SparseContainer>;
using SparseFixedCellList = CellListType<FixedGrid, SparseContainer>;

extern template class SASACellList<DensePeriodicCellList>;
extern template class SASACellList<DenseFixedCellList>;
extern template class SASACellList<SparsePeriodicCellList>;
extern template class SASACellList<SparseFixedCellList>;

} // namespace SASA
} // namespace Faunus

#endif // FAUNUS_SASA_H
