#pragma once

#ifndef FAUNUS_SASA_H
#define FAUNUS_SASA_H

#include "celllistimpl.h"
#include "particle.h"

namespace Faunus {

class Space;
class Change;

namespace SASA {

class SASABase {

  public:
    using AtomIndex = AtomData::index_type;

    struct Neighbours {
        std::vector<size_t> indices; //!< indices of neighbouring particles in ParticleVector
        PointVector points;          //!< vectors to neighbouring particles
        AtomIndex index;             //!< index of particle whose neighbours are in indices
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
    inline AtomIndex indexOf(const Particle& particle) {
        return static_cast<AtomIndex>(std::addressof(particle) - first_particle);
    }

    /**
     * @brief Calcuates SASA of a single particle defined by NeighbourData object
     * @param neighbours NeighbourData object of given particle
     */
    double calcSASAOfParticle(const Neighbours& neighbour) const;

    /**
     * @brief Calcuates total arc length in radians of overlapping arcs defined by two angles
     * @param vector of arcs, defined by a pair (first angle, second angle)
     */
    double exposedArcLength(std::vector<std::pair<double, double>>& arcs) const;

  public:
    /**
     * @brief updates sasa of target particles
     * @param neighbours_data array of NeighbourData objects
     * @param target_indices absolute indicies of target particles in ParticleVector
     */
    void updateSASA(const std::vector<SASABase::Neighbours>& neighbours_data,
                    const std::vector<size_t>& target_indices);

    virtual void init(Space& spc) = 0;

    virtual std::vector<SASABase::Neighbours> calcNeighbourData(Space& spc,
                                                                const std::vector<size_t>& target_indices) = 0;

    virtual SASABase::Neighbours calcNeighbourDataOfParticle(Space& spc, const size_t target_index) = 0;

    virtual void update(Space& spc, const Change& change) = 0;

    const std::vector<double>& getAreas() const;

    /**
     * @param spc
     * @param probe_radius in angstrom
     * @param slices_per_atom number of slices of each sphere in sasa calculations
     */
    SASABase(Space& spc, double probe_radius, int slices_per_atom);
    virtual ~SASABase() {}
};

class SASA : public SASABase {

  public:
    /**
     * @brief resizes areas buffer to size of ParticleVector and fill radii buffer with radii
     * @param space
     */
    void init(Space& spc) override;

    /**
     * @brief calculates neighbourData objects of particles specified by target indices in ParticleVector
     * @param space
     * @param target_indices absolute indicies of target particles in ParticleVector
     */
    std::vector<SASABase::Neighbours> calcNeighbourData(Space& spc, const std::vector<size_t>& target_indices) override;

    /**
     * @brief calculates neighbourData object of a target particle specified by target indiex in ParticleVector
     * @brief using the naive O(N) neighbour search
     * @param space
     * @param target_index indicex of target particle in ParticleVector
     */
    SASA::Neighbours calcNeighbourDataOfParticle(Space& spc, const AtomIndex target_index);

    /**
     * @brief updates radii vector in case of matter change
     * @param space
     * @param change
     */
    void update([[maybe_unused]] Space& spc, [[maybe_unused]] const Change& change) override;

    /**
     * @param spc
     * @param probe_radius in angstrom
     * @param slices_per_atom number of slices of each sphere in sasa calculations
     */
    SASA(Space& spc, double probe_radius, int slices_per_atom);
    SASA(const json& j, Space& spc);
};

//!< TODO update function does perhaps unnecessary containsMember(Member&) checks
//!< TODO finish proper test case
//!< TODO create a wrapper class for cell_list so that the Space dependence is in there and not here
template <typename CellList_T> class SASACellList : public SASABase {

  private:
    using CellCoord = typename CellList_T::Grid::CellCoord;

    std::unique_ptr<CellList_T> cell_list; //!< pointer to cell list
    double cell_length;                    //!< dimension of a single cell
    std::vector<CellCoord> cell_offsets;   //!< holds offsets which define a 3x3x3 cube around central cell

  public:
    /**
     * @param spc
     * @param probe_radius in angstrom
     * @param slices_per_atom number of slices of each sphere in sasa calculations
     */
    SASACellList(Space& spc, double probe_radius, int slices_per_atom);
    SASACellList(const json& j, Space& spc);
    virtual ~SASACellList() {}

    /**
     * @brief constructs cell_list with appropriate cell_length and fills it with particles from space
     * @param space
     */
    void init(Space& spc) override;
    /**
     * @brief calculates neighbourData object of a target particle
     * @brief specified by target index in ParticleVector using cell list
     * @param space
     * @param target_index indicex of target particle in ParticleVector
     */
    SASABase::Neighbours calcNeighbourDataOfParticle(Space& spc, const size_t target_index) override;

    /**
     * @brief calculates neighbourData objects of particles
     * @brief specified by target indices in ParticleVector using cell lists
     * @param space
     * @param target_indices absolute indicies of target particles in ParticleVector
     */
    std::vector<SASABase::Neighbours> calcNeighbourData(Space& spc, const std::vector<size_t>& target_indices) override;

    /**
     * @brief updates cell_list according to change, if the volume changes the cell_list gets rebuilt
     * @brief also updates radii in case of matter change
     * @param space
     * @param change
     */
    void update(Space& spc, const Change& change) override;

  private:
    /**
     * @brief creates cell_list if it does not exist, otherwise resets its and fills it with active particles
     * @param space
     */
    void createCellList(Space& spc);

    /**
     * @brief updates cell_list when particles got activated or desactivated
     * @param space
     * @param change
     */
    void updateMatterChange(Space& spc, const Change& change);
};

using PeriodicCellList = CellList::CellListSpatial<CellList::CellListType<size_t, CellList::Grid::Grid3DPeriodic>>;
using FixedCellList = CellList::CellListSpatial<CellList::CellListType<size_t, CellList::Grid::Grid3DFixed>>;
extern template class SASACellList<PeriodicCellList>;
extern template class SASACellList<FixedCellList>;

} // namespace SASA

} // namespace Faunus

#endif // FAUNUS_SASA_H