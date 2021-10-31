#pragma once

#ifndef FAUNUS_SASA_H
#define FAUNUS_SASA_H

#include "atomdata.h"
#include "particle.h"

namespace Faunus {

class Space;
struct Change;

class SASA {
  protected:
    using Index = AtomData::index_type;

  public:
    struct Neighbours {
        std::vector<size_t> indices; //!< indices of neighbouring particles in ParticleVector
        PointVector points;          //!< vectors to neighbouring particles
        Index index;                 //!< index of particle which corresponds to the object
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
    inline size_t indexOf(const Particle& particle) {
        return static_cast<size_t>(std::addressof(particle) - first_particle);
    }

    /**
     * @brief Calcuates SASA of a single particle defined by NeighbourData object
     * @param neighbours NeighbourData object of given particle
     * @param radii of the particles in the ParticleVector
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
     * @param radii of the particles in the ParticleVector
     * @param target_indices absolute indicies of target particles in ParticleVector
     */
    void updateSASA(const std::vector<SASA::Neighbours>& neighbours_data, const std::vector<size_t>& target_indices);

    /**
     * @brief resizes areas buffer to size of ParticleVector
     * @param space
     */
    void init(ParticleVector& particles);

    /**
     * @brief calculates neighbourData objects of particles specified by target indices in ParticleVector
     * @param space
     * @param target_indices absolute indicies of target particles in ParticleVector
     */
    std::vector<SASA::Neighbours> calcNeighbourData(Space& spc, const std::vector<size_t>& target_indices);

    /**
     * @brief calculates neighbourData object of a target particle specified by target indiex in ParticleVector
     * @param space
     * @param target_index indicex of target particle in ParticleVector
     */
    SASA::Neighbours calcNeighbourDataOfParticle(Space& spc, const size_t target_index);

    void update([[maybe_unused]] Space& spc, [[maybe_unused]] const Change& change);

    const std::vector<double>& getAreas() const;

    /**
     * @param spc
     * @param probe_radius in angstrom
     * @param slices_per_atom number of slices of each sphere in sasa calculations
     */
    SASA(Space& spc, double probe_radius, int slices_per_atom);

    SASA(const json& j, Space& spc);

};

} // namespace Faunus

#endif // FAUNUS_SASA_H