#pragma once

#ifndef FAUNUS_SASA_H
#define FAUNUS_SASA_H

#include "space.h"
#include "celllistimpl.h"
#include <range/v3/view/iota.hpp>
#include <unordered_set>

namespace Faunus {

class SASA {

  public:
    class NeighboursData {
      public:
        std::vector<int> neighbour_inds; // indices of neighbouring particles in ParticleVector
        PointVector points;              // vectors to neighbouring particles
        int ind;                         // index of particle which corresponds to the object
    };

  protected:
    double probe_radius;
    std::vector<double> areas;
    int n_slices_per_atom;
    const double TWOPI = 2. * M_PI;

  private:
    /**
     * @brief Calcuates SASA of a single particle defined by NeighbourData object
     * @param NeighbourData object of given particle
     * @param radii of ALL the particles
     */
    double calcSASAOfParticle(const NeighboursData& neighbours, const std::vector<double>& radii) const;

    /**
     * @brief Calcuates total arc length in radians of overlapping arcs
     * @param vector of arcs, defined by a pair (beginning angle, end angle)
     */
    double exposedArcLength(std::vector<std::pair<double, double>>& arcs) const;

  public:
    /**
     * @brief updates sasa of target particles
     * @param absolute indicies of target particles in ParticleVector
     * @param neighbourData object of all particles in ParticleVector
     */
    void updateSASA(const std::vector<SASA::NeighboursData>& neighbours, const std::vector<double>& radii,
                    const std::vector<int>& target_inds);

    void init(Space& spc);

    /**
     * @brief calculates neighbourData objects of particles specified by target indices in ParticleVector
     * @param space
     * @param absolute indicies of target particles in ParticleVector
     */
    const std::vector<SASA::NeighboursData> calcNeighbourData(Space& spc, const std::vector<int>& target_inds);

    /**
     * @brief calculates neighbourData object of a target particle specified by target indiex in ParticleVector
     * @param space
     * @param absolute indicies of target particles in ParticleVector
     */
    SASA::NeighboursData calcNeighbourDataOfParticle(Space& spc, const int target_ind);

    void update([[maybe_unused]] Space& spc, [[maybe_unused]] const Change& change) {}

    void sync([[maybe_unused]] Space& other_space, [[maybe_unused]] const Change& change) {}

    const std::vector<double>& getAreas() const { return areas; }

    SASA(const json& j, [[maybe_unused]] Space& spc)
        : SASA(spc, j.value("radius", 1.4) * 1.0_angstrom, j.value("slices", 20)) {}
    SASA([[maybe_unused]] Space& spc, double probe_radius, int n_slices_per_atom)
        : probe_radius(probe_radius), n_slices_per_atom(n_slices_per_atom) {}
};

// not used yet
class SASACellList : public SASA {

  private:
    using CellList_t = CellList::CellListSpatial<CellList::CellListType<int, CellList::Grid::Grid3DPeriodic>>;
    std::unique_ptr<CellList_t> cell_list; // pointer to cell list
    double cell_size;                      // dimension of a single cell
    using CellCoord = CellList::GridOf<CellList_t>::CellCoord;
    std::vector<CellCoord> cell_offsets; // holds offsets which define a 3x3x3 cube around central cell

  public:
    SASACellList(const json& j, Space& spc)
        : SASACellList(spc, j.value("radius", 1.4) * 1.0_angstrom, j.value("slices", 20)) {}

    SASACellList(Space& spc, double probe_radius, int n_slices_per_atom) : SASA(spc, probe_radius, n_slices_per_atom) {}

    void init(Space& spc) {

        cell_offsets.clear();
        for (auto i = -1; i <= 1; ++i) {
            for (auto j = -1; j <= 1; ++j) {
                for (auto k = -1; k <= 1; ++k) {
                    cell_offsets.push_back(CellCoord(i, j, k));
                }
            }
        }

        double max_radius = 0.0;
        for (const auto& particle : spc.particles) {
            if (particle.traits().sigma > max_radius) {
                max_radius = particle.traits().sigma;
            }
        }
        cell_size = max_radius + 2.0 * probe_radius;
        cell_list = std::make_unique<CellList_t>(spc.geometry.getLength(), cell_size);
        for (auto& particle : spc.activeParticles()) {
            const int i = static_cast<int>(std::addressof(particle) - std::addressof(spc.particles.at(0)));
            spc.geometry.boundary(particle.pos);
            cell_list->insertMember(i, particle.pos);
        }

        areas.resize(spc.particles.size());
    }

    SASA::NeighboursData calcNeighbourDataOfParticle(Space& spc, const int target_ind) {

        SASA::NeighboursData neighbour;

        const int i = target_ind;
        const auto& particle_i = spc.particles.at(i);
        neighbour.ind = i;
        const auto rc_i = particle_i.traits().sigma * 0.5 + probe_radius;
        const auto& center_cell = cell_list->getGrid().coordinatesAt(particle_i.pos);

        for (const auto& offset : cell_offsets) {
            const auto& neighbour_particle_inds = cell_list->getNeighborMembers(center_cell, offset);
            for (const int j : neighbour_particle_inds) {

                const auto& particle_j = spc.particles.at(j);
                const double rc_j = particle_j.traits().sigma * 0.5 + probe_radius;
                const auto sq_cutoff = (rc_i + rc_j) * (rc_i + rc_j);
                if (i != j && spc.geometry.sqdist(particle_i.pos, particle_j.pos) <= sq_cutoff) {

                    Point dr = particle_i.pos - particle_j.pos;
                    spc.geometry.boundary(dr);

                    neighbour.neighbour_inds.push_back(j);
                    neighbour.points.push_back(dr);
                }
            }
        }
        return neighbour;
    }

    const std::vector<SASA::NeighboursData> calcNeighbourData(Space& spc, const std::vector<int>& target_inds) {

        std::vector<SASA::NeighboursData> neighbours(spc.particles.size());

        for (const int i : target_inds) {
            const auto& particle_i = spc.particles.at(i);
            neighbours[i].ind = i;
            const auto rc_i = particle_i.traits().sigma * 0.5 + probe_radius;
            const auto& center_cell = cell_list->getGrid().coordinatesAt(particle_i.pos);

            for (const auto& offset : cell_offsets) {
                const auto& neighbour_particle_inds = cell_list->getNeighborMembers(center_cell, offset);
                for (const int j : neighbour_particle_inds) {

                    const auto& particle_j = spc.particles.at(j);
                    const double rc_j = particle_j.traits().sigma * 0.5 + probe_radius;
                    const auto sq_cutoff = (rc_i + rc_j) * (rc_i + rc_j);
                    if (i != j && spc.geometry.sqdist(particle_i.pos, particle_j.pos) <= sq_cutoff) {

                        Point dr = particle_i.pos - particle_j.pos;
                        spc.geometry.boundary(dr);

                        neighbours[i].neighbour_inds.push_back(j);
                        neighbours[i].points.push_back(dr);
                    }
                }
            }
        }
        return neighbours;
    }

    void sync(Space& other_space, const Change& change) { this->update(other_space, change); }

    void update(Space& spc, const Change& change) {

        if (change.everything) {
            cell_list.reset(new CellList_t(spc.geometry.getLength(), cell_size));
            for (auto& particle : spc.activeParticles()) {
                const int i = static_cast<int>(std::addressof(particle) - std::addressof(spc.particles.at(0)));
                spc.geometry.boundary(particle.pos);
                cell_list->insertMember(i, particle.pos);
            }
            return;
        }
        /*else if (!change.matter_change) {
            for (const auto& group_change : change.groups) {
                const auto& group = spc.groups.at(group_change.group_index);
                const auto offset = spc.getFirstParticleIndex(group);
                if( group_change.relative_atom_indices.empty() ){
                    auto indices = ranges::cpp20::views::iota(offset, offset + group.size());
                    for( const auto index : indices){
                        cell_list->updateMemberAt(index, spc.particles.at(index).pos);
                    }
                }
                auto absolute_atom_index = group_change.relative_atom_indices |
                                           ranges::cpp20::views::transform([offset](auto i) { return i + offset; });
                for (auto i : absolute_atom_index) {
                    cell_list->updateMemberAt(i, spc.particles.at(i).pos);
                }
            }
        } else {
            auto isActive = getGroupFilter<Group<Particle>::ACTIVE>();
            for (const auto& group_change : change.groups) {
                const auto& group = spc.groups.at(group_change.group_index);
                const auto offset = spc.getFirstParticleIndex(group);
                if (isActive(group)) { // if the group is active, there is at least one active particle in it
                    for (auto i : group_change.relative_atom_indices) {
                        if (i > std::distance(group.begin(), group.end()) // if index lies behind last active index
                            && cell_list->containsMember(i + offset)) {
                            cell_list->removeMember(i + offset);
                        } else if( !cell_list->containsMember(i + offset) ) {
                            cell_list->insertMember(i + offset, spc.particles.at(i + offset).pos);
                        }
                    }
                } else { // if everything in the group is inactive remove all members
                    for (auto i : group_change.relative_atom_indices) {
                        if (cell_list->containsMember(i + offset)) {
                            cell_list->removeMember(i + offset);
                        }
                    }
                }
            }
        }*/
    }
};

} // namespace Faunus

#endif // FAUNUS_SASA_H