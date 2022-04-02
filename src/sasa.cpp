
#include "sasa.h"
#include "space.h"
#include <range/v3/view/zip.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <doctest/doctest.h>

namespace Faunus {

namespace SASA {

/**
 * @brief updates sasa of target particles
 * @param neighbours_data array of NeighbourData objects
 * @param target_indices absolute indicies of target particles in ParticleVector
 */
void SASABase::updateSASA(const std::vector<SASA::Neighbours>& neighbours,
                          const std::vector<index_type>& target_indices) {
    // here is a potential place for parallelization?
    //#pragma OMP parallel num_threads(2)
    //{
    for (const auto& [neighbour, index] : ranges::views::zip(neighbours, target_indices)) {
        areas.at(index) = calcSASAOfParticle(neighbour);
    }
    //}
}

/**
 * @brief Calcuates SASA of a single particle defined by NeighbourData object
 * @details cuts a sphere in z-direction, for each slice, radius of circle_i in the corresponding z-plane is
   calculated
 * then for each neighbour, calculate the overlaping part of circle_i with neighbouring circle_j and add these
 * arcs into vector, finally from this vector, calculate the exposed part of circle_i
 * @param neighbours NeighbourData object of given particle
 */
double SASABase::calcSASAOfParticle(const SASABase::Neighbours& neighbour) const {
    const auto sasa_radius_i = sasa_radii.at(neighbour.index);
    double area(0.);

    auto slice_height = 2. * sasa_radius_i / slices_per_atom;
    auto z = -sasa_radius_i - 0.5 * slice_height;

    for (int islice = 0; islice != slices_per_atom; ++islice) {
        z += slice_height;
        const auto sqrd_circle_radius_i = sasa_radius_i * sasa_radius_i - z * z;
        if (sqrd_circle_radius_i < 0.) {
            continue;
        }

        const auto circle_radius_i = std::sqrt(sqrd_circle_radius_i);
        if (circle_radius_i <= 0.) {
            continue;
        } /* round-off errors */

        std::vector<std::pair<double, double>> arcs;
        bool is_buried = false;
        for (const auto& [d_r, neighbour_index] : ranges::views::zip(neighbour.points, neighbour.indices)) {
            const auto sasa_radius_j = sasa_radii.at(neighbour_index);
            const auto z_distance = std::fabs(d_r.z() - z);

            if (z_distance < sasa_radius_j) {
                const auto sqrd_circle_radius_j = sasa_radius_j * sasa_radius_j - z_distance * z_distance;
                const auto circle_radius_j = std::sqrt(sqrd_circle_radius_j);
                const auto xy_distance = std::sqrt(d_r.x() * d_r.x() + d_r.y() * d_r.y());
                if (xy_distance >= circle_radius_i + circle_radius_j) { /* atoms aren't in contact */
                    continue;
                }
                if (xy_distance + circle_radius_i < circle_radius_j) { /* circle i is completely inside j */
                    is_buried = true;
                    break;
                }
                if (xy_distance + circle_radius_j < circle_radius_i) { /* circle j is completely inside i */
                    continue;
                }
                /* arc of circle i intersected by circle j */
                auto intersected_arc_halfsize =
                    std::acos((sqrd_circle_radius_i + xy_distance * xy_distance - sqrd_circle_radius_j) /
                              (2.0 * circle_radius_i * xy_distance));
                /* position of mid-point of intersection along circle i */
                auto intersection_midpoint_angle = std::atan2(d_r.y(), d_r.x()) + M_PI;
                auto beginning_arc_angle = intersection_midpoint_angle - intersected_arc_halfsize;
                auto end_arc_angle = intersection_midpoint_angle + intersected_arc_halfsize;
                if (beginning_arc_angle < 0) {
                    beginning_arc_angle += two_pi;
                }
                if (end_arc_angle > two_pi) {
                    end_arc_angle -= two_pi;
                }
                /* store the arc, if arc passes 2*PI split into two */
                if (end_arc_angle < beginning_arc_angle) {
                    /* store arcs as pairs of angles */
                    arcs.insert(arcs.end(), {{0.0, end_arc_angle}, {beginning_arc_angle, two_pi}});
                } else {
                    arcs.emplace_back(beginning_arc_angle, end_arc_angle);
                }
            }
        }
        if (!is_buried) {
            area += slice_height * sasa_radius_i * exposedArcLength(arcs);
        }
    }
    return area;
}

/**
 * @brief Calcuates total arc length in radians of overlapping arcs defined by two angles
 * @param vector of arcs, defined by a pair (first angle, second angle)
 */
double SASABase::exposedArcLength(std::vector<std::pair<double, double>>& arcs) const {
    if (arcs.empty()) {
        return two_pi;
    }

    std::sort(arcs.begin(), arcs.end(), [](auto& a, auto& b) { return a.first < b.first; });

    auto total_arc_angle = arcs.at(0).first;
    auto end_arc_angle = arcs.at(0).second;

    std::for_each(arcs.begin() + 1, arcs.end(), [&](const auto& arc) {
        if (end_arc_angle < arc.first) {
            total_arc_angle += arc.first - end_arc_angle;
        }
        if (arc.second > end_arc_angle) {
            end_arc_angle = arc.second;
        }
    });
    return total_arc_angle + two_pi - end_arc_angle;
}

const std::vector<double>& SASABase::getAreas() const { return areas; }

/**
 * @param spc
 * @param probe_radius in angstrom
 * @param slices_per_atom number of slices of each sphere in sasa calculations
 */
SASABase::SASABase(const Space& spc, double probe_radius, int slices_per_atom)
    : probe_radius(probe_radius)
    , slices_per_atom(slices_per_atom)
    , first_particle(std::addressof(spc.particles.at(0))) {}

/**
 * @brief resizes areas buffer to size of ParticleVector and fill radii buffer with radii
 * @param space
 */
void SASA::init(const Space& spc) {
    using namespace ranges::cpp20::views;
    auto sasa_radius_l = [this](const Particle& particle) { return 0.5 * particle.traits().sigma + probe_radius; };
    sasa_radii = spc.particles | ranges::cpp20::views::transform(sasa_radius_l) | ranges::to<std::vector>;
    areas.resize(spc.particles.size());
}

/**
 * @brief calculates neighbourData object of a target particle specified by target indiex in ParticleVector
 * @brief using the naive O(N) neighbour search for a given target particle
 * @param space
 * @param target_index indicex of target particle in ParticleVector
 */
SASA::Neighbours SASA::calcNeighbourDataOfParticle(const Space& spc, const index_type target_index) const {

    SASA::Neighbours neighbour;

    const auto& particle_i = spc.particles.at(target_index);
    const auto sasa_radius_i = sasa_radii.at(target_index);
    neighbour.index = target_index;

    for (const auto& particle_j : spc.activeParticles()) {
        const auto neighbour_index = indexOf(particle_j);
        const auto sasa_radius_j = sasa_radii.at(neighbour_index);
        const auto sq_cutoff = (sasa_radius_i + sasa_radius_j) * (sasa_radius_i + sasa_radius_j);

        if (target_index != neighbour_index && spc.geometry.sqdist(particle_i.pos, particle_j.pos) < sq_cutoff) {
            const auto dr = spc.geometry.vdist(particle_i.pos, particle_j.pos);
            neighbour.points.push_back(dr);
            neighbour.indices.push_back(neighbour_index);
        }
    }
    return neighbour;
}

/**
 * @brief calculates neighbourData objects of particles specified by target indices in ParticleVector
 * @param space
 * @param target_indices absolute indicies of target particles in ParticleVector
 */
std::vector<SASA::Neighbours> SASA::calcNeighbourData(const Space& spc,
                                                      const std::vector<index_type>& target_indices) const {

    return target_indices |
           ranges::views::transform([&](auto index) { return calcNeighbourDataOfParticle(spc, index); }) |
           ranges::to<std::vector>;
}

/**
 * @param spc
 * @param probe_radius in angstrom
 * @param slices_per_atom number of slices of each sphere in sasa calculations
 */
SASA::SASA(const Space& spc, double probe_radius, int slices_per_atom)
    : SASABase(spc, probe_radius, slices_per_atom) {}
SASA::SASA(const json& j, Space& spc) : SASABase(spc, j.value("radius", 1.4) * 1.0_angstrom, j.value("slices", 20)) {}

/**
 * @brief updates radii vector in case of matter change
 * @param space
 * @param change
 */
void SASA::update([[maybe_unused]] const Space& spc, [[maybe_unused]] const Change& change) {}

TEST_CASE("[Faunus] SASAPBC") {
    using doctest::Approx;
    Change change; // change object telling that a full energy calculation
    change.everything = true;
    pc::temperature = 300.0_K;
    atoms = R"([
        { "A": { "sigma": 4.0, "tfe": 1.0 } },
        { "B": { "sigma": 2.4, "tfe": 1.0 } }
    ])"_json.get<decltype(atoms)>();
    molecules = R"([
        { "M": { "atoms": ["A", "B"], "atomic": true } }
    ])"_json.get<decltype(molecules)>();
    json j = R"({
        "geometry": {"type": "cuboid", "length": [200.0, 200.0, 200.0] },
        "insertmolecules": [ { "M": { "N": 1 } } ]
    })"_json;
    Space spc = j;

    SASA sasa(spc, 1.4_angstrom, 20);
    const std::vector<double> radii = {4.0 * 0.5, 2.4 * 0.5};
    sasa.init(spc);

    SUBCASE("not intersecting") {

        spc.particles.at(0).pos = {30.0, 0.0, 0.0};
        spc.particles.at(1).pos = {5.0, 0.0, 0.0};

        const auto& neighbours = sasa.calcNeighbourData(spc, {0, 1});
        sasa.updateSASA(neighbours, {0, 1});
        const auto& areas = sasa.getAreas();

        CHECK(areas[0] == Approx(3.4 * 3.4 * M_PI * 4.));
    }

    SUBCASE("intersecting") {

        spc.particles.at(0).pos = {7.0, 0.0, 0.0};
        spc.particles.at(1).pos = {5.0, 0.0, 0.0};

        const auto& neighbours = sasa.calcNeighbourData(spc, {0, 1});
        sasa.updateSASA(neighbours, {0, 1});
        const auto& areas = sasa.getAreas();

        CHECK(areas[0] == Approx(119.48260171150575));
    }

    SUBCASE("intersecting accross boundary") {

        spc.particles.at(0).pos = {199., 0.0, 0.0};
        spc.particles.at(1).pos = {1., 0.0, 199.};

        const auto& neighbours = sasa.calcNeighbourData(spc, {0, 1});
        sasa.updateSASA(neighbours, {0, 1});
        const auto& areas = sasa.getAreas();

        CHECK(areas[0] == Approx(118.99710056237043));
    }

    SUBCASE("smaller buried in larger") {

        spc.particles.at(0).pos = {1.0, 0.0, 0.0};
        spc.particles.at(1).pos = {1.1, 0.0, 0.0};

        const auto& neighbours = sasa.calcNeighbourData(spc, {0, 1});
        sasa.updateSASA(neighbours, {0, 1});
        const auto& areas = sasa.getAreas();

        CHECK(areas[0] == Approx(3.4 * 3.4 * M_PI * 4));
        CHECK(areas[1] == Approx(0.));
    }
}

/**
 * @param spc
 * @param probe_radius in angstrom
 * @param slices_per_atom number of slices of each sphere in sasa calculations
 */
template <typename CellList>
SASACellList<CellList>::SASACellList(const Space& spc, double probe_radius, int slices_per_atom)
    : SASABase(spc, probe_radius, slices_per_atom) {}

template <typename CellList>
SASACellList<CellList>::SASACellList(const json& j, Space& spc)
    : SASABase(spc, j.value("radius", 1.4) * 1.0_angstrom, j.value("slices", 20)) {}

/**
 * @brief constructs cell_list with appropriate cell_length and fills it with particles from space
 * @param space
 */
template <typename CellList> void SASACellList<CellList>::init(const Space& spc) {
    cell_offsets.clear();
    for (auto i = -1; i <= 1; ++i) {
        for (auto j = -1; j <= 1; ++j) {
            for (auto k = -1; k <= 1; ++k) {
                cell_offsets.emplace_back(i, j, k);
            }
        }
    }

    auto get_sasa_radius = [this](auto& i) { return 0.5 * i.traits().sigma + probe_radius; };
    sasa_radii = spc.particles | ranges::cpp20::views::transform(get_sasa_radius) | ranges::to<std::vector>;
    const auto max_sasa_radius = ranges::cpp20::max(sasa_radii);

    cell_length = 2.0 * (max_sasa_radius);
    const auto active_particles = spc.activeParticles();
    createCellList(active_particles.begin(), active_particles.end(), spc.geometry);

    areas.resize(spc.particles.size(), 0.);
}

/**
 * @brief calculates neighbourData object of a target particle
 * @brief specified by target index in ParticleVector using cell list
 * @param space
 * @param target_index indicex of target particle in ParticleVector
 */
template <typename CellList>
SASABase::Neighbours SASACellList<CellList>::calcNeighbourDataOfParticle(const Space& spc,
                                                                         const index_type target_index) const {

    SASABase::Neighbours neighbours;

    const auto& particle_i = spc.particles.at(target_index);
    neighbours.index = target_index;
    const auto sasa_radius_i = sasa_radii.at(target_index);
    const auto& center_cell = cell_list->getGrid().coordinatesAt(particle_i.pos + 0.5 * spc.geometry.getLength());

    auto neighour_particles_at = [&](const CellCoord& offset) {
        return cell_list->getNeighborMembers(center_cell, offset);
    };

    for (const auto& cell_offset : cell_offsets) {
        const auto& neighbour_particle_indices = neighour_particles_at(cell_offset);
        for (const auto neighbour_particle_index : neighbour_particle_indices) {
            const auto& particle_j = spc.particles.at(neighbour_particle_index);
            const auto sasa_radius_j = sasa_radii.at(neighbour_particle_index);
            const auto sq_cutoff = (sasa_radius_i + sasa_radius_j) * (sasa_radius_i + sasa_radius_j);

            if (target_index != neighbour_particle_index &&
                spc.geometry.sqdist(particle_i.pos, particle_j.pos) <= sq_cutoff) {
                const auto dr = spc.geometry.vdist(particle_i.pos, particle_j.pos);
                neighbours.points.push_back(dr);
                neighbours.indices.push_back(neighbour_particle_index);
            }
        }
    }
    return neighbours;
}

/**
 * @brief calculates neighbourData objects of particles
 * @brief specified by target indices in ParticleVector using cell lists
 * @param space
 * @param target_indices absolute indicies of target particles in ParticleVector
 */
template <typename CellList>
std::vector<SASABase::Neighbours>
SASACellList<CellList>::calcNeighbourData(const Space& spc, const std::vector<index_type>& target_indices) const {
    return target_indices |
           ranges::views::transform([&](auto index) { return calcNeighbourDataOfParticle(spc, index); }) |
           ranges::to<std::vector>;
}

/**
 * @brief updates cell_list according to change, if the volume changes the cell_list gets rebuilt
 * @brief also updates radii in case of matter change
 * @param space
 * @param change
 */
template <typename CellList> void SASACellList<CellList>::update(const Space& spc, const Change& change) {
    if (change.everything || change.volume_change) {
        const auto active_particles = spc.activeParticles();
        createCellList(active_particles.begin(), active_particles.end(), spc.geometry);
    } else if (change.matter_change) {
        updateMatterChange(spc, change);
    } else {
        updatePositionsChange(spc, change);
    }
}

/**
 * @brief creates cell_list if it does not exist, otherwise resets its and fills it with active particles
 * @param space
 */
template <typename CellList>
template <typename TBegin, typename TEnd>
void SASACellList<CellList>::createCellList(TBegin begin, TEnd end, const GeometryType& geometry) {
    if (cell_list.get()) {
        delete cell_list.release();
        cell_list = std::make_unique<CellList>(geometry.getLength(), cell_length);
    } else {
        cell_list = std::make_unique<CellList>(geometry.getLength(), cell_length);
    }
    std::for_each(begin, end, [&, half_box = 0.5 * geometry.getLength()](const Particle& particle) {
        cell_list->insertMember(indexOf(particle), particle.pos + half_box);
    });
}

/**
 * @brief updates cell_list when particles got activated or desactivated
 * @param space
 * @param change
 */
template <typename CellList> void SASACellList<CellList>::updateMatterChange(const Space& spc, const Change& change) {
    for (const auto& group_change : change.groups) {
        const auto& group = spc.groups.at(group_change.group_index);
        const auto offset = spc.getFirstParticleIndex(group);
        for (const auto relative_index : group_change.relative_atom_indices) {
            const auto absolute_index = relative_index + offset;
            if (relative_index >= group.size()) { // if index lies behind last active index
                cell_list->removeMember(absolute_index);
            } else if (relative_index < group.size()) {
                cell_list->insertMember(absolute_index,
                                        spc.particles.at(absolute_index).pos + 0.5 * spc.geometry.getLength());
            }
        }
    }
}

/**
 * @brief updates cell_list when particles moved only
 * @param space
 * @param change
 */
template <typename CellList>
void SASACellList<CellList>::updatePositionsChange(const Space& spc, const Change& change) {
    for (const auto& group_change : change.groups) {
        const auto& group = spc.groups.at(group_change.group_index);
        const auto offset = spc.getFirstParticleIndex(group);

        auto update = [&, half_box = 0.5 * spc.geometry.getLength()](auto index) {
            cell_list->updateMemberAt(index, spc.particles.at(index).pos + half_box);
        };

        if (group_change.relative_atom_indices.empty()) {
            const auto changed_atom_indices = ranges::cpp20::views::iota(offset, offset + group.size());
            ranges::cpp20::for_each(changed_atom_indices, update);
        } else {
            const auto changed_atom_indices = group_change.relative_atom_indices |
                                              ranges::cpp20::views::transform([offset](auto i) { return i + offset; });
            ranges::cpp20::for_each(changed_atom_indices, update);
        }
    }
}

template class SASACellList<DensePeriodicCellList>;
template class SASACellList<DenseFixedCellList>;
template class SASACellList<SparsePeriodicCellList>;
template class SASACellList<SparseFixedCellList>;
} // namespace SASA

TEST_CASE("[Faunus] SASA_CellList") {
    using namespace SASA;
    using doctest::Approx;
    Change change;
    change.everything = true;
    pc::temperature = 300.0_K;
    atoms = R"([
        { "A": { "sigma": 4.0, "tfe": 1.0 } },
        { "B": { "sigma": 2.4, "tfe": 1.0 } }
    ])"_json.get<decltype(atoms)>();
    molecules = R"([
        { "M": { "atoms": ["A", "B"], "atomic": true } }
    ])"_json.get<decltype(molecules)>();

    SUBCASE("periodic boundary") {
        json j = R"({
        "geometry": {"type": "cuboid", "length": [50.0, 50.0, 50.0] },
        "insertmolecules": [ { "M": { "N": 1 } } ]
        })"_json;
        Space spc = j;

        spc.particles.at(0).pos = {49.0, 0.0, 0.0};
        spc.particles.at(1).pos = {1.0, 0.0, 0.0};

        SASACellList<SparsePeriodicCellList> sasa(spc, 1.4_angstrom, 20);
        sasa.init(spc);

        auto neighbours = sasa.calcNeighbourData(spc, {0, 1});
        CHECK_EQ(neighbours[0].indices[0], 1);
        CHECK_EQ(neighbours[1].indices[0], 0);

        spc.groups[0].deactivate(spc.groups[0].begin(), spc.groups[0].begin() + 1);
        Change::GroupChange changed_data;
        changed_data.group_index = 0;
        changed_data.relative_atom_indices = {1};
        change.groups.push_back(changed_data);
        change.matter_change = 1;
        change.everything = 0;

        sasa.update(spc, change);

        neighbours = sasa.calcNeighbourData(spc, {0, 1});
        CHECK(neighbours[0].indices.empty());
        CHECK_EQ(neighbours[1].indices[0], 0);
    }
    SUBCASE("fixed boundary") {
        json j = R"({
        "geometry": {"type": "sphere", "radius": 50.0 },
        "insertmolecules": [ { "M": { "N": 1 } } ]
        })"_json;
        Space spc = j;

        spc.particles.at(0).pos = {49.0, 0.0, 0.0};
        spc.particles.at(1).pos = {1.0, 0.0, 0.0};

        SASACellList<SparseFixedCellList> sasa(spc, 1.4_angstrom, 20);
        sasa.init(spc);

        auto neighbours = sasa.calcNeighbourData(spc, {0, 1});
        CHECK(neighbours[0].indices.empty());
        CHECK(neighbours[1].indices.empty());
    }
}

} // namespace Faunus
