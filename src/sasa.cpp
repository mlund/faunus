
#include "sasa.h"
#include <range/v3/view/zip.hpp>

namespace Faunus {

namespace SASA {

void SASABase::updateSASA(const std::vector<SASA::Neighbours>& neighbours,
                          const std::vector<AtomIndex>& target_indices) {

    // here is a potential place for parallelization?
    // #pragma OMP parallel num_threads(2)
    // {
    for (const auto& [neighbour, index] : ranges::views::zip(neighbours, target_indices)) {
        areas[index] = calcSASAOfParticle(neighbour);
    }
    //  }
}

//!< slices a sphere in z-direction, for each slice, radius of circle_i in the corresponding z-plane is calculated
//!< then for each neighbour, calculate the overlaping part of circle_i with neighbouring circle_j and add these
//!< arcs into vector, finally from this vector, calculate the exposed part of circle_i
double SASABase::calcSASAOfParticle(const SASABase::Neighbours& neighbour) const {
    const auto particle_radius_i = radii[neighbour.index] + probe_radius;
    double area(0.);

    auto slice_height = 2. * particle_radius_i / slices_per_atom;
    auto z = -particle_radius_i - 0.5 * slice_height;

    for (int islice = 0; islice != slices_per_atom; ++islice) {
        z += slice_height;
        const auto sqrd_circle_radius_i = particle_radius_i * particle_radius_i - z * z;
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
            const auto particle_radius_j = radii[neighbour_index] + probe_radius;
            const auto z_distance = std::fabs(d_r.z() - z);

            if (z_distance < particle_radius_j) {
                const auto sqrd_circle_radius_j = particle_radius_j * particle_radius_j - z_distance * z_distance;
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
                    beginning_arc_angle += TWOPI;
                }
                if (end_arc_angle > TWOPI) {
                    end_arc_angle -= TWOPI;
                }
                /* store the arc, if arc passes 2*PI split into two */
                if (end_arc_angle < beginning_arc_angle) {
                    /* store arcs as pairs of angles */
                    arcs.insert(arcs.end(), {{0.0, end_arc_angle}, {beginning_arc_angle, TWOPI}});
                } else {
                    arcs.insert(arcs.end(), {beginning_arc_angle, end_arc_angle});
                }
            }
        }
        if (!is_buried) {
            area += slice_height * particle_radius_i * exposedArcLength(arcs);
        }
    }
    return area;
}

double SASABase::exposedArcLength(std::vector<std::pair<double, double>>& arcs) const {
    if (arcs.empty()) {
        return TWOPI;
    }

    std::sort(arcs.begin(), arcs.end(), [](auto& a, auto& b) { return a.first < b.first; });

    auto total_arc_angle = arcs[0].first;
    auto end_arc_angle = arcs[0].second;

    std::for_each(arcs.begin() + 1, arcs.end(), [&](const auto& arc) {
        if (end_arc_angle < arc.first) {
            total_arc_angle += arc.first - end_arc_angle;
        }
        if (arc.second > end_arc_angle) {
            end_arc_angle = arc.second;
        }
    });
    return total_arc_angle + TWOPI - end_arc_angle;
}

const std::vector<double>& SASABase::getAreas() const { return areas; }

SASABase::SASABase(Space& spc, double probe_radius, int slices_per_atom)
    : probe_radius(probe_radius), slices_per_atom(slices_per_atom),
      first_particle(std::addressof(spc.particles.at(0U))) {}

void SASA::init(Space& spc) {
    radii.clear();
    std::for_each(spc.particles.begin(), spc.particles.end(),
                  [&](const Particle& particle) { radii.push_back(particle.traits().sigma * 0.5); });
    areas.resize(spc.particles.size());
}

SASA::Neighbours SASA::calcNeighbourDataOfParticle(Space& spc, const AtomIndex target_index) {

    // O(N^2) search for neighbours
    SASA::Neighbours neighbour;

    const auto& particle_i = spc.particles.at(target_index);
    const auto sasa_radius_i = radii[target_index] + probe_radius;
    neighbour.index = target_index;

    for (const auto& particle_j : spc.activeParticles()) {
        const auto neighbour_index = indexOf(particle_j);
        const auto sasa_radius_j = radii[neighbour_index] + probe_radius;
        const auto sq_cutoff = (sasa_radius_i + sasa_radius_j) * (sasa_radius_i + sasa_radius_j);

        if (target_index != neighbour_index && spc.geometry.sqdist(particle_i.pos, particle_j.pos) < sq_cutoff) {

            const auto dr = spc.geometry.vdist(particle_i.pos, particle_j.pos);
            neighbour.points.push_back(dr);
            neighbour.indices.push_back(neighbour_index);
        }
    }
    return neighbour;
}

std::vector<SASA::Neighbours> SASA::calcNeighbourData(Space& spc, const std::vector<AtomIndex>& target_indices) {

    // O(N^2) search for neighbours ... will be done using Cell-Lists
    const auto number_of_indices = static_cast<AtomIndex>(target_indices.size());
    std::vector<SASA::Neighbours> neighbour(number_of_indices);

    for (size_t i = 0U; i != number_of_indices; ++i) {
        neighbour.at(i) = calcNeighbourDataOfParticle(spc, target_indices.at(i));
    }

    return neighbour;
}

SASA::SASA(Space& spc, double probe_radius, int slices_per_atom) : SASABase(spc, probe_radius, slices_per_atom) {}
SASA::SASA(const json& j, Space& spc) : SASABase(spc, j.value("radius", 1.4) * 1.0_angstrom, j.value("slices", 20)) {}

void SASA::update([[maybe_unused]] Space& spc, [[maybe_unused]] const Change& change) {}

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

template <typename CellList_T>
SASACellList<CellList_T>::SASACellList(Space& spc, double probe_radius, int slices_per_atom)
    : SASABase(spc, probe_radius, slices_per_atom) {}

template <typename CellList_T>
SASACellList<CellList_T>::SASACellList(const json& j, Space& spc)
    : SASABase(spc, j.value("radius", 1.4) * 1.0_angstrom, j.value("slices", 20)) {}

template <typename CellList_T> void SASACellList<CellList_T>::init(Space& spc) {

    cell_offsets.clear();
    for (auto i = -1; i <= 1; ++i) {
        for (auto j = -1; j <= 1; ++j) {
            for (auto k = -1; k <= 1; ++k) {
                cell_offsets.push_back(CellCoord(i, j, k));
            }
        }
    }

    radii.clear();
    std::for_each(spc.particles.begin(), spc.particles.end(),
                  [&](const Particle& particle) { radii.push_back(particle.traits().sigma * 0.5); });
    auto max_radius = ranges::max(radii);

    cell_length = 2.0 * (max_radius + probe_radius);
    createCellList(spc);

    areas.resize(spc.particles.size(), 0.);
}

template <typename CellList_T>
SASABase::Neighbours SASACellList<CellList_T>::calcNeighbourDataOfParticle(Space& spc, const AtomIndex target_index) {

    SASABase::Neighbours neighbours;

    const auto& particle_i = spc.particles.at(target_index);
    neighbours.index = target_index;
    const auto sasa_radius_i = radii[target_index] + probe_radius;
    const auto& center_cell = cell_list->getGrid().coordinatesAt(particle_i.pos + spc.geometry.getLength() / 2.);

    auto neighour_particles_at = [&](const CellCoord& offset) {
        return cell_list->getNeighborMembers(center_cell, offset);
    };

    for (const auto& cell_offset : cell_offsets) {

        const auto& neighbour_particle_indices = neighour_particles_at(cell_offset);

        for (const auto neighbour_particle_index : neighbour_particle_indices) {

            const auto& particle_j = spc.particles.at(neighbour_particle_index);
            const auto sasa_radius_j = radii[neighbour_particle_index] + probe_radius;
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

template <typename CellList_T>
std::vector<SASABase::Neighbours>
SASACellList<CellList_T>::calcNeighbourData(Space& spc, const std::vector<AtomIndex>& target_indices) {

    const auto number_of_indices = target_indices.size();
    std::vector<SASA::Neighbours> neighbours(number_of_indices);

    for (size_t i = 0; i != number_of_indices; ++i) {
        neighbours.at(i) = calcNeighbourDataOfParticle(spc, target_indices.at(i));
    }

    return neighbours;
}
template <typename CellList_T> void SASACellList<CellList_T>::update(Space& spc, const Change& change) {

    if (change.everything || change.volume_change) {
        createCellList(spc);
    } else if (change.matter_change) {
        updateMatterChange(spc, change);
    } else {
        for (const auto& group_change : change.groups) {
            const auto& group = spc.groups.at(group_change.group_index);
            const auto offset = spc.getFirstParticleIndex(group);
            if (group_change.relative_atom_indices.empty()) {
                const auto changed_atom_indices = ranges::cpp20::views::iota(offset, offset + group.size());
                for (const auto changed_index : changed_atom_indices) {
                    cell_list->updateMemberAt(changed_index,
                                              spc.particles.at(changed_index).pos + spc.geometry.getLength() / 2.);
                }
            } else {
                const auto changed_atom_indices =
                    group_change.relative_atom_indices |
                    ranges::cpp20::views::transform([offset](auto i) { return i + offset; });
                for (const auto changed_index : changed_atom_indices) {
                    cell_list->updateMemberAt(changed_index,
                                              spc.particles.at(changed_index).pos + spc.geometry.getLength() / 2.);
                }
            }
        }
    }
}

template <typename CellList_T> void SASACellList<CellList_T>::createCellList(Space& spc) {

    cell_list = std::make_unique<CellList_T>(spc.geometry.getLength(), cell_length);
    for (const auto& particle : spc.activeParticles()) {
        const auto particle_index = indexOf(particle);
        cell_list->insertMember(particle_index, particle.pos + spc.geometry.getLength() / 2.);
    }
}

template <typename CellList_T> void SASACellList<CellList_T>::updateMatterChange(Space& spc, const Change& change) {
    const auto is_active = getGroupFilter<Group::Selectors::ACTIVE>();
    for (const auto& group_change : change.groups) {
        const auto& group = spc.groups.at(group_change.group_index);
        const auto offset = spc.getFirstParticleIndex(group);
        for (const auto relative_index : group_change.relative_atom_indices) {
            const auto absolute_index = relative_index + offset;
            if (relative_index >= group.size() &&
                cell_list->containsMember(absolute_index)) { // if index lies behind last active index
                cell_list->removeMember(absolute_index);
            } else if (relative_index < group.size() && !cell_list->containsMember(absolute_index)) {
                cell_list->insertMember(absolute_index,
                                        spc.particles.at(absolute_index).pos + spc.geometry.getLength() / 2.);
            }
        }
    }
}

template class SASACellList<PeriodicCellList>;
template class SASACellList<FixedCellList>;
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

        SASACellList<PeriodicCellList> sasa(spc, 1.4_angstrom, 20);
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

        SASACellList<FixedCellList> sasa(spc, 1.4_angstrom, 20);
        sasa.init(spc);

        auto neighbours = sasa.calcNeighbourData(spc, {0, 1});
        CHECK(neighbours[0].indices.empty());
        CHECK(neighbours[1].indices.empty());
    }
}

} // namespace Faunus
