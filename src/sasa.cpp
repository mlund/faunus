
#include "sasa.h"
#include <range/v3/view/zip.hpp>

namespace Faunus {

void SASABase::updateSASA(const std::vector<SASA::Neighbours>& neighbours, const std::vector<size_t>& target_indices) {

    // here is a potential place for parallelization?
    // #pragma OMP parallel num_threads(2)
    // {
    for (const auto& [neighbour, index] : ranges::views::zip(neighbours, target_indices)) {
        areas[index] = calcSASAOfParticle(neighbour);
    }
    //  }
}

void SASA::init(Space& spc) {
    radii.clear();
    std::for_each(spc.particles.begin(), spc.particles.end(),
                  [&](const Particle& particle) { radii.push_back(particle.traits().sigma * 0.5); });
    areas.resize(spc.particles.size());
}

SASA::Neighbours SASA::calcNeighbourDataOfParticle(Space& spc, const size_t target_index) {

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

std::vector<SASA::Neighbours> SASA::calcNeighbourData(Space& spc, const std::vector<size_t>& target_indices) {

    // O(N^2) search for neighbours ... will be done using Cell-Lists
    const auto number_of_indices = static_cast<size_t>(target_indices.size());
    std::vector<SASA::Neighbours> neighbour(number_of_indices);

    for (size_t i = 0U; i != number_of_indices; ++i) {
        neighbour.at(i) = calcNeighbourDataOfParticle(spc, target_indices.at(i));
    }

    return neighbour;
}

SASA::SASA(Space& spc, double probe_radius, int slices_per_atom) : SASABase(spc, probe_radius, slices_per_atom) {}
SASA::SASA(const json& j, Space& spc) : SASABase(spc, j.value("radius", 1.4) * 1.0_angstrom, j.value("slices", 20)) {}

//!< slices a sphere in z-direction, for each slice, radius of circle_i in the corresponding z-plane is calculated
//!< then for each neighbour, calculate the overlaping part of circle_i with neighbouring circle_j and add these
//!< arcs into vector, finally from this vector, calculate the exposed part of circle_i
double SASABase::calcSASAOfParticle(const SASA::Neighbours& neighbour) const {

    const double particle_radius_i = radii[neighbour.index] + probe_radius;
    double area(0.);

    double slice_height = 2. * particle_radius_i / slices_per_atom;
    double z = -particle_radius_i - 0.5 * slice_height;

    for (int islice = 0; islice != slices_per_atom; ++islice) {
        z += slice_height;
        const double sqrd_circle_radius_i = particle_radius_i * particle_radius_i - z * z;
        if (sqrd_circle_radius_i < 0.) {
            continue;
        }

        const double circle_radius_i = std::sqrt(sqrd_circle_radius_i);
        if (circle_radius_i <= 0.) {
            continue;
        } /* round-off errors */

        std::vector<std::pair<double, double>> arcs;
        bool is_buried = false;
        for (const auto& [d_r, neighbour_index] : ranges::views::zip(neighbour.points, neighbour.indices)) {
            const double particle_radius_j = radii[neighbour_index] + probe_radius;
            const double z_distance = std::fabs(d_r.z() - z);

            if (z_distance < particle_radius_j) {
                const double sqrd_circle_radius_j = particle_radius_j * particle_radius_j - z_distance * z_distance;
                const double circle_radius_j = std::sqrt(sqrd_circle_radius_j);
                const double xy_distance = std::sqrt(d_r.x() * d_r.x() + d_r.y() * d_r.y());
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
                double intersected_arc_halfsize =
                    std::acos((sqrd_circle_radius_i + xy_distance * xy_distance - sqrd_circle_radius_j) /
                              (2.0 * circle_radius_i * xy_distance));
                /* position of mid-point of intersection along circle i */
                double intersection_midpoint_angle = std::atan2(d_r.y(), d_r.x()) + M_PI;
                double beginning_arc_angle = intersection_midpoint_angle - intersected_arc_halfsize;
                double end_arc_angle = intersection_midpoint_angle + intersected_arc_halfsize;
                if (beginning_arc_angle < 0) {
                    beginning_arc_angle += TWOPI;
                }
                if (end_arc_angle > TWOPI) {
                    end_arc_angle -= TWOPI;
                }
                /* store the arc, if arc passes 2*PI split into two */
                if (end_arc_angle < beginning_arc_angle) {
                    /* store arcs as pairs of angles */
                    arcs.insert(arcs.end(),
                                {std::make_pair(0.0, end_arc_angle), std::make_pair(beginning_arc_angle, TWOPI)});
                } else {
                    arcs.insert(arcs.end(), std::make_pair(beginning_arc_angle, end_arc_angle));
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

    auto sortByFirst = [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first;
    };
    std::sort(arcs.begin(), arcs.end(), sortByFirst);

    double total_arc_angle = arcs[0].first;
    double end_arc_angle = arcs[0].second;

    for (size_t i = 1; i < arcs.size(); ++i) {
        if (end_arc_angle < arcs[i].first) {
            total_arc_angle += arcs[i].first - end_arc_angle;
        }
        if (arcs[i].second > end_arc_angle) {
            end_arc_angle = arcs[i].second;
        }
    }
    return total_arc_angle + TWOPI - end_arc_angle;
}

SASABase::SASABase(Space& spc, double probe_radius, int slices_per_atom)
    : probe_radius(probe_radius), slices_per_atom(slices_per_atom),
      first_particle(std::addressof(spc.particles.at(0U))) {}

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

/*
TEST_CASE("[Faunus] SASAPBC_CellList") {
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

    spc.particles.at(0).pos = {30.0, 0.0, 0.0};
    spc.particles.at(1).pos = {5.0, 0.0, 0.0};
    SASACellList sasa(spc, 1.4_angstrom, 20);

    spc.particles.at(1).pos = {55.0, 0.0, 0.0};
    Change::GroupChange changed_data;
    changed_data.group_index = 0;
    changed_data.relative_atom_indices = {1};
    change.groups.push_back(changed_data);

    sasa.update(spc, change);

}
*/

} // namespace Faunus
