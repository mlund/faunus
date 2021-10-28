
#include "sasa.h"
#include <range/v3/view/zip.hpp>

namespace Faunus {

void SASA::init(Space& spc) {
    // make areas vector with a size of ALL particles in ParticleVector
    areas.resize(spc.particles.size());
}

void SASA::updateSASA(const std::vector<SASA::NeighboursData>& neighbours_data, const std::vector<double>& radii,
                      const std::vector<int>& target_indices) {

    // here is a potential place for parallelization?
    // #pragma OMP parallel num_threads(2)
    // {
    for (const auto& [neighbour_data, index] : ranges::views::zip(neighbours_data, target_indices)) {
        areas[index] = calcSASAOfParticle(neighbour_data, radii);
    }
    //  }
}

SASA::NeighboursData SASA::calcNeighbourDataOfParticle(Space& spc, const int target_index) {

    // O(N^2) search for neighbours
    SASA::NeighboursData neighbour_data;

    const auto& particle_i = spc.particles.at(target_index);
    const double rc_i = particle_i.traits().sigma * 0.5 + probe_radius;
    neighbour_data.index = target_index;

    for (const auto& particle_j : spc.activeParticles()) {
        const double rc_j = particle_j.traits().sigma * 0.5 + probe_radius;
        const auto sq_cutoff = (rc_i + rc_j) * (rc_i + rc_j);
        const auto neighbour_index = std::addressof(particle_j) - std::addressof(spc.particles[0]);

        if (spc.geometry.sqdist(particle_i.pos, particle_j.pos) < sq_cutoff && target_index != neighbour_index) {
            neighbour_data.neighbour_indices.push_back(neighbour_index);

            Point dr = particle_i.pos - particle_j.pos;
            spc.geometry.boundary(dr);
            neighbour_data.points.push_back(dr);
        }
    }
    return neighbour_data;
}

const std::vector<SASA::NeighboursData> SASA::calcNeighbourData(Space& spc, const std::vector<int>& target_indices) {

    // O(N^2) search for neighbours ... will be done using Cell-Lists
    const auto number_of_indices = target_indices.size();
    std::vector<SASA::NeighboursData> neighbour_data(number_of_indices);

    for (size_t i = 0; i != number_of_indices; ++i) {
        neighbour_data.at(i) = calcNeighbourDataOfParticle(spc, target_indices.at(i));
    }

    return neighbour_data;
}

//!< slices a sphere in z-direction, for each slice, radius of circle_i in the corresponding z-plane is calculated
//!< then for each neighbour, calculate the overlaping part of circle_i with neighbouring circle_j and add these
//!< arcs into vector, finally from this vector, calculate the exposed part of circle_i
double SASA::calcSASAOfParticle(const SASA::NeighboursData& neighbour_data, const std::vector<double>& radii) const {

    const double particle_radius_i = radii[neighbour_data.index] + probe_radius;
    double area(0.);

    double slice_height = 2. * particle_radius_i / n_slices_per_atom;
    double z = -particle_radius_i - 0.5 * slice_height;

    for (int islice = 0; islice != n_slices_per_atom; ++islice) {
        z += slice_height;
        const double sqrd_circle_radius_i = particle_radius_i * particle_radius_i - z * z;
        if (sqrd_circle_radius_i < 0) {
            continue;
        }

        const double circle_radius_i = std::sqrt(sqrd_circle_radius_i);
        if (circle_radius_i <= 0) {
            continue;
        } /* round-off errors */

        std::vector<std::pair<double, double>> arcs;
        bool is_buried = false;
        for (const auto& [d_r, neighbour_index] :
             ranges::views::zip(neighbour_data.points, neighbour_data.neighbour_indices)) {
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
                double alpha = std::acos((sqrd_circle_radius_i + xy_distance * xy_distance - sqrd_circle_radius_j) /
                                         (2.0 * circle_radius_i * xy_distance));
                /* position of mid-point of intersection along circle i */
                double beta = std::atan2(d_r.y(), d_r.x()) + M_PI;
                double inf = beta - alpha;
                double sup = beta + alpha;
                if (inf < 0) {
                    inf += TWOPI;
                }
                if (sup > TWOPI) {
                    sup -= TWOPI;
                }
                /* store the arc, if arc passes 2*PI split into two */
                if (sup < inf) {
                    /* store arcs as pairs of angles */
                    arcs.insert(arcs.end(), {std::make_pair(0.0, sup), std::make_pair(inf, TWOPI)});
                } else {
                    arcs.insert(arcs.end(), std::make_pair(inf, sup));
                }
            }
        }
        if (!is_buried) {
            area += slice_height * particle_radius_i * exposedArcLength(arcs);
        }
    }
    return area;
}

double SASA::exposedArcLength(std::vector<std::pair<double, double>>& arcs) const {

    if (arcs.empty()) {
        return TWOPI;
    }

    auto sortByFirst = [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first;
    };
    std::sort(arcs.begin(), arcs.end(), sortByFirst);

    double sum = arcs[0].first;
    double sup = arcs[0].second;
    double tmp;

    for (size_t i = 1; i < arcs.size(); ++i) {
        if (sup < arcs[i].first) {
            sum += arcs[i].first - sup;
        }
        tmp = arcs[i].second;
        if (tmp > sup) {
            sup = tmp;
        }
    }
    return sum + TWOPI - sup;
}

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
        sasa.updateSASA(neighbours, radii, {0, 1});
        const auto& areas = sasa.getAreas();

        CHECK(areas[0] == Approx(3.4 * 3.4 * M_PI * 4.));
    }

    SUBCASE("intersecting") {

        spc.particles.at(0).pos = {7.0, 0.0, 0.0};
        spc.particles.at(1).pos = {5.0, 0.0, 0.0};

        const auto& neighbours = sasa.calcNeighbourData(spc, {0, 1});
        sasa.updateSASA(neighbours, radii, {0, 1});
        const auto& areas = sasa.getAreas();

        CHECK(areas[0] == Approx(119.48260171150575));
    }

    SUBCASE("intersecting accross boundary") {

        spc.particles.at(0).pos = {199., 0.0, 0.0};
        spc.particles.at(1).pos = {1., 0.0, 199.};

        const auto& neighbours = sasa.calcNeighbourData(spc, {0, 1});
        sasa.updateSASA(neighbours, radii, {0, 1});
        const auto& areas = sasa.getAreas();

        CHECK(areas[0] == Approx(118.99710056237043));
    }

    SUBCASE("smaller buried in larger") {

        spc.particles.at(0).pos = {1.0, 0.0, 0.0};
        spc.particles.at(1).pos = {1.1, 0.0, 0.0};

        const auto& neighbours = sasa.calcNeighbourData(spc, {0, 1});
        sasa.updateSASA(neighbours, radii, {0, 1});
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
