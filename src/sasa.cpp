
#include "sasa.h"
#include <range/v3/view/zip.hpp>

namespace Faunus {

void SASA::init(Space& spc) {
    // make areas vector with a size of ALL particles in ParticleVector
    areas.resize(spc.particles.size());
}

void SASA::updateSASA(const std::vector<SASA::NeighboursData>& neighbours, const std::vector<double>& radii,
                      const std::vector<int>& target_inds) {

    // here is a potential place for parallelization?
    // #pragma OMP parallel num_threads(2)
    // {
    for (const auto [neighbour, ind] : ranges::views::zip(neighbours, target_inds)) {
        double sasa = calcSASAOfParticle(neighbour, radii);
        areas[ind] = sasa;
    }
    //  }
}

SASA::NeighboursData SASA::calcNeighbourDataOfParticle(Space& spc, const int target_ind) {

    // O(N^2) search for neighbours
    SASA::NeighboursData neighbour_data;

    const auto& p_i = spc.particles.at(target_ind);
    const double rc_i = p_i.traits().sigma * 0.5 + probe_radius;
    neighbour_data.ind = target_ind;

    for (const auto& p_j : spc.activeParticles()) {
        const double rc_j = p_j.traits().sigma * 0.5 + probe_radius;
        const auto sq_cutoff = (rc_i + rc_j) * (rc_i + rc_j);
        const int j = static_cast<int>(std::addressof(p_j) - std::addressof(spc.particles[0]));

        if (spc.geometry.sqdist(p_i.pos, p_j.pos) < sq_cutoff && target_ind != j) {
            neighbour_data.neighbour_inds.push_back(j);

            Point dr = p_i.pos - p_j.pos;
            spc.geometry.boundary(dr);
            neighbour_data.points.push_back(dr);
        }
    }
    return neighbour_data;
}
const std::vector<SASA::NeighboursData> SASA::calcNeighbourData(Space& spc, const std::vector<int>& target_inds) {

    // O(N^2) search for neighbours ... will be done using Cell-Lists
    const int n_inds = target_inds.size();
    std::vector<SASA::NeighboursData> neighbour_data(n_inds);

    for (int i = 0; i != n_inds; ++i) {
        neighbour_data.at(i) = calcNeighbourDataOfParticle(spc, target_inds.at(i));
    }

    return neighbour_data;
}

double SASA::calcSASAOfParticle(const SASA::NeighboursData& neighbours, const std::vector<double>& radii) const {

    const int number_of_neighbours = neighbours.points.size();

    double dj, dij, ri_prime2, ri_prime, rj_prime, rj_prime2;
    const double ri = radii[neighbours.ind] + probe_radius;
    double area(0.);

    double delta = 2. * ri / n_slices_per_atom;
    double z = -ri - 0.5 * delta;

    for (int islice = 0; islice != n_slices_per_atom; ++islice) {
        z += delta;
        ri_prime2 = ri * ri - z * z;
        if (ri_prime2 < 0)
            continue;
        ri_prime = std::sqrt(ri_prime2);
        if (ri_prime <= 0)
            continue; /* round-off errors */
        std::vector<std::pair<double, double>> arc;
        bool is_buried = false;
        for (int j = 0; j != number_of_neighbours; ++j) {
            const Point d_r = neighbours.points[j];
            const double rj = radii[neighbours.neighbour_inds[j]] + probe_radius;
            dj = std::fabs(d_r.z() - z);

            if (dj < rj) {
                rj_prime2 = rj * rj - dj * dj;
                rj_prime = std::sqrt(rj_prime2);
                dij = std::sqrt(d_r.x() * d_r.x() + d_r.y() * d_r.y());
                if (dij >= ri_prime + rj_prime) { /* atoms aren't in contact */
                    continue;
                }
                if (dij + ri_prime < rj_prime) { /* circle i is completely inside j */
                    is_buried = true;
                    break;
                }
                if (dij + rj_prime < ri_prime) { /* circle j is completely inside i */
                    continue;
                }
                /* arc of circle i intersected by circle j */
                double alpha = std::acos((ri_prime2 + dij * dij - rj_prime2) / (2.0 * ri_prime * dij));
                /* position of mid-point of intersection along circle i */
                double beta = std::atan2(d_r.y(), d_r.x()) + M_PI;
                double inf = beta - alpha;
                double sup = beta + alpha;
                if (inf < 0)
                    inf += TWOPI;
                if (sup > TWOPI)
                    sup -= TWOPI;
                /* store the arc, if arc passes 2*PI split into two */
                if (sup < inf) {
                    /* store arcs as pairs of angles */
                    arc.insert(arc.end(), {std::make_pair(0.0, sup), std::make_pair(inf, TWOPI)});
                } else {
                    arc.insert(arc.end(), std::make_pair(inf, sup));
                }
            }
        }
        if (!is_buried) {
            area += delta * ri * exposedArcLength(arc);
        }
    }
    return area;
}

double SASA::exposedArcLength(std::vector<std::pair<double, double>>& arc_f) const {

    if (arc_f.empty()) {
        return TWOPI;
    }

    auto sortByFirst = [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first;
    };
    std::sort(arc_f.begin(), arc_f.end(), sortByFirst);

    double sum = arc_f[0].first;
    double sup = arc_f[0].second;
    double tmp;

    for (int i = 1; i < arc_f.size(); ++i) {
        if (sup < arc_f[i].first) {
            sum += arc_f[i].first - sup;
        }
        tmp = arc_f[i].second;
        if (tmp > sup) {
            sup = tmp;
        }
    }

    return sum + TWOPI - sup;
}

TEST_CASE("[Faunus] SASAPBC2") {
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
