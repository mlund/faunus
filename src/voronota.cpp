#include "analysis.h"
#include "mpicontroller.h"
#include <voronotalt/voronotalt.h>
#include <numeric>

namespace Faunus::analysis {

Voronota::Voronota(double probe_radius, const Faunus::Space& spc)
    : Analysis(spc, "voronota")
    , probe_radius(probe_radius) {
    cite = "doi:10/mq8k";
    auto n_pbc = spc.geometry.asSimpleGeometry()->boundary_conditions.isPeriodic().count();
    switch (n_pbc) {
    case 0:
        faunus_logger->debug("{}: No PBC detected", name);
        use_pbc = false;
        break;
    case 3:
        faunus_logger->debug("{}: 3D PBC detected", name);
        use_pbc = true;
        break;
    default:
        faunus_logger->warn("{}: Non-uniform PBC is currently ignored - be careful!", name);
        use_pbc = false;
    }
}

Voronota::Voronota(const Faunus::json& input, const Faunus::Space& spc)
    : Voronota(input.value("radius", 1.4_angstrom), spc) {
    from_json(input);
}

void Voronota::_from_json(const Faunus::json& input) {
    if (filename = input.value("file", ""s); !filename.empty()) {
        output_stream = IO::openCompressedOutputStream(MPI::prefix + filename);
        *output_stream << "# step SASA\n";
    }
}

void Voronota::_to_json(json& json_output) const {
    if (!average_data.area.empty()) {
        json_output = {{"⟨SASA⟩", average_data.area.avg()},
                      {"⟨SASA²⟩-⟨SASA⟩²", average_data.area_squared.avg() - std::pow(average_data.area.avg(), 2)}};
    }
    json_output["radius"] = probe_radius;
}

void Voronota::_to_disk() {
    if (output_stream) {
        output_stream->flush();
    }
}

void Voronota::_sample() {
    using voronotalt::SimplePoint;
    using namespace ranges::cpp20::views;

    // Convert single `Particle` to Voronota's `SimpleSphere`
    auto to_sphere = [&](const Particle& p) -> voronotalt::SimpleSphere {
        return {{p.pos.x(), p.pos.y(), p.pos.z()}, (0.5 * p.traits().sigma) + probe_radius};
    };
    const auto spheres = spc.activeParticles() | transform(to_sphere) | ranges::to_vector;

    voronotalt::RadicalTessellation::Result result;

    if (use_pbc) {
        auto to_point = [](const Point& p) -> SimplePoint { return {p.x(), p.y(), p.z()}; };
        const Point corner = 0.5 * spc.geometry.getLength();
        const std::vector<SimplePoint> box_corners = {to_point(-corner), to_point(corner)};
        voronotalt::RadicalTessellation::construct_full_tessellation(spheres, box_corners, result);
    } else {
        voronotalt::RadicalTessellation::construct_full_tessellation(spheres, result);
    }

    auto areas = result.cells_summaries | transform([](auto& summary) { return summary.sas_area; });
    const auto total_area = std::accumulate(areas.begin(), areas.end(), 0.0);
    average_data.area.add(total_area);
    average_data.area_squared.add(total_area * total_area);

    if (output_stream) {
        *output_stream << this->getNumberOfSteps() << " " << total_area << "\n";
    }
}

} // namespace Faunus::analysis
