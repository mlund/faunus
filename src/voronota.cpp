#include "analysis.h"
#include "mpicontroller.h"
#include <voronotalt/voronotalt.h>
#include <numeric>

namespace Faunus::analysis {

Voronota::Voronota(double probe_radius, const Faunus::Space& spc)
    : Analysis(spc, "voronota")
    , probe_radius(probe_radius) {
    cite = "doi:10.1093/bioinformatics/btab448";
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

void Voronota::_to_json(json& json_ouput) const {
    if (!average_data.area.empty()) {
        json_ouput = {{"⟨SASA⟩", average_data.area.avg()},
                      {"⟨SASA²⟩-⟨SASA⟩²", average_data.area_squared.avg() - std::pow(average_data.area.avg(), 2)}};
    }
    json_ouput["radius"] = probe_radius;
}

void Voronota::_to_disk() {
    if (output_stream) {
        output_stream->flush();
    }
}

void Voronota::_sample() {
    using namespace ranges::cpp20::views;

    // Convert single `Particle` to Voronota's `SimpleSphere`
    auto to_sphere = [&](const Particle& p) -> voronotalt::SimpleSphere {
        return {{p.pos.x(), p.pos.y(), p.pos.z()}, 0.5 * p.traits().sigma + probe_radius};
    };
    const auto spheres = spc.activeParticles() | transform(to_sphere) | ranges::to_vector;

    voronotalt::RadicalTessellation::Result result;
    voronotalt::RadicalTessellation::construct_full_tessellation(spheres, result);

    auto areas = result.cells_summaries | transform([](auto& s) { return s.sas_area; });
    const auto total_area = std::accumulate(areas.begin(), areas.end(), 0.0);
    average_data.area.add(total_area);
    average_data.area_squared.add(total_area * total_area);

    if (output_stream) {
        *output_stream << this->getNumberOfSteps() << " " << total_area << "\n";
    }
}

} // namespace Faunus::analysis