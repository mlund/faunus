#include "analysis.h"
#include "core.h"
#include "mpicontroller.h"
#include <voronotalt/voronotalt.h>
#include <numeric>

namespace Faunus::analysis {

voronotalt::PeriodicBox get_periodic_box_from_space(const Faunus::Space& spc)
{
    const Point corner = 0.5 * spc.geometry.getLength();
    const std::vector<voronotalt::SimplePoint> box_corners = {
        voronotalt::SimplePoint(0.0 - corner.x(), 0.0 - corner.y(), 0.0 - corner.z()),
        voronotalt::SimplePoint(corner.x(), corner.y(), corner.z())};
    return voronotalt::PeriodicBox::create_periodic_box_from_corners(box_corners);
}

Voronota::Modes get_voronota_mode_from_name(const std::string& name)
{
    if (name == "full") {
        return Voronota::Modes::FULL;
    }
    else if (name == "interchain") {
        return Voronota::Modes::INTERCHAIN;
    }
    else if (name == "updateable") {
        return Voronota::Modes::UPDATEABLE;
    }
    return Voronota::Modes::INVALID;
}

std::string get_name_from_voronota_mode(const Voronota::Modes mode)
{
    if (mode == Voronota::Modes::FULL) {
        return "full"s;
    }
    else if (mode == Voronota::Modes::INTERCHAIN) {
        return "interchain"s;
    }
    else if (mode == Voronota::Modes::UPDATEABLE) {
        return "updateable"s;
    }
    return "invalid"s;
}

class Voronota::Impl
{
  public:
    Impl(const Voronota::Modes mode, const double probe_radius, const bool use_pbc)
        : mode(mode)
        , probe_radius(probe_radius)
        , use_pbc(use_pbc)
        , updateable_tessellation(false)
        , num_of_chains(0)
    {
    }

    void reset_all(const Faunus::Space& spc)
    {
        periodic_box = voronotalt::PeriodicBox();
        spheres.clear();
        grouping_of_spheres.clear();
        num_of_chains = 0;

        if (use_pbc) {
            periodic_box = get_periodic_box_from_space(spc);
        }

        for (const auto& group : spc.groups) {
            bool new_chain_active = false;
            for (const auto& particle : group) {
                if (!new_chain_active) {
                    num_of_chains++;
                    new_chain_active = true;
                }
                spheres.push_back(
                    voronotalt::SimpleSphere(particle.pos.x(), particle.pos.y(), particle.pos.z(),
                                             (0.5 * particle.traits().sigma) + probe_radius));
                if (mode == Voronota::Modes::INTERCHAIN) {
                    grouping_of_spheres.push_back(num_of_chains);
                }
            }
        }

        if (mode == Voronota::Modes::UPDATEABLE) {
            updateable_tessellation.init(spheres, periodic_box);
        }
    }

    void reset_spheres(const Faunus::Space& spc)
    {
        if ((use_pbc != periodic_box.enabled()) ||
            (use_pbc && periodic_box.enabled() &&
             !periodic_box.equals(get_periodic_box_from_space(spc)))) {
            reset_all(spc);
            return;
        }

        int current_num_of_chains = 0;
        std::size_t sphere_index = 0;
        for (const auto& group : spc.groups) {
            bool new_chain_active = false;
            for (const auto& particle : group) {
                if (!new_chain_active) {
                    current_num_of_chains++;
                    new_chain_active = true;
                }
                if ((sphere_index >= spheres.size()) ||
                    (mode == Voronota::Modes::INTERCHAIN &&
                     (sphere_index >= grouping_of_spheres.size() ||
                      grouping_of_spheres[sphere_index] != current_num_of_chains))) {
                    reset_all(spc);
                    return;
                }
                voronotalt::SimpleSphere& s = spheres[sphere_index];
                s.p.x = particle.pos.x();
                s.p.y = particle.pos.y();
                s.p.z = particle.pos.z();
                s.r = (0.5 * particle.traits().sigma) + probe_radius;
                sphere_index++;
            }
        }
    }

    Voronota::Modes mode;
    double probe_radius;
    bool use_pbc;
    voronotalt::PeriodicBox periodic_box;
    std::vector<voronotalt::SimpleSphere> spheres;
    std::vector<int> grouping_of_spheres;
    voronotalt::UpdateableRadicalTessellation updateable_tessellation;
    int num_of_chains;
};

Voronota::Voronota(Voronota::Modes mode, double probe_radius, const Faunus::Space& spc)
    : Analysis(spc, "voronota")
    , mode(mode)
    , probe_radius(probe_radius)
    , use_pbc(false)
{
    cite = "doi:10/mq8k";

    if (mode == Voronota::Modes::INVALID) {
        throw ConfigurationError("invalid running mode");
    }

    auto n_pbc = spc.geometry.asSimpleGeometry()->boundary_conditions.isPeriodic().count();
    if (n_pbc == 3) {
        use_pbc = true;
        faunus_logger->debug("{}: 3D PBC detected", name);
    }
    else if (n_pbc == 0) {
        faunus_logger->debug("{}: No PBC detected", name);
    }
    else {
        faunus_logger->warn("{}: Non-uniform PBC is currently ignored - be careful!", name);
    }
}

Voronota::Voronota(const Faunus::json& input, const Faunus::Space& spc)
    : Voronota(get_voronota_mode_from_name(input.value("mode", "full"s)),
               input.value("radius", 1.4_angstrom), spc)
{
    from_json(input);
}

Voronota::~Voronota() = default;

void Voronota::_from_json(const Faunus::json& input)
{
    if (filename = input.value("file", ""s); !filename.empty()) {
        output_stream = IO::openCompressedOutputStream(MPI::prefix + filename, true);
        *output_stream << "# step chains spheres ";
        if (mode == Voronota::Modes::FULL) {
            *output_stream << "collisions relevant_collisions area_of_contacts volume_inside_sas "
                              "area_of_sas area_of_sas_per_chain";
        }
        else if (mode == Voronota::Modes::INTERCHAIN) {
            *output_stream << "collisions relevant_collisions area_of_contacts";
        }
        else if (mode == Voronota::Modes::UPDATEABLE) {
            *output_stream
                << " area_of_contacts volume_inside_sas area_of_sas area_of_sas_per_chain";
        }
        *output_stream << "\n";
    }
}

void Voronota::_to_json(json& json_output) const
{
    json_output = {{"⟨area_of_contacts⟩", average_data.area_of_contacts.avg()},
                   {"⟨volume_inside_sas⟩", average_data.volume_inside_sas.avg()},
                   {"⟨area_of_sas⟩", average_data.area_of_sas.avg()},
                   {"⟨area_of_sas_per_chain⟩", average_data.area_of_sas_per_chain.avg()}};
    json_output["mode"] = get_name_from_voronota_mode(mode);
    json_output["radius"] = probe_radius;
}

void Voronota::_to_disk()
{
    if (output_stream) {
        output_stream->flush();
    }
}

void Voronota::_sample()
{
    if (!pimpl) {
        pimpl = std::make_unique<Impl>(mode, probe_radius, use_pbc);
        pimpl->reset_all(spc);
    }
    else {
        pimpl->reset_spheres(spc);
    }

    if (mode == Voronota::Modes::FULL) {
        voronotalt::RadicalTessellation::Result result;
        voronotalt::RadicalTessellation::construct_full_tessellation(pimpl->spheres,
                                                                     pimpl->periodic_box, result);

        const double area_of_contacts = result.total_contacts_summary.area;
        const double volume_inside_sas = result.total_cells_summary.sas_inside_volume;
        const double area_of_sas = result.total_cells_summary.sas_area;
        const double area_of_sas_per_chain =
            (area_of_sas / static_cast<double>(pimpl->num_of_chains));

        average_data.area_of_contacts.add(area_of_contacts);
        average_data.volume_inside_sas.add(volume_inside_sas);
        average_data.area_of_sas.add(area_of_sas);
        average_data.area_of_sas_per_chain.add(area_of_sas_per_chain);

        if (output_stream) {
            *output_stream << this->getNumberOfSteps() << " " << pimpl->num_of_chains << " "
                           << pimpl->spheres.size() << " " << result.total_collisions << " "
                           << result.total_relevant_collisions << " " << area_of_contacts << " "
                           << volume_inside_sas << " " << area_of_sas << " "
                           << area_of_sas_per_chain << "\n";
        }
    }
    else if (mode == Voronota::Modes::INTERCHAIN) {
        voronotalt::RadicalTessellation::Result result;
        voronotalt::RadicalTessellation::construct_full_tessellation(
            pimpl->spheres, pimpl->grouping_of_spheres, pimpl->periodic_box, result);

        const double area_of_contacts = result.total_contacts_summary.area;

        average_data.area_of_contacts.add(area_of_contacts);

        if (output_stream) {
            *output_stream << this->getNumberOfSteps() << " " << pimpl->num_of_chains << " "
                           << pimpl->spheres.size() << " " << result.total_collisions << " "
                           << result.total_relevant_collisions << " " << area_of_contacts << "\n";
        }
    }
    else if (mode == Voronota::Modes::UPDATEABLE) {
        pimpl->updateable_tessellation.update(pimpl->spheres);

        voronotalt::UpdateableRadicalTessellation::ResultSummary result =
            pimpl->updateable_tessellation.result_summary();

        const double area_of_contacts = result.total_contacts_summary.area;
        const double volume_inside_sas = result.total_cells_summary.sas_inside_volume;
        const double area_of_sas = result.total_cells_summary.sas_area;
        const double area_of_sas_per_chain =
            (area_of_sas / static_cast<double>(pimpl->num_of_chains));

        average_data.area_of_contacts.add(area_of_contacts);
        average_data.volume_inside_sas.add(volume_inside_sas);
        average_data.area_of_sas.add(area_of_sas);
        average_data.area_of_sas_per_chain.add(area_of_sas_per_chain);

        if (output_stream) {
            *output_stream << this->getNumberOfSteps() << " " << pimpl->num_of_chains << " "
                           << pimpl->spheres.size() << " " << " " << area_of_contacts << " "
                           << volume_inside_sas << " " << area_of_sas << " "
                           << area_of_sas_per_chain << "\n";
        }
    }
}

} // namespace Faunus::analysis
