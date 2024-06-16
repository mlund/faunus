#include "actions.h"
#include "energy.h"
#include "potentials.h"
#include "io.h"
#include "aux/arange.h"
#include <progress_tracker.h>
#include <range/v3/view/concat.hpp>

namespace Faunus {

AngularScan::AngularScan(const json& j, const Space& spc)
{
    const auto angular_resolution = j.at("angular_resolution").get<double>();
    angles = Geometry::TwobodyAnglesState(angular_resolution);
    max_energy = j.value("max_energy", pc::infty) * 1.0_kJmol;
    std::tie(zmin, zmax, dz) = j.at("zrange").get<std::tuple<double, double, double>>();

    const auto [ndx1, ndx2] = j.at("indices").get<std::tuple<size_t, size_t>>();
    molecules.first.initialize(spc.groups, ndx1);
    molecules.second.initialize(spc.groups, ndx2);

    if (j.contains("traj")) {
        trajectory = std::make_unique<XTCWriter>(j["traj"]);
        if (angles.size() > 1e5) {
            faunus_logger->warn("{}: large trajectory with {} frames will be generated", name, angles.size());
        }
    }

    stream = IO::openCompressedOutputStream(j.value("file", "poses.gz"s), true);
    *stream << fmt::format("# Δ⍺ = {:.1f}° -> {} x {} x {} = {} poses 💃🏽🕺🏼\n", angular_resolution / 1.0_deg,
                           angles.quaternions_1.size(), angles.quaternions_2.size(), angles.dihedrals.size(),
                           angles.size())
            << fmt::format("# Mass center separation range along z: [{}, {}) with dz = {} Å\n", zmin, zmax, dz)
            << fmt::format("# Max energy threshold: {} kJ/mol\n", max_energy / 1.0_kJmol)
            << "# Column 0-3: Quaternion for molecule 1 (x y z w)\n"
            << "# Column 4-7: Quaternion for molecule 2 (x y z w)\n"
            << "# Column 8:   Mass center z displacement of molecule 2\n"
            << "# Column 9:   Energy (kJ/mol)\n";

    faunus_logger->info("{}: COM distance range = [{:.1f}, {:.1f}) max energy = {:.1f} kJ/mol", name, zmin, zmax,
                        max_energy / 1.0_kJmol);
}

/**
 * @brief Set molecular index; back up original positions, centered on the origin; save reference structure to disk
 */
void AngularScan::Molecule::initialize(const Space::GroupVector& groups, int molecule_index)
{
    namespace rv = ranges::cpp20::views;
    index = molecule_index;
    const auto& group = groups.at(index);
    if (group.isAtomic()) {
        throw ConfigurationError("{}: group {} is not molecular", name, index);
    }
    auto as_centered_position = [&](auto& particle) -> Point { return particle.pos - group.mass_center; };
    ref_positions = group | rv::transform(as_centered_position) | ranges::to_vector;
    XYZWriter().save(fmt::format("molecule{}_reference.xyz", index), group.begin(), group.end(), Point::Zero());
}

ParticleVector AngularScan::Molecule::getRotatedReference(const Space::GroupVector& groups, const Eigen::Quaterniond& q)
{
    namespace rv = ranges::cpp20::views;
    const auto& group = groups.at(index);
    auto particles = ParticleVector(group.begin(), group.end()); // copy particles from Space
    auto positions = ref_positions | rv::transform([&](const auto& pos) -> Point { return q * pos; });
    std::copy(positions.begin(), positions.end(), (particles | rv::transform(&Particle::pos)).begin());
    return particles;
}

void AngularScan::report(const Group& group1, const Group& group2, const Eigen::Quaterniond& q1,
                         const Eigen::Quaterniond& q2, Energy::NonbondedBase& nonbonded)
{
    const auto energy = nonbonded.groupGroupEnergy(group1, group2);
    if (energy >= max_energy) {
        return;
    }

    auto format = [](const auto& q) { return fmt::format("{:8.4f}{:8.4f}{:8.4f}{:8.4f}", q.x(), q.y(), q.z(), q.w()); };

#pragma omp critical
    {
        energy_analysis.add(energy);
        *stream << format(q1) << format(q2)
                << fmt::format("{:8.4f} {:>10.3E}\n", group2.mass_center.z(), energy / 1.0_kJmol);
        if (trajectory) {
            auto positions = ranges::views::concat(group1, group2) | ranges::cpp20::views::transform(&Particle::pos);
            trajectory->writeNext({500, 500, 500}, positions.begin(), positions.end());
        }
    }
}

void AngularScan::operator()(Space& spc, Energy::Hamiltonian& hamiltonian)
{
    auto nonbonded = hamiltonian.findFirstOf<Energy::NonbondedBase>();
    if (!nonbonded) {
        throw ConfigurationError("{}: at least one nonbonded energy term required", name);
    }

    for (const auto z_pos : arange(zmin, zmax, dz)) {
        faunus_logger->info("{}: separation = {}", name, z_pos);
        auto translate = [=](auto& particle) { particle.pos.z() += z_pos; };
        auto progress_tracker =
            std::make_unique<ProgressIndicator::ProgressBar>(angles.quaternions_1.size() * angles.quaternions_2.size());
        assert(progress_tracker);
        energy_analysis.clear();

#pragma omp parallel for
        for (const auto& q1 : angles.quaternions_1) {
            auto particles1 = molecules.second.getRotatedReference(spc.groups, q1);
            auto group1 = Group(0, particles1.begin(), particles1.end());
            group1.updateMassCenter(spc.geometry.getBoundaryFunc(), {0, 0, 0});

            for (const auto& q_body2 : angles.quaternions_2) {
                for (const auto& q_dihedral : angles.dihedrals) {
                    const auto q2 = q_dihedral * q_body2; // simulataneous rotations (noncummutative)
                    auto particles2 = molecules.second.getRotatedReference(spc.groups, q2);
                    ranges::cpp20::for_each(particles2, translate);
                    auto group2 = Group(0, particles2.begin(), particles2.end());
                    group2.mass_center = {0.0, 0.0, z_pos};
                    report(group1, group2, q1, q2, *nonbonded);
                }
#pragma omp critical
                if (progress_tracker && (++(*progress_tracker) % 10 == 0)) {
                    progress_tracker->display();
                }
            }
        }
        energy_analysis.printLog();
    } // separation loop
}

std::unique_ptr<SystemAction> createAction(std::string_view name, const json& j, Space& spc)
{
    try {
        if (name == "angular_scan") {
            return std::make_unique<AngularScan>(j, spc);
        }
        throw ConfigurationError("'{}' unknown", name);
    }
    catch (std::exception& e) {
        throw std::runtime_error("error creating action -> "s + e.what());
    }
}

std::vector<std::unique_ptr<SystemAction>> createActionList(const json& input, Space& spc)
{
    if (!input.is_array()) {
        throw ConfigurationError("json array expected");
    }

    std::vector<std::unique_ptr<SystemAction>> actions;

    for (const auto& j : input) {
        try {
            const auto& [name, params] = jsonSingleItem(j);
            actions.emplace_back(createAction(name, params, spc));
        }
        catch (std::exception& e) {
            throw ConfigurationError("actions: {}", e.what()).attachJson(j);
        }
    }
    return actions;
}

void AngularScan::EnergyAnalysis::clear()
{
    mean_exp_energy.clear();
    partition_sum = 0;
    energy_sum = 0;
}

void AngularScan::EnergyAnalysis::add(const double energy)
{
    const auto exp_energy = std::exp(-energy);
    mean_exp_energy += exp_energy;
    partition_sum += exp_energy;
    energy_sum += energy * exp_energy;
}

double AngularScan::EnergyAnalysis::getFreeEnergy() const { return -std::log(mean_exp_energy.avg()); }

double AngularScan::EnergyAnalysis::getMeanEnergy() const { return energy_sum / partition_sum; }

void AngularScan::EnergyAnalysis::printLog() const
{
    faunus_logger->info("{}: free energy <w/kT> = {:.3f} mean energy <u/kT> = {:.3f}", name, getFreeEnergy(),
                        getMeanEnergy());
}

} // namespace Faunus
