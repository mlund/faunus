#include "actions.h"
#include "energy.h"
#include "potentials.h"
#include "io.h"
#include "progress_tracker.h"
#include "aux/arange.h"
#include <range/v3/view/concat.hpp>

namespace Faunus {

AngularScan::AngularScan(const json& j, const Space& spc) {
    angles = Geometry::TwobodyAnglesState(j.value("angular_resolution", 0.1));
    stream = IO::openCompressedOutputStream(j.value("file", "poses.gz"s), true);
    *stream << fmt::format("# Rotation points per molecule = {}\n", angles.quaternions_1.size())
            << "# column 0-3: Quaternion for molecule 1 (x y z w)\n"
            << "# column 4-7: Quaternion for molecule 2 (x y z w)\n"
            << "# column 8:   mass center z displacement of molecule 2\n"
            << "# column 9:   Energy (kJ/mol)\n";

    molecules.first.initialize(spc.groups, j.at("index1").get<size_t>());
    molecules.second.initialize(spc.groups, j.at("index2").get<size_t>());

    zmin = j.at("zmin").get<double>();
    zmax = j.at("zmax").get<double>();
    dz = j.at("dz").get<double>();
    max_energy = j.value("max_energy", pc::infty) * 1.0_kJmol;

    if (auto it = j.find("traj"); it != j.end()) {
        trajectory = std::make_unique<XTCWriter>(*it);
    }

    faunus_logger->debug("angular scan: COM distance range = [{:.1f}, {:.1f}) max energy = {:.1f} kJ/mol", zmin, zmax,
                         max_energy / 1.0_kJmol);
}

/**
 * @brief Set molecular index; back up original positions, centered on the origin; save reference structure to disk
 */
void AngularScan::Molecule::initialize(const Space::GroupVector& groups, int molecule_index) {
    namespace rv = ranges::cpp20::views;
    index = molecule_index;
    const auto& group = groups.at(index);
    if (group.isAtomic()) {
        throw ConfigurationError("group {} is not molecular", index);
    }
    auto as_centered_position = [&](auto& particle) -> Point { return particle.pos - group.mass_center; };
    ref_positions = group | rv::transform(as_centered_position) | ranges::to_vector;
    XYZWriter().save(fmt::format("molecule{}_reference.xyz", index), group.begin(), group.end(), Point::Zero());
}

ParticleVector AngularScan::Molecule::getRotatedReference(const Space::GroupVector& groups,
                                                          const Eigen::Quaterniond& q) {
    namespace rv = ranges::cpp20::views;
    const auto& group = groups.at(index);
    auto as_mut_position = [](auto& particle) -> Point& { return particle.pos; };
    auto particles = ParticleVector(group.begin(), group.end()); // copy particles from Space
    auto positions = ref_positions | rv::transform([&](const auto& pos) -> Point { return q * pos; });
    std::copy(positions.begin(), positions.end(), (particles | rv::transform(as_mut_position)).begin());
    return particles;
}

void AngularScan::operator()(Space& spc, Energy::Hamiltonian& hamiltonian) {
    namespace rv = ranges::cpp20::views;

    auto nonbonded = hamiltonian.findFirstOf<Energy::NonbondedBase>();
    if (!nonbonded) {
        throw std::runtime_error("No nonbonded energy!");
    }

    auto format = [](const auto& q) { return fmt::format("{:9.4f}{:9.4f}{:9.4f}{:9.4f}", q.x(), q.y(), q.z(), q.w()); };

    for (const auto z_pos : arange(zmin, zmax, dz)) {

        faunus_logger->info("angular scan: separation = {}", z_pos);
        auto translate = [=](auto& particle) { particle.pos.z() += z_pos; };
        auto progress_tracker =
            std::make_unique<ProgressIndicator::ProgressBar>(angles.quaternions_1.size() * angles.quaternions_2.size());
        assert(progress_tracker);

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
                    const auto energy = nonbonded->groupGroupEnergy(group1, group2);
#pragma omp critical
                    {
                        if (energy < max_energy) {
                            *stream << format(q1) << format(q2) << fmt::format("{:9.4f} {:9.4E}\n", z_pos, energy);
                            if (trajectory) {
                                auto positions = ranges::views::concat(group1, group2) | rv::transform(&Particle::pos);
                                trajectory->writeNext(spc.geometry.getLength(), positions.begin(), positions.end());
                            }
                        }
                    }
                } // dihedral loop
#pragma omp critical
                if (progress_tracker && (++(*progress_tracker) % 10 == 0)) {
                    progress_tracker->display();
                }
            } // second body quaternion loop
        }     // first body quaternions loop
    }         // separation loop
}

std::unique_ptr<SystemAction> createAction(std::string_view name, const json& j, Space& spc) {
    try {
        if (name == "angular_scan") {
            return std::make_unique<AngularScan>(j, spc);
        }
        throw ConfigurationError("'{}' unknown", name);
    } catch (std::exception& e) {
        throw std::runtime_error("error creating action -> "s + e.what());
    }
}

std::vector<std::unique_ptr<SystemAction>> createActionList(const json& input, Space& spc) {
    if (!input.is_array()) {
        throw ConfigurationError("json array expected");
    }

    std::vector<std::unique_ptr<SystemAction>> actions;

    for (const auto& j : input) {
        try {
            const auto& [name, params] = jsonSingleItem(j);
            actions.emplace_back(createAction(name, params, spc));
        } catch (std::exception& e) {
            throw ConfigurationError("actions: {}", e.what()).attachJson(j);
        }
    }
    return actions;
}

} // namespace Faunus
