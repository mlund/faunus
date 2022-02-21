#include "analysis.h"
#include "move.h"
#include "energy.h"
#include "penalty.h"
#include "reactioncoordinate.h"
#include "multipole.h"
#include "potentials.h"
#include "aux/iteratorsupport.h"
#include "aux/eigensupport.h"
#include "aux/arange.h"
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/cache1.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <zstr.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>

#include <iomanip>
#include <iostream>
#include <memory>

namespace Faunus::Analysis {

void to_json(json& j, const Analysisbase& base) { base.to_json(j); }

/**
 * This is virtual and can be overridden in derived classes
 */
void Analysisbase::_to_disk() {}

/**
 * Wrapper function used for non-intrusive future padding around _to_disk().
 */
void Analysisbase::to_disk() { _to_disk(); }

/**
 * For each call:
 *
 * 1. increment step count
 * 2. check if interval and skipped steps matches
 * 3. if so, call `_sample()` and increment number of samples
 *
 * The call to the sampling function is timed.
 */
void Analysisbase::sample() {
    try {
        number_of_steps++;
        if (sample_interval > 0 && number_of_steps > number_of_skipped_steps) {
            if ((number_of_steps % sample_interval) == 0) {
                number_of_samples++;
                timer.start();
                _sample();
                timer.stop();
            }
        }
    } catch (std::exception& e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}

void Analysisbase::from_json(const json& j) {
    try {
        number_of_skipped_steps = j.value("nskip", 0);
        sample_interval = j.value("nstep", 0);
        _from_json(j);
    } catch (std::exception& e) {
        throw ConfigurationError("{}: {}", name, e.what());
    }
}

void Analysisbase::to_json(json& json_output) const {
    try {
        auto& j = json_output[name];
        _to_json(j); // fill in info from derived classes
        if (number_of_samples > 0) {
            j["nstep"] = sample_interval;
            j["samples"] = number_of_samples;
            if (number_of_skipped_steps > 0) {
                j["nskip"] = number_of_skipped_steps;
            }
            if (timer.result() > 0.01) { // only print if more than 1% of the time
                j["relative time"] = roundValue(timer.result());
            }
        }
        if (!cite.empty()) {
            j["reference"] = cite;
        }
    } catch (std::exception& e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}

void Analysisbase::_to_json(json&) const {}

void Analysisbase::_from_json(const json&) {}

int Analysisbase::getNumberOfSteps() const { return number_of_steps; }

Analysisbase::Analysisbase(const Space& spc, std::string_view name) : spc(spc), name(name) { assert(!name.empty()); }

Analysisbase::Analysisbase(const Space& spc, std::string_view name, int sample_interval, int number_of_skipped_steps)
    : number_of_skipped_steps(number_of_skipped_steps)
    , spc(spc)
    , sample_interval(sample_interval)
    , name(name) {}

/**
 * @brief Factory function for generating analysis based on name
 * @param name name of analysis to create
 * @param j configuration for analysis
 * @param spc space to operate on
 * @param pot Hamiltonian
 * @return shared pointer to created analysis base class
 */
std::unique_ptr<Analysisbase> createAnalysis(const std::string& name, const json& j, Space& spc,
                                             Energy::Hamiltonian& pot) {
    try {
        if (name == "atomprofile") {
            return std::make_unique<AtomProfile>(j, spc);
        } else if (name == "atomrdf") {
            return std::make_unique<AtomRDF>(j, spc);
        } else if (name == "atomdipdipcorr") {
            return std::make_unique<AtomDipDipCorr>(j, spc);
        } else if (name == "density") {
            faunus_logger->warn("`density` is replaced by `molecule_density` and `atom_density`. "
                                "Activating the latter...");
            return std::make_unique<AtomDensity>(j, spc);
        } else if (name == "molecule_density") {
            return std::make_unique<MoleculeDensity>(j, spc);
        } else if (name == "atom_density") {
            return std::make_unique<AtomDensity>(j, spc);
        } else if (name == "electricpotential") {
            return std::make_unique<ElectricPotential>(j, spc);
        } else if (name == "chargefluctuations") {
            return std::make_unique<ChargeFluctuations>(j, spc);
        } else if (name == "molrdf") {
            return std::make_unique<MoleculeRDF>(j, spc);
        } else if (name == "multipole") {
            return std::make_unique<Multipole>(j, spc);
        } else if (name == "atominertia") {
            return std::make_unique<AtomInertia>(j, spc);
        } else if (name == "inertia") {
            return std::make_unique<InertiaTensor>(j, spc);
        } else if (name == "moleculeconformation") {
            return std::make_unique<MolecularConformationID>(j, spc);
        } else if (name == "multipolemoments") {
            return std::make_unique<MultipoleMoments>(j, spc);
        } else if (name == "multipoledist") {
            return std::make_unique<MultipoleDistribution>(j, spc);
        } else if (name == "polymershape") {
            return std::make_unique<PolymerShape>(j, spc);
        } else if (name == "qrfile") {
            return std::make_unique<QRtraj>(j, spc);
        } else if (name == "psctraj") {
            return std::make_unique<PatchySpheroCylinderTrajectory>(j, spc);
        } else if (name == "reactioncoordinate") {
            return std::make_unique<FileReactionCoordinate>(j, spc);
        } else if (name == "sanity") {
            return std::make_unique<SanityCheck>(j, spc);
        } else if (name == "savestate") {
            return std::make_unique<SaveState>(j, spc);
        } else if (name == "penaltyfunction") {
            return std::make_unique<SavePenaltyEnergy>(j, spc, pot);
        } else if (name == "scatter") {
            return std::make_unique<ScatteringFunction>(j, spc);
        } else if (name == "sliceddensity") {
            return std::make_unique<SlicedDensity>(j, spc);
        } else if (name == "systemenergy") {
            return std::make_unique<SystemEnergy>(j, spc, pot);
        } else if (name == "virtualvolume") {
            return std::make_unique<VirtualVolumeMove>(j, spc, pot);
        } else if (name == "virtualtranslate") {
            return std::make_unique<VirtualTranslate>(j, spc, pot);
        } else if (name == "widom") {
            return std::make_unique<WidomInsertion>(j, spc, pot);
        } else if (name == "xtcfile") {
            return std::make_unique<XTCtraj>(j, spc);
        } else if (name == "spacetraj") {
            return std::make_unique<SpaceTrajectory>(j, spc);
        }
        // append more analysis here...
        throw ConfigurationError("unknown analysis");
    } catch (std::exception& e) {
        usageTip.pick(name);
        throw ConfigurationError("{}: {}", name, e.what());
    }
}

void CombinedAnalysis::sample() {
    for (auto& analysis : this->vec) {
        analysis->sample();
    }
}

void CombinedAnalysis::to_disk() {
    for (auto& analysis : this->vec) {
        analysis->to_disk();
    }
}

CombinedAnalysis::CombinedAnalysis(const json& json_array, Space& spc, Energy::Hamiltonian& pot) {
    if (!json_array.is_array()) {
        throw ConfigurationError("json array expected");
    }
    for (const auto& j : json_array) {
        try {
            const auto& [key, json_parameters] = jsonSingleItem(j);
            vec.emplace_back(createAnalysis(key, json_parameters, spc, pot));
        } catch (std::exception& e) {
            throw ConfigurationError("analysis: {}", e.what()).attachJson(j);
        }
    }
}

void SystemEnergy::normalize() {
    const auto sum = energy_histogram.sumy();
    for (auto& i : energy_histogram.getMap()) {
        i.second = i.second / sum;
    }
}

void SystemEnergy::_sample() {
    const auto energies = calculateEnergies(); // current energy from all terms in Hamiltonian
    const auto total_energy = ranges::accumulate(energies, 0.0);
    if (std::isfinite(total_energy)) {
        mean_energy += total_energy;
        mean_squared_energy += total_energy * total_energy;
    }
    *output_stream << fmt::format("{:10d}{}{:.6E}", getNumberOfSteps(), separator, total_energy);
    for (auto energy : energies) {
        *output_stream << fmt::format("{}{:.6E}", separator, energy);
    }
    *output_stream << "\n";
    // ehist(tot)++;
}

void SystemEnergy::_to_json(json& j) const {
    j = {{"file", file_name}, {"init", initial_energy}, {"final", calculateEnergies()}};
    if (!mean_energy.empty()) {
        j["mean"] = mean_energy.avg();
        j["Cv/kB"] = mean_squared_energy.avg() - std::pow(mean_energy.avg(), 2);
    }
    roundJSON(j, 5);
    // normalize();
    // ehist.save( "distofstates.dat" );
}

void SystemEnergy::_from_json(const json& j) { file_name = MPI::prefix + j.at("file").get<std::string>(); }
void SystemEnergy::createOutputStream() {
    output_stream = IO::openCompressedOutputStream(file_name, true);
    if (auto suffix = file_name.substr(file_name.find_last_of('.') + 1); suffix == "csv") {
        separator = ",";
    } else {
        separator = " ";
        *output_stream << "#";
    }
    *output_stream << fmt::format("{:>9}{}{:12}", "step", separator, "total");
    ranges::for_each(hamiltonian,
                     [&](auto& energy) { *output_stream << fmt::format("{}{:12}", separator, energy->name); });
    *output_stream << "\n";
}

SystemEnergy::SystemEnergy(const json& j, const Space& spc, Energy::Hamiltonian& hamiltonian)
    : Analysisbase(spc, "systemenergy"), hamiltonian(hamiltonian) {
    from_json(j);
    energy_histogram.setResolution(0.25);
    initial_energy = ranges::accumulate(calculateEnergies(), 0.0);
    createOutputStream();
}
std::vector<double> SystemEnergy::calculateEnergies() const {
    Change change;
    change.everything = true; // trigger full energy calculation
    return hamiltonian | ranges::views::transform([&](auto& i) { return i->energy(change); }) | ranges::to_vector;
}

void SystemEnergy::_to_disk() {
    output_stream->flush(); // empty buffer
}

// --------------------------------

void SaveState::_to_json(json& j) const { j["file"] = filename; }

void SaveState::_sample() {
    assert(sample_interval >= 0);
    if (use_numbered_files) { // tag filename with step number
        auto numbered_filename = filename;
        numbered_filename.insert(filename.find_last_of('.'), fmt::format("_{}", getNumberOfSteps()));
        writeFunc(numbered_filename);
    } else {
        writeFunc(filename);
    }
}

SaveState::~SaveState() {
    if (sample_interval == -1) { // writes data just before destruction
        writeFunc(filename);
    }
}

SaveState::SaveState(json j, const Space& spc) : Analysisbase(spc, "savestate") {
    if (j.count("nstep") == 0) { // by default, disable _sample() and
        j["nstep"] = -1;         // store only when _to_disk() is called
    }
    from_json(j);
    save_random_number_generator_state = j.value("saverandom", false);
    filename = MPI::prefix + j.at("file").get<std::string>();
    use_numbered_files = !j.value("overwrite", false);
    convert_hexagonal_prism_to_cuboid = j.value("convert_hexagon", false);

    setWriteFunction(spc);
}
void SaveState::setWriteFunction(const Space& spc) {
    const auto suffix = filename.substr(filename.find_last_of('.') + 1);
    if (std::shared_ptr<StructureFileWriter> writer = createStructureFileWriter(suffix)) {
        writeFunc = [&, w = writer](auto& file) {
            if (convert_hexagonal_prism_to_cuboid) {
                saveAsCuboid(file, spc, *w);
            } else {
                w->save(file, spc.groups, spc.geometry.getLength());
            }
        };
    } else if (suffix == "json") { // JSON state file
        writeFunc = [&](auto& file) { saveJsonStateFile(file, spc); };
    } else if (suffix == "ubj") { // Universal Binary JSON state file
        writeFunc = [&](auto& file) { saveBinaryJsonStateFile(file, spc); };
    } else {
        throw ConfigurationError("unknown file extension for '{}'", filename);
    }
}

void SaveState::saveBinaryJsonStateFile(const std::string& filename, const Space& spc) const {
    if (std::ofstream f(filename, std::ios::binary); f) {
        json j;
        Faunus::to_json(j, spc);
        if (save_random_number_generator_state) {
            j["random-move"] = Move::MoveBase::slump;
            j["random-global"] = random;
        }
        auto buffer = json::to_ubjson(j); // json --> binary
        f.write((const char*)buffer.data(), buffer.size() * sizeof(decltype(buffer)::value_type));
    }
}
void SaveState::saveJsonStateFile(const std::string& filename, const Space& spc) const {
    if (std::ofstream f(filename); f) {
        json j;
        Faunus::to_json(j, spc);
        j.erase("reactionlist");
        if (save_random_number_generator_state) {
            j["random-move"] = Move::MoveBase::slump;
            j["random-global"] = random;
        }
        f << std::setw(1) << j;
    }
}

void SaveState::saveAsCuboid(const std::string& filename, const Space& spc, StructureFileWriter& writer) const {
    auto hexagonal_prism = std::dynamic_pointer_cast<Geometry::HexagonalPrism>(spc.geometry.asSimpleGeometry());
    if (hexagonal_prism) {
        faunus_logger->debug("creating cuboid from hexagonal prism");
        const auto& [cuboid, particles] = Geometry::hexagonalPrismToCuboid(*hexagonal_prism, spc.activeParticles());
        writer.save(filename, particles.begin(), particles.end(), cuboid.getLength());
    } else {
        throw std::runtime_error("hexagonal prism required for `convert_to_hexagon`");
    }
}

PairFunctionBase::PairFunctionBase(const Space& spc, const json& j, const std::string_view name) : Analysisbase(spc, name) {
    from_json(j);
}

void PairFunctionBase::_to_json(json& j) const {
    j = {{"dr", dr / 1.0_angstrom}, {"name1", name1},       {"name2", name2},        {"file", file},
         {"dim", dimensions},       {"slicedir", slicedir}, {"thickness", thickness}};
}

void PairFunctionBase::_from_json(const json& j) {
    file = j.at("file").get<std::string>();
    name1 = j.at("name1").get<std::string>();
    name2 = j.at("name2").get<std::string>();
    dimensions = j.value("dim", 3);
    dr = j.value("dr", 0.1) * 1.0_angstrom;
    slicedir = j.value("slicedir", slicedir);
    thickness = j.value("thickness", 0);
    histogram.setResolution(dr, 0);
}
void PairFunctionBase::_to_disk() {
    if (std::ofstream f(MPI::prefix + file); f) {
        histogram.stream_decorator = [&](std::ostream& o, double r, double N) {
            const auto volume_at_r = volumeElement(r);
            if (volume_at_r > 0.0) {
                const auto total_number_of_samples = histogram.sumy();
                o << fmt::format("{:.6E} {:.6E}\n", r, N * mean_volume.avg() / (volume_at_r * total_number_of_samples));
            }
        };
        f << histogram;
    }
}
double PairFunctionBase::volumeElement(double r) const {
    switch (dimensions) {
    case 3:
        return 4.0 * pc::pi * r * r * dr;
    case 2:
        if (auto hypersphere = std::dynamic_pointer_cast<Geometry::Hypersphere2d>(spc.geometry.asSimpleGeometry())) {
            faunus_logger->trace("{}: hypersphere detected; radius set", name);
            const auto radius = hypersphere->getRadius();
            return 2.0 * pc::pi * radius * std::sin(r / radius) * dr;
        }
        return 2.0 * pc::pi * r * dr;
    case 1:
        return dr;
    default:
        throw std::runtime_error("only 1-3 dimensions supported");
    }
}

PairAngleFunctionBase::PairAngleFunctionBase(const Space& spc, const json& j, const std::string& name)
    : PairFunctionBase(spc, j, name) {
    from_json(j);
    correlation_filename = MPI::prefix + file;
    average_correlation_vs_distance.stream_decorator = [&](auto& stream, double distance, Average<double> N) {
        if (!N.empty() && std::isfinite(N.avg())) {
            stream << distance << " " << N.avg() << "\n";
        }
    };
    // rename `file` so that it will not be overwritten by base class destructor(!)
    // @todo Call The Wolf!
    file = fmt::format("{}{}.dummy", MPI::prefix, file);
}

void PairAngleFunctionBase::_to_disk() {
    if (auto f = std::ofstream(correlation_filename)) {
        f << average_correlation_vs_distance;
    }
}

void PairAngleFunctionBase::_from_json(const json&) { average_correlation_vs_distance.setResolution(dr, 0); }

void PerturbationAnalysisBase::_to_disk() {
    if (stream) {
        stream->flush(); // empty buffer
    }
}
PerturbationAnalysisBase::PerturbationAnalysisBase(const std::string& name, Energy::Energybase& pot, Space& spc,
                                                   const std::string& filename)
    : Analysisbase(spc, name), mutable_space(spc), pot(pot), filename(filename) {
    if (!filename.empty()) {
        this->filename = MPI::prefix + filename;
        stream = IO::openCompressedOutputStream(this->filename, true); // throws if error
    }
}

bool PerturbationAnalysisBase::collectWidomAverage(const double energy_change) {
    if (-energy_change > pc::max_exp_argument) {
        faunus_logger->warn(
            "{}: skipping sample event due to too negative energy; consider decreasing the perturbation", name);
        number_of_samples--; // update_counter is incremented by sample() so we need to decrease
        return false;
    }
    mean_exponentiated_energy_change += std::exp(-energy_change);
    return true;
}
double PerturbationAnalysisBase::meanFreeEnergy() const { return -std::log(mean_exponentiated_energy_change.avg()); }

void VirtualVolumeMove::_sample() {
    if (std::fabs(volume_displacement) <= pc::epsilon_dbl) {
        return;
    }
    const auto old_volume = mutable_space.geometry.getVolume(); // store old volume
    const auto old_energy = pot.energy(change);                 // ...and energy
    const auto scale = mutable_space.scaleVolume(old_volume + volume_displacement,
                                                 volume_scaling_method); // scale entire system to new volume
    const auto new_energy = pot.energy(change);                          // energy after scaling
    mutable_space.scaleVolume(old_volume, volume_scaling_method);        // restore saved system

    const auto energy_change = new_energy - old_energy; // system energy change
    if (collectWidomAverage(energy_change)) {
        writeToFileStream(scale, energy_change);
        sanityCheck(old_energy);
    }
}

/**
 * Check if volume and particle positions are properly restored.
 * Expensive and one would normally not perform this test and we trigger it
 * only when using log-level "debug" or lower
 */
void VirtualVolumeMove::sanityCheck(const double old_energy) {
    if (faunus_logger->level() <= spdlog::level::debug and old_energy != 0.0) {
        const auto should_be_small = 1.0 - pot.energy(change) / old_energy;
        if (std::fabs(should_be_small) > 1e-6) {
            faunus_logger->error("{} failed to restore system", name);
        }
    }
}
void VirtualVolumeMove::writeToFileStream(const Point& scale, const double energy_change) const {
    if (stream) {
        const auto mean_excess_pressure = -meanFreeEnergy() / volume_displacement; // units of kT/Å³
        *stream << fmt::format("{:d} {:.3E} {:.6E} {:.6E} {:.6E}", getNumberOfSteps(), volume_displacement,
                               energy_change, std::exp(-energy_change), mean_excess_pressure);

        // if anisotropic scaling, add an extra column with area or length perturbation
        const auto box_length = spc.geometry.getLength();
        if (volume_scaling_method == Geometry::VolumeMethod::XY) {
            const auto area_change = box_length.x() * box_length.y() * (scale.x() * scale.y() - 1.0);
            *stream << fmt::format(" {:.6E}", area_change);
        } else if (volume_scaling_method == Geometry::VolumeMethod::Z) {
            const auto length_change = box_length.z() * (scale.z() - 1.0);
            *stream << fmt::format(" {:.6E}", length_change);
        }
        *stream << "\n"; // trailing newline
    }
}

void VirtualVolumeMove::_from_json(const json& j) {
    volume_displacement = j.at("dV").get<double>();
    volume_scaling_method = j.value("scaling", Geometry::VolumeMethod::ISOTROPIC);
    if (volume_scaling_method == Geometry::VolumeMethod::ISOCHORIC) {
        throw ConfigurationError("isochoric volume scaling not allowed");
    }
}

void VirtualVolumeMove::_to_json(json& j) const {
    if (number_of_samples > 0) {
        const auto excess_pressure = -meanFreeEnergy() / volume_displacement;
        j = {{"dV", volume_displacement},
             {"scaling", volume_scaling_method},
             {"-ln\u27e8exp(-dU)\u27e9", meanFreeEnergy()},
             {"Pex/mM", excess_pressure / 1.0_millimolar},
             {"Pex/Pa", excess_pressure / 1.0_Pa},
             {"Pex/kT/" + u8::angstrom + u8::cubed, excess_pressure}};
        roundJSON(j, 5);
    }
}

VirtualVolumeMove::VirtualVolumeMove(const json& j, Space& spc, Energy::Energybase& pot)
    : PerturbationAnalysisBase("virtualvolume", pot, spc, j.value("file", ""s)) {
    cite = "doi:10.1063/1.472721";
    from_json(j);
    change.volume_change = true;
    change.everything = true;
    if (stream) {
        *stream << "# steps dV/" + u8::angstrom + u8::cubed + " du/kT exp(-du/kT) <Pex>/kT/" + u8::angstrom + u8::cubed;
        // if non-isotropic scaling, add another column with dA or dL
        if (volume_scaling_method == Geometry::VolumeMethod::XY) {
            *stream << " dA/" + u8::angstrom + u8::squared;
        } else if (volume_scaling_method == Geometry::VolumeMethod::Z) {
            *stream << " dL/" + u8::angstrom;
        }
        *stream << "\n"; // trailing newline
    }
}

void MolecularConformationID::_sample() {
    auto molecules = spc.findMolecules(molid, Space::Selection::ACTIVE);
    for (const auto& group : molecules) {
        histogram[group.conformation_id]++;
    }
}
void MolecularConformationID::_to_json(json& j) const { j["histogram"] = histogram; }

MolecularConformationID::MolecularConformationID(const json& j, const Space& spc)
    : Analysisbase(spc, "moleculeconformation") {
    from_json(j);
    const auto molname = j.at("molecule").get<std::string>();
    molid = Faunus::findMoleculeByName(molname).id();
}

void QRtraj::_sample() { write_to_file(); }

void QRtraj::_to_json(json& j) const { j = {{"file", filename}}; }

QRtraj::QRtraj(const json& j, const Space& spc, const std::string &name) : Analysisbase(spc, name) {
    from_json(j);
    filename = MPI::prefix + j.value("file", "qrtraj.dat"s);
    stream = IO::openCompressedOutputStream(filename, true);
    write_to_file = [&groups = spc.groups, &stream = stream]() {
        for (const auto& group : groups) {
            for (auto it = group.begin(); it != group.trueend(); ++it) { // loop over *all* particles
                if (it < group.end()) {                                  // active particles...
                    *stream << fmt::format("{} {} ", it->charge, it->traits().sigma * 0.5);
                } else {               // inactive particles...
                    *stream << "0 0 "; // ... have zero charge and size
                }
            }
        }
        *stream << "\n"; // newline for every frame
    };
}

void QRtraj::_to_disk() {
    if (*stream) {
        stream->flush(); // empty buffer
    }
}

PatchySpheroCylinderTrajectory::PatchySpheroCylinderTrajectory(const json& j, const Space& spc)
    : QRtraj(j, spc, "psctraj") {
    if (!j.contains("file")) {
        throw ConfigurationError("missing filename");
    }
    write_to_file = [&]() {
        assert(stream);
        *stream << fmt::format("{}\nsweep {}; box ", spc.particles.size(), getNumberOfSteps())
                << spc.geometry.getLength().transpose() << "\n";

        for (const auto& group : spc.groups) {
            for (auto it = group.begin(); it != group.trueend(); ++it) {
                const auto particle_is_active = it < group.end();
                const auto scale = static_cast<double>(particle_is_active);
                *stream << scale * it->pos.transpose() << " " << scale * it->getExt().scdir.transpose() << " "
                        << scale * it->getExt().patchdir.transpose() << "\n";
            }
        }
    };
}

void FileReactionCoordinate::_to_json(json& j) const {
    json rcjson = static_cast<json>(*reaction_coordinate);
    if (rcjson.size() != 1) {
        throw std::runtime_error("error writing json for reaction coordinate");
    }
    j = rcjson.begin().value();
    j["type"] = rcjson.begin().key();
    j.erase("range");      // these are for penalty function
    j.erase("resolution"); // use only, so no need to show
    if (stream) {
        j["file"] = MPI::prefix + filename;
    }
    if (mean_reaction_coordinate) {
        j["average"] = mean_reaction_coordinate.avg();
    }
}

void FileReactionCoordinate::_sample() {
    const auto value = reaction_coordinate->operator()();
    mean_reaction_coordinate += value;
    if (stream) {
        (*stream) << fmt::format("{} {:.6f} {:.6f}\n", getNumberOfSteps(), value, mean_reaction_coordinate.avg());
    }
}

FileReactionCoordinate::FileReactionCoordinate(
    const Space& spc, const std::string& filename,
    std::unique_ptr<ReactionCoordinate::ReactionCoordinateBase> reaction_coordinate)
    : Analysisbase(spc, "reactioncoordinate")
    , filename(filename)
    , reaction_coordinate(std::move(reaction_coordinate)) {
    if (filename.empty()) {
        faunus_logger->warn("{}: no filename given - only the mean coordinate will be saved", name);
    } else {
        stream = IO::openCompressedOutputStream(MPI::prefix + filename, true);
    }
}

FileReactionCoordinate::FileReactionCoordinate(const json& j, const Space& spc)
    : FileReactionCoordinate(
          spc, j.value("filename", ""s),
          ReactionCoordinate::createReactionCoordinate({{j.at("type").get<std::string>(), j}}, spc)) {
    from_json(j);
}

void FileReactionCoordinate::_to_disk() {
    if (stream) {
        stream->flush(); // empty buffer
    }
}

/**
 * This searches for an inactive group of type `molid`
 * and prepares the `change` object for energy evaluation. If no
 * ghost molecule is found, `change` is left empty.
 */
void WidomInsertion::selectGhostGroup() {
    change.clear();
    auto inactive_groups = spc.findMolecules(molid, Space::Selection::INACTIVE);
    if (!ranges::cpp20::empty(inactive_groups)) {
        const auto& group = *inactive_groups.begin(); // select first group
        if (group.empty() && group.capacity() > 0) {  // must be inactive and have a non-zero capacity
            auto& group_changes = change.groups.emplace_back();
            group_changes.group_index = spc.getGroupIndex(group); // group index in space
            group_changes.all = true;
            group_changes.internal = group.isAtomic(); // internal energy only if non-molecular
        }
    }
}

void WidomInsertion::_sample() {
    selectGhostGroup(); // will prepare `change`
    if (change.empty()) {
        faunus_logger->warn("{}: no inactive {} groups available", name, Faunus::molecules[molid].name);
    } else {
        auto& group = mutable_space.groups.at(change.groups.at(0).group_index); // inactive "ghost" group
        group.resize(group.capacity());                                         // activate ghost
        ParticleVector particles;                                               // particles to insert
        for (int cnt = 0; cnt < number_of_insertions; ++cnt) {
            particles =
                inserter->operator()(mutable_space.geometry, Faunus::molecules[molid], spc.particles); // random pos&orientation
            updateGroup(group, particles);
            const auto energy_change = pot.energy(change); // in kT
            collectWidomAverage(energy_change);
        }
        group.resize(0); // de-activate ghost
    }
}

void WidomInsertion::updateGroup(Space::GroupType& group, const ParticleVector& particles) {
    assert(particles.size() == group.size());
    std::copy(particles.begin(), particles.end(), group.begin()); // copy to ghost group
    if (absolute_z_coords) {
        std::for_each(group.begin(), group.end(), [](Particle& i) { i.pos.z() = std::fabs(i.pos.z()); });
    }
    if (auto mass_center = group.massCenter()) { // update molecular mass-center for molecular groups
        (*mass_center).get() =
            Geometry::massCenter(group.begin(), group.end(), this->spc.geometry.getBoundaryFunc(), -group.begin()->pos);
    }
}

void WidomInsertion::_to_json(json& j) const {
    if (!mean_exponentiated_energy_change.empty()) {
        const double excess_chemical_potential = meanFreeEnergy();
        j = {{"molecule", Faunus::molecules[molid].name},
             {"insertions", mean_exponentiated_energy_change.size()},
             {"absz", absolute_z_coords},
             {"insertscheme", *inserter},
             {u8::mu + "/kT", {{"excess", excess_chemical_potential}}}};
    }
}

void WidomInsertion::_from_json(const json& j) {
    number_of_insertions = j.at("ninsert").get<int>();
    absolute_z_coords = j.value("absz", false);
    if (auto ptr = std::dynamic_pointer_cast<RandomInserter>(inserter); ptr) {
        ptr->dir = j.value("dir", Point({1, 1, 1}));
    } // set insert directions for RandomInserter

    const auto molecule_name = j.at("molecule").get<std::string>();
    molid = findMoleculeByName(molecule_name).id();
}

WidomInsertion::WidomInsertion(const json& j, Space& spc, Energy::Hamiltonian& pot)
    : PerturbationAnalysisBase("widom", pot, spc) {
    cite = "doi:10/dkv4s6";
    inserter = std::make_shared<RandomInserter>();
    from_json(j);
}

double DensityBase::updateVolumeStatistics() {
    const auto volume = spc.geometry.getVolume();
    mean_volume += volume;
    mean_cubic_root_of_volume += std::cbrt(volume);
    mean_inverse_volume += 1.0 / volume;
    return volume;
}

void DensityBase::_sample() {
    const auto volume = updateVolumeStatistics();
    for (auto [id, number] : count()) {
        mean_density[id] += number / volume;
        probability_density.at(id)(number)++;
    }
}

void DensityBase::_to_json(json& j) const {
    j["<V>"] = mean_volume.avg();
    j["<∛V>"] = mean_cubic_root_of_volume.avg();
    j["∛<V>"] = std::cbrt(mean_volume.avg());
    j["<1/V>"] = mean_inverse_volume.avg();

    auto& densities = j["densities"] = json::object();
    for (auto [id, density] : mean_density) {
        if (!density.empty() && density.avg() > pc::epsilon_dbl) {
            densities[std::string{names.at(id)}] = json({{"c/M", density.avg() / 1.0_molar}});
        }
    }
    roundJSON(j, 4);
}

/**
 * Write histograms to disk
 */
void DensityBase::writeTable(std::string_view name, Table& table) {
    if (table.size() <= 1) { // do not save non-existent or non-fluctuating groups
        return;
    }
    table.stream_decorator = [&table](auto& stream, int N, double counts) {
        if (counts > 0) {
            stream << fmt::format("{} {} {:.3E}\n", N, counts, counts / table.sumy());
        }
    };
    const auto filename = fmt::format("{}rho-{}.dat", MPI::prefix, name);
    if (std::ofstream file(filename); file) {
        file << "# N counts probability\n" << table;
    } else {
        throw std::runtime_error("could not writeKeyValuePairs "s + filename);
    }
}

void DensityBase::_to_disk() {
    for (auto [id, table] : probability_density) { // atomic molecules
        writeTable(names.at(id), table);
    }
}

void AtomDensity::_sample() {
    DensityBase::_sample();
    std::set<id_type> unique_reactive_atoms;
    for (const auto& reaction : reactions) { // in case of reactions involving atoms (swap moves)
        const auto& atomic_ids = reaction.getReactantsAndProducts().first;
        unique_reactive_atoms.insert(atomic_ids.begin(), atomic_ids.end());
    }
    ranges::cpp20::for_each(unique_reactive_atoms, [&](auto id) {
        const auto count = spc.countAtoms(id);
        atomswap_probability_density[id](count)++;
    });
}

/**
 * @brief Counts atoms in atomic groups
 */
std::map<DensityBase::id_type, int> AtomDensity::count() const {
    namespace rv = ranges::cpp20::views;

    // All ids incl. inactive are counted; std::vector ensures constant lookup (index = id)
    std::vector<int> atom_count(names.size(), 0);

    // Count number of active atoms in atomic groups
    auto particle_ids_in_atomic_groups =
        spc.groups | rv::filter(&Group::isAtomic) | rv::join | rv::transform(&Particle::id);
    ranges::cpp20::for_each(particle_ids_in_atomic_groups, [&](auto id) { atom_count.at(id)++; });

    // Copy vector --> map
    id_type id = 0U;
    std::map<id_type, int> map;
    ranges::cpp20::for_each(atom_count, [&id, &map](auto count) { map.emplace_hint(map.end(), id++, count); });
    return map;
}

AtomDensity::AtomDensity(const json& j, Space& spc) : DensityBase(spc, Faunus::atoms, "atom_density") {
    from_json(j);
    for (const auto& reaction : Faunus::reactions) { // in case of reactions involving atoms (swap moves)
        const auto reactive_atomic_species = reaction.getReactantsAndProducts().first;
        for (auto atomid : reactive_atomic_species) {
            atomswap_probability_density[atomid].setResolution(1, 0);
        }
    }
}

void AtomDensity::_to_disk() {
    DensityBase::_to_disk();
    for (const auto& reaction : Faunus::reactions) {
        const auto reactive_atomic_species = reaction.getReactantsAndProducts().first;
        for (auto atomid : reactive_atomic_species) {
            writeTable(names.at(atomid), atomswap_probability_density.at(atomid));
        }
    }
}

/**
 * @brief Counts active, molecular groups
 */
std::map<DensityBase::id_type, int> MoleculeDensity::count() const {
    using namespace ranges::cpp20;
    std::map<id_type, int> molecular_group_count;

    // ensure that also inactive groups are registered (as zero)
    for_each(Faunus::molecules | views::filter(&MoleculeData::isMolecular),
             [&](auto& moldata) { molecular_group_count[moldata.id()] = 0; });

    auto non_empty_molecular = [](const Group& group) { return group.isMolecular() && !group.empty(); };
    auto molecular_group_ids = spc.groups | views::filter(non_empty_molecular) | views::transform(&Group::id);

    for_each(molecular_group_ids, [&](auto id) { molecular_group_count[id]++; });
    return molecular_group_count;
}

MoleculeDensity::MoleculeDensity(const json& j, Space& spc) : DensityBase(spc, Faunus::molecules, "molecule_density") {
    from_json(j);
}

void SanityCheck::_sample() {
    try {
        checkGroupsCoverParticles();
        for (const auto& group : spc.groups) {
            checkWithinContainer(group);
            checkMassCenter(group);
        }
    } catch (std::exception& e) {
        PQRWriter().save(fmt::format("{}step{}-error.pqr", MPI::prefix, getNumberOfSteps()), spc.groups,
                         spc.geometry.getLength());
        throw std::runtime_error(e.what());
    }
}
void SanityCheck::checkGroupsCoverParticles() {
    size_t particle_index = 0;
    for (const auto& group : spc.groups) {
        for (auto it = group.begin(); it != group.trueend(); ++it) {
            const auto address_of_particle = std::addressof(*it);
            if (address_of_particle != std::addressof(spc.particles.at(particle_index++))) {
                throw std::runtime_error("group vector out of sync");
            }
        }
    }
    if (particle_index != spc.particles.size()) {
        throw std::runtime_error("particle <-> group mismatch");
    }
}

/**
 * Reports all particles that lies outside the simulation cell boundaries.
 *
 * @param group Group to check
 * @throw if any particle is outside simulation cell
 */
void SanityCheck::checkWithinContainer(const Space::GroupType& group) {
    bool outside_simulation_cell = false;
    auto outside_particles = group | ranges::cpp20::views::filter([&](const Particle& particle) {
                                 return spc.geometry.collision(particle.pos);
                             });

    std::for_each(outside_particles.begin(), outside_particles.end(), [&](const auto& particle) {
        outside_simulation_cell = true;
        auto group_str = fmt::format("{}{}", group.traits().name, spc.getGroupIndex(group));
        if (group.traits().numConformations() > 1) {
            group_str += fmt::format(" (conformation {})", group.conformation_id);
        }
        faunus_logger->error("step {}: atom {}{} in molecule {}", getNumberOfSteps(), particle.traits().name,
                             group.getParticleIndex(particle), group_str);
        faunus_logger->error("  (x,y,z) = {:.3f} {:.3f} {:.3f}", particle.pos.x(), particle.pos.y(), particle.pos.z());
    });

    if (outside_simulation_cell) {
        throw std::runtime_error("particle(s) outside simulation cell");
    }
}

void SanityCheck::checkMassCenter(const Space::GroupType& group) {
    if (group.isMolecular() && !group.empty()) {
        const auto mass_center =
            Geometry::massCenter(group.begin(), group.end(), spc.geometry.getBoundaryFunc(), -group.mass_center);
        const auto distance = spc.geometry.vdist(group.mass_center, mass_center).norm();
        if (distance > mass_center_tolerance) {
            throw std::runtime_error(fmt::format(
                "step {}: {}{} mass center out of sync by {:.3f} Å. This *may* be due to a "
                "molecule being longer than half the box length: consider increasing the simulation cell size.",
                getNumberOfSteps(), group.traits().name, spc.getGroupIndex(group), distance));
        }
    }
}

SanityCheck::SanityCheck(const json& j, const Space& spc) : Analysisbase(spc, "sanity") {
    from_json(j);
    sample_interval = j.value("nstep", -1);
}

void AtomRDF::sampleDistance(const Particle& particle1, const Particle& particle2) {
    const auto distance = spc.geometry.vdist(particle1.pos, particle2.pos);
    if (slicedir.sum() > 0) {
        if (distance.cwiseProduct((Point::Ones() - slicedir.cast<double>()).cwiseAbs()).norm() < thickness) {
            histogram(distance.cwiseProduct(slicedir.cast<double>()).norm())++;
        }
    } else {
        histogram(distance.norm())++;
    }
}

void AtomRDF::_sample() {
    mean_volume += spc.geometry.getVolume(dimensions);
    if (id1 != id2) {
        sampleDifferent();
    } else {
        sampleIdentical();
    }
}
void AtomRDF::sampleIdentical() {
    auto particles = spc.findAtoms(id1); // (id1 == id2)
    for (auto i = particles.begin(); i != particles.end(); ++i) {
        for (auto j = i; ++j != particles.end();) {
            sampleDistance(*i, *j);
        }
    }
}
void AtomRDF::sampleDifferent() {
    auto particles1 = spc.findAtoms(id1);
    auto particles2 = spc.findAtoms(id2);
    for (const auto& i : particles1) {
        for (const auto& j : particles2) {
            sampleDistance(i, j);
        }
    }
}
AtomRDF::AtomRDF(const json& j, const Space& spc) : PairFunctionBase(spc, j, "atomrdf") {
    id1 = Faunus::findAtomByName(name1).id();
    id2 = Faunus::findAtomByName(name2).id();
}

void MoleculeRDF::_sample() {
    mean_volume += spc.geometry.getVolume(dimensions);
    if (id1 != id2) {
        sampleDifferent();
    } else {
        sampleIdentical();
    }
}
void MoleculeRDF::sampleIdentical() {
    auto groups = spc.findMolecules(id1);
    for (auto i = groups.begin(); i != groups.end(); ++i) {
        for (auto j = i; ++j != groups.end();) {
            sampleDistance(*i, *j);
        }
    }
}
void MoleculeRDF::sampleDifferent() {
    auto ids1 = spc.findMolecules(id1);
    auto ids2 = spc.findMolecules(id2);
    for (const auto& i : ids1) {
        for (const auto& j : ids2) {
            sampleDistance(i, j);
        }
    }
}

void MoleculeRDF::sampleDistance(const Group& group_i, const Group& group_j) {
    const auto distance = sqrt(spc.geometry.sqdist(group_i.mass_center, group_j.mass_center));
    histogram(distance)++;
}

MoleculeRDF::MoleculeRDF(const json& j, const Space& spc) : PairFunctionBase(spc, j, "molrdf") {
    id1 = findMoleculeByName(name1).id();
    id2 = findMoleculeByName(name2).id();
}

void AtomDipDipCorr::_sample() {
    mean_volume += spc.geometry.getVolume(dimensions);
    const auto particles = spc.activeParticles();

    auto sample_distance = [&](const auto& particle1, const auto& particle2, const Point& distance) {
        if (particle1.hasExtension() && particle2.hasExtension()) {
            const auto cosine_angle = particle1.getExt().mu.dot(particle2.getExt().mu);
            const auto r = distance.norm();
            average_correlation_vs_distance(r) += cosine_angle;
            histogram(r)++; // get g(r) for free
        }
    };

    for (auto i = particles.begin(); i != particles.end(); ++i) {
        for (auto j = i; ++j != particles.end();) {
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                const Point distance = spc.geometry.vdist(i->pos, j->pos);
                if (slicedir.sum() > 0) {
                    if (distance.cwiseProduct(slicedir.cast<double>()).norm() < thickness) {
                        sample_distance(*i, *j, distance);
                    }
                } else {
                    sample_distance(*i, *j, distance);
                }
            }
        }
    }
}
AtomDipDipCorr::AtomDipDipCorr(const json& j, const Space& spc) : PairAngleFunctionBase(spc, j, "atomdipdipcorr") {
    id1 = findAtomByName(name1).id();
    id2 = findAtomByName(name2).id();
}

// =============== XTCtraj ===============

/**
 * @param spc Space to operate on
 * @param filename Output xtc file name
 * @param molecule_names Save only a subset of molecules with matching names (default empty = all)
 */
XTCtraj::XTCtraj(Space& spc, const std::string& filename, const std::vector<std::string>& molecule_names)
    : Analysisbase(spc, "xtcfile") {
    namespace rv = ranges::cpp20::views;
    writer = std::make_unique<XTCWriter>(filename);
    if (!molecule_names.empty()) {
        group_ids = Faunus::names2ids(Faunus::molecules, molecule_names); // molecule types to save
        group_indices = group_ids |
                        rv::transform([&](auto id) { return spc.findMolecules(id, Space::Selection::ALL); }) |
                        ranges::views::cache1 | rv::join |
                        rv::transform([&](const Group& group) { return spc.getGroupIndex(group); }) | ranges::to_vector;
        if (group_indices.empty()) {
            throw ConfigurationError("xtc selection is empty - nothing to sample");
        }
    }
}

XTCtraj::XTCtraj(const json& j, Space& spc)
    : XTCtraj(spc, MPI::prefix + j.at("file").get<std::string>(), j.value("molecules", std::vector<std::string>())) {
    Analysisbase::from_json(j);
}

void XTCtraj::_to_json(json& j) const {
    j["file"] = writer->filename;
    if (!group_ids.empty()) {
        j["molecules"] = group_ids |
                         ranges::cpp20::views::transform([](auto id) { return Faunus::molecules.at(id).name; }) |
                         ranges::to_vector;
    }
}

/**
 * Write one frame to the xtc file containing both active AND inactive particles.
 */
void XTCtraj::_sample() {
    namespace rv = ranges::cpp20::views;
    if (group_ids.empty()) {
        auto positions = spc.particles | rv::transform(&Particle::pos);
        writer->writeNext(spc.geometry.getLength(), positions.begin(), positions.end());
    } else {
        auto positions = group_indices | rv::transform([&](auto i) { return spc.groups.at(i).all(); }) |
                         ranges::views::cache1 | rv::join | rv::transform(&Particle::pos);
        writer->writeNext(spc.geometry.getLength(), positions.begin(), positions.end());
    }
}

// =============== MultipoleDistribution ===============

double MultipoleDistribution::g2g(const Space::GroupType& group1, const Space::GroupType& group2) {
    double energy = 0.0;
    for (const auto& particle_i : group1) {
        for (const auto& particle_j : group2) {
            energy += particle_i.charge * particle_j.charge / spc.geometry.vdist(particle_i.pos, particle_j.pos).norm();
        }
    }
    return energy;
}

/**
 * @note `fmt` is currently included w. spdlog but has been accepted into c++20.
 */
void MultipoleDistribution::save() const {
    if (number_of_samples > 0) {
        if (std::ofstream stream(MPI::prefix + filename.c_str()); stream) {
            stream << "# Multipolar energies (kT/lB)\n"
                   << fmt::format("# {:>8}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n", "R", "exact", "tot", "ii",
                                  "id", "dd", "iq", "mucorr");
            for (auto [r_bin, energy] : mean_energy) {
                const auto distance = r_bin * dr;
                const auto u_tot = energy.ion_ion.avg() + energy.ion_dipole.avg() + energy.dipole_dipole.avg() +
                                   energy.ion_quadrupole.avg();
                stream << fmt::format("{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}\n", distance,
                                      energy.exact.avg(), u_tot, energy.ion_ion.avg(), energy.ion_dipole.avg(),
                                      energy.dipole_dipole.avg(), energy.ion_quadrupole.avg(),
                                      energy.dipole_dipole_correlation.avg());
            }
        }
    }
}

void MultipoleDistribution::_sample() {
    for (auto& group1 : spc.findMolecules(ids[0])) {     // find active molecules
        for (auto& group2 : spc.findMolecules(ids[1])) { // find active molecules
            if (&group1 != &group2) {
                const auto a = Faunus::toMultipole(group1, spc.geometry.getBoundaryFunc());
                const auto b = Faunus::toMultipole(group2, spc.geometry.getBoundaryFunc());
                const auto distance = spc.geometry.vdist(group1.mass_center, group2.mass_center);
                auto& data = mean_energy[to_bin(distance.norm(), dr)];
                data.exact += g2g(group1, group2);
                data.ion_ion += a.charge * b.charge / distance.norm();
                data.ion_dipole += q2mu(a.charge * b.getExt().mulen, b.getExt().mu, b.charge * a.getExt().mulen,
                                        a.getExt().mu, distance);
                data.dipole_dipole +=
                    mu2mu(a.getExt().mu, b.getExt().mu, a.getExt().mulen * b.getExt().mulen, distance);
                data.ion_quadrupole += q2quad(a.charge, b.getExt().Q, b.charge, a.getExt().Q, distance);
                data.dipole_dipole_correlation += a.getExt().mu.dot(b.getExt().mu);
            }
        }
    }
}

void MultipoleDistribution::_to_json(json& j) const { j = {{"molecules", names}, {"file", filename}, {"dr", dr}}; }

MultipoleDistribution::MultipoleDistribution(const json& j, const Space& spc) : Analysisbase(spc, "Multipole Distribution") {
    from_json(j);
    dr = j.at("dr").get<double>();
    filename = j.at("file").get<std::string>();
    names = j.at("molecules").get<decltype(names)>();
    ids = Faunus::names2ids(molecules, names);
    if (ids.size() != 2) {
        throw std::runtime_error("specify exactly two molecules");
    }
}

void MultipoleDistribution::_to_disk() { save(); }

// =============== AtomInertia ===============

void AtomInertia::_to_json(json& j) const {
    j["index"] = atom_id; // atom id
}
Point AtomInertia::compute() {
    auto slice = spc.findAtoms(atom_id);
    const auto cm = Geometry::massCenter(slice.begin(), slice.end(), spc.geometry.getBoundaryFunc());
    const auto I = Geometry::inertia(slice.begin(), slice.end(), cm, spc.geometry.getBoundaryFunc());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(I);
    return esf.eigenvalues();
}
void AtomInertia::_sample() {
    if (output_stream) {
        output_stream << getNumberOfSteps() << " " << compute().transpose() << "\n";
    }
}
AtomInertia::AtomInertia(const json& j, const Space& spc) : Analysisbase(spc, "Atomic Inertia Eigenvalues") {
    from_json(j);
    filename = MPI::prefix + j.at("file").get<std::string>();
    output_stream.open(filename);
    atom_id = j.at("index").get<size_t>();
}
void AtomInertia::_to_disk() {
    if (output_stream) {
        output_stream.flush(); // empty buffer
    }
}

// =============== InertiaTensor ===============

void InertiaTensor::_to_json(json& j) const {
    j["indexes"] = particle_range; // range of indexes within the group
    j["index"] = group_index;      // group index
    j["file"] = filename;
}

std::pair<Point, Point> InertiaTensor::compute() const {
    const auto& group = spc.groups.at(group_index);
    auto subgroup =
        ranges::make_subrange(group.begin() + particle_range.at(0), group.begin() + particle_range.at(1) + 1);

    auto I = Geometry::inertia(subgroup.begin(), subgroup.end(), group.mass_center, spc.geometry.getBoundaryFunc());

    std::ptrdiff_t index; // index of eigenvector w. minimum eigenvalue
    const auto esf = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(I);
    Point eigen_values = esf.eigenvalues();
    eigen_values.minCoeff(&index);
    Point principle_axis = esf.eigenvectors().col(index).real(); // eigenvector corresponding to the smallest eigenvalue
    return {eigen_values, principle_axis};
}

void InertiaTensor::_sample() {
    const auto [eigen_values, principle_axis] = compute();
    *stream << fmt::format("{} {} {}\n", getNumberOfSteps(), eigen_values.transpose(), principle_axis.transpose());
}

InertiaTensor::InertiaTensor(const json& j, const Space& spc)
    : Analysisbase(spc, "Inertia Tensor") {
    from_json(j);
    group_index = j.at("index").get<size_t>();
    const auto& group = spc.groups.at(group_index);
    if (!group.massCenter()) {
        throw ConfigurationError("group must have a well-defined mass center (e.g. molecular groups)");
    }
    particle_range = j.value("indexes", std::vector<size_t>({0, group.size()})); // whole molecule by default
    if (particle_range.size() != 2) {
        throw ConfigurationError("{}: either two or no indices expected", name);
    }
    filename = MPI::prefix + j.at("file").get<std::string>();
    stream = IO::openCompressedOutputStream(filename, true);
    *stream << "# step eigen_values_xyz principal_axis_xyz\n";
}

void InertiaTensor::_to_disk() {
    stream->flush(); // empty buffer
}

// =============== MultipoleMoments ===============

void MultipoleMoments::_to_json(json& j) const {
    const auto& group = spc.groups.at(group_index);
    const auto particle1 = group.begin() + particle_range[0];
    const auto particle2 = group.begin() + particle_range[1];
    j["particles"] = fmt::format("{}{} {}{}", particle1->traits().name, particle_range[0], particle2->traits().name,
                                 particle_range[1]);
    j["molecule"] = group.traits().name;
}
MultipoleMoments::Data MultipoleMoments::calculateMultipoleMoment() const {
    const auto& group = spc.groups.at(group_index);
    Space::GroupType subgroup(group.id, group.begin() + particle_range[0], group.begin() + particle_range[1] + 1);
    const auto mass_center = use_molecular_mass_center ? group.mass_center
                                                       : Geometry::massCenter(subgroup.begin(), subgroup.end(),
                                                                              spc.geometry.getBoundaryFunc());

    MultipoleMoments::Data multipole;
    Tensor quadrupole; // quadrupole tensor
    quadrupole.setZero();
    for (const auto& particle : subgroup) {
        Point position = particle.pos - mass_center;
        spc.geometry.boundary(position);
        multipole.charge += particle.charge;
        multipole.dipole_moment += particle.charge * position;
        quadrupole += particle.charge * (3.0 * position * position.transpose() -
                                         Eigen::Matrix<double, 3, 3>::Identity() * position.squaredNorm());
    }
    quadrupole = 0.5 * quadrupole;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(quadrupole);
    multipole.eivals = esf.eigenvalues();
    std::ptrdiff_t i_eival;
    multipole.eivals.minCoeff(&i_eival);
    multipole.eivec = esf.eigenvectors().col(i_eival).real(); // eigenvector corresponding to the smallest eigenvalue
    multipole.center = mass_center;
    return multipole;
}
void MultipoleMoments::_sample() {
    if (output_stream) {
        const auto multipole = calculateMultipoleMoment();
        output_stream << getNumberOfSteps() << " " << multipole.charge << " " << multipole.dipole_moment.transpose()
                      << " " << multipole.center.transpose() << " " << multipole.eivals.transpose() << " "
                      << multipole.eivec.transpose() << "\n";
    }
}
MultipoleMoments::MultipoleMoments(const json& j, const Space& spc) : Analysisbase(spc, "Multipole Moments") {
    from_json(j);
    try {
        use_molecular_mass_center = j.value("mol_cm", true); // use the mass center of the whole molecule
        filename = MPI::prefix + j.at("file").get<std::string>();
        output_stream.open(filename); // output file

        group_index = j.at("index").get<decltype(group_index)>();
        const auto& group = spc.groups.at(group_index);

        particle_range = j.value("indexes", decltype(particle_range)({0, group.size() - 1}));

        if (particle_range.size() != 2) {
            throw ConfigurationError("exactly two `indexes` must be given");
        }
        if (particle_range[0] > particle_range[1] || particle_range[1] >= group.size()) {
            throw ConfigurationError("`indexes [{},{}]` outside {}", particle_range[0], particle_range[1],
                                     group.traits().name);
        }
    } catch (std::exception& e) {
        throw ConfigurationError("{}: {}", name, e.what());
    }
}
void MultipoleMoments::_to_disk() {
    if (output_stream) {
        output_stream.flush(); // empty buffer
    }
}

// =============== PolymerShape ===============

void PolymerShape::_to_json(json& j) const {
    if (!data.gyration_radius.empty()) {
        j = {{"molecule", Faunus::molecules[molid].name},
             {"⟨s²⟩-⟨s⟩²", data.gyration_radius_squared.avg() - std::pow(data.gyration_radius.avg(), 2)},
             {"⟨r²⟩∕⟨s²⟩", data.end_to_end_squared.avg() / data.gyration_radius_squared.avg()},
             {"Rg = √⟨s²⟩", std::sqrt(data.gyration_radius_squared.avg())},
             {"Re = √⟨r²⟩", std::sqrt(data.end_to_end_squared.avg())},
             {"asphericity (b)", data.aspherity.avg()},
             {"acylindricity (c)", data.acylindricity.avg()},
             {"relative shape anisotropy (𝜅²)", data.relative_shape_anisotropy.avg()}};
    }
}

void PolymerShape::_sample() {
    auto molecules = spc.findMolecules(molid, Space::Selection::ACTIVE);
    const auto num_molecules = std::distance(molecules.begin(), molecules.end());

    if (num_molecules > 1 && tensor_output_stream) {
        throw std::runtime_error("tensor output `file` cannot be used with multiple molecules");
    }

    for (const auto& group : molecules) {
        if (group.size() >= 2) { // two or more particles required to form a polymer
            const auto gyration_tensor =
                Geometry::gyration(group.begin(), group.end(), group.mass_center, spc.geometry.getBoundaryFunc());
            const auto principal_moment = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(gyration_tensor).eigenvalues();
            const auto gyration_radius_squared = gyration_tensor.trace();
            const auto end_to_end_squared = spc.geometry.sqdist(group.begin()->pos, std::prev(group.end())->pos);

            data.end_to_end_squared += end_to_end_squared;
            data.gyration_radius += std::sqrt(gyration_radius_squared);
            data.gyration_radius_squared += gyration_radius_squared;
            gyration_radius_histogram(std::sqrt(gyration_radius_squared))++;

            const double asphericity = 3.0 / 2.0 * principal_moment.z() - gyration_radius_squared / 2.0;
            const double acylindricity = principal_moment.y() - principal_moment.x();
            const double relative_shape_anisotropy =
                (asphericity * asphericity + 3.0 / 4.0 * acylindricity * acylindricity) /
                (gyration_radius_squared * gyration_radius_squared);

            data.aspherity += asphericity;
            data.acylindricity += acylindricity;
            data.relative_shape_anisotropy += relative_shape_anisotropy;

            if (tensor_output_stream) {
                const auto& t = gyration_tensor;
                *tensor_output_stream << fmt::format("{} {:.2f} {:5e} {:5e} {:5e} {:5e} {:5e} {:5e}\n",
                                                     this->getNumberOfSteps(), std::sqrt(gyration_radius_squared),
                                                     t(0, 0), t(0, 1), t(0, 2), t(1, 1), t(1, 2), t(2, 2));
            }
        }
    }
}

void PolymerShape::_to_disk() {
    const std::string filename = fmt::format("{}gyration_{}.dat", MPI::prefix, Faunus::molecules[molid].name);
    if (auto stream = std::ofstream(filename); stream) {
        gyration_radius_histogram.stream_decorator = [](auto& stream, auto Rg, auto observations) {
            if (observations > 0) {
                stream << fmt::format("{} {}\n", Rg, observations);
            }
        };
        stream << "# Rg N\n" << gyration_radius_histogram;
    }
    if (tensor_output_stream) {
        *tensor_output_stream << std::flush;
    }
}

PolymerShape::PolymerShape(const json& j, const Space& spc) : Analysisbase(spc, "Polymer Shape") {
    from_json(j);
    cite = "https://dx.doi.org/10/d6ff";
    if (j.count("molecules") > 0) {
        throw ConfigurationError("{}: 'molecules' is deprecated, use a single 'molecule' instead.", name);
    }
    const std::string molname = j.at("molecule");
    molid = findMoleculeByName(molname).id();
    if (Faunus::molecules[molid].atomic) {
        faunus_logger->warn("polymer shape analysis on an atomic group is unadvisable");
    }
    gyration_radius_histogram.setResolution(j.value("histogram_resolution", 0.2), 0.0);

    if (auto filename = j.value("file", ""s); !filename.empty()) {
        tensor_output_stream = IO::openCompressedOutputStream(MPI::prefix + filename, true);
        *tensor_output_stream << "# step Rg xx xy xz xy yy yz xz yz zz\n";
    }
}

void AtomProfile::_from_json(const json& j) {
    origin = j.value("origo", Point(0, 0, 0));
    dir = j.value("dir", dir);
    file = j.at("file").get<std::string>();
    names = j.at("atoms").get<decltype(names)>();                            // atom names
    const auto vec_of_ids = Faunus::names2ids(Faunus::atoms, names);         // names --> molids
    atom_id_selection = std::set<AtomData::index_type>(vec_of_ids.begin(), vec_of_ids.end()); // copy vector to set
    dr = j.value("dr", 0.1);
    table.setResolution(dr, 0);
    count_charge = j.value("charge", false);

    if (j.count("atomcom") == 1) {
        const auto atom_com = j.at("atomcom").get<std::string>();
        center_of_mass_atom_id = findAtomByName(atom_com).id();
    }
}
void AtomProfile::_to_json(json& j) const {
    j = {{"origo", origin}, {"dir", dir}, {"atoms", names}, {"file", file}, {"dr", dr}, {"charge", count_charge}};
    if (center_of_mass_atom_id >= 0) {
        j["atomcom"] = Faunus::atoms.at(center_of_mass_atom_id).name;
    }
}

void AtomProfile::_sample() {
    if (center_of_mass_atom_id >= 0) { // calc. mass center of selected atoms
        auto mass_center_particles = spc.findAtoms(center_of_mass_atom_id);
        origin = Geometry::massCenter(mass_center_particles.begin(), mass_center_particles.end(),
                                      spc.geometry.getBoundaryFunc());
    }

    auto selected_particles = spc.activeParticles() | ranges::cpp20::views::filter([&](const Particle& particle) {
                                  return atom_id_selection.count(particle.id) > 0;
                              });

    for (const auto& particle : selected_particles) {
        const auto distance = distanceToOrigin(particle.pos);
        table(distance) += count_charge ? particle.charge : 1.0;
    }
}

double AtomProfile::distanceToOrigin(const Point& position) const {
    const Point distance = spc.geometry.vdist(position, origin);
    return distance.cwiseProduct(dir.cast<double>()).norm();
}

AtomProfile::AtomProfile(const json& j, const Space& spc) : Analysisbase(spc, "atomprofile") { from_json(j); }

void AtomProfile::_to_disk() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        table.stream_decorator = [&](std::ostream& o, double r, double N) {
            double Vr = 1.0;
            auto dimensions = dir.sum();
            switch (dimensions) {
            case 3:
                Vr = 4.0 * pc::pi * std::pow(r, 2) * dr;
                break;
            case 2:
                Vr = 2.0 * pc::pi * r * dr;
                break;
            case 1:
                Vr = dr;
                break;
            default:
                throw std::runtime_error("bad dimension; sum of `dir` components must be one, two, or three");
            }
            if (Vr < dr) { // take care of the case where Vr=0
                Vr = dr;   // then the volume element is simply dr
            }

            N = N / double(number_of_samples); // average number of particles/charges
            o << fmt::format("{:.6E} {:.6E} {:.6E}\n", r, N, N / Vr * 1e27 / pc::Nav); // ... and molar concentration
        };
        f << "# r N rho/M\n" << table;
    }
}
void SlicedDensity::_from_json(const json& j) {
    file = j.at("file").get<std::string>();
    atom_names = j.at("atoms").get<decltype(atom_names)>();
    atom_ids = names2ids(atoms, atom_names);
    dz = j.value("dz", 0.1);
    if (j.count("atomcom") == 1) {
        const auto atom_com = j.at("atomcom").get<std::string>();
        center_of_mass_atom_id = Faunus::findAtomByName(atom_com).id();
    }
    histogram.setResolution(dz);
}
void SlicedDensity::_to_json(json& j) const {
    j = {{"atoms", atom_names}, {"file", file}, {"dz", dz}};
    if (center_of_mass_atom_id >= 0) {
        j["atomcom"] = Faunus::atoms.at(center_of_mass_atom_id).name;
    }
}

void SlicedDensity::_sample() {
    double z_offset = 0.0;

    if (center_of_mass_atom_id >= 0) { // calc. mass center of selected atoms
        auto mass_center_particles = spc.findAtoms(center_of_mass_atom_id);
        z_offset = Geometry::massCenter(mass_center_particles.begin(), mass_center_particles.end(),
                                        spc.geometry.getBoundaryFunc())
                       .z();
    }

    auto filtered_particles = spc.activeParticles() | ranges::cpp20::views::filter([&](const auto& particle) {
                                  return std::find(atom_ids.begin(), atom_ids.end(), particle.id) != atom_ids.end();
                              });

    for (const auto& particle : filtered_particles) {
        histogram(particle.pos.z() - z_offset)++;
    }
}
SlicedDensity::SlicedDensity(const json& j, const Space& spc) : Analysisbase(spc, "sliceddensity") { from_json(j); }

void SlicedDensity::_to_disk() {
    if (number_of_samples == 0) {
        return;
    }
    if (std::ofstream stream(MPI::prefix + file); stream) {
        const auto box = spc.geometry.getLength();
        const auto slice_volume = box.x() * box.y() * dz;
        const auto normalize = slice_volume * number_of_samples * 1.0_molar;
        const auto zhalf = 0.5 * box.z();
        stream << "# z rho/M\n";
        for (auto z : arange(-zhalf, zhalf + dz, dz)) { // interval [-Lz/2, Lz/2]
            stream << fmt::format("{:.6E} {:.6E}\n", z, histogram(z) / normalize);
        }
    }
}
void ChargeFluctuations::_sample() {
    auto filtered_molecules = spc.findMolecules(mol_iter->id(), Space::Selection::ACTIVE);
    for (const auto& group : filtered_molecules) {
        size_t particle_index = 0;
        for (const auto& particle : group) {
            atom_histograms[particle_index][particle.id]++;
            atom_mean_charges[particle_index] += particle.charge;
            particle_index++;
        }
    }
}

void ChargeFluctuations::_to_json(json& j) const {
    j["molecule"] = mol_iter->name;
    if (not filename.empty()) {
        j["pqrfile"] = filename;
    }
    if (verbose && number_of_samples > 0) {
        j = {{"dominant atoms", getPredominantParticleNames()},
             {"<q>", getMeanCharges()},
             {"std", getChargeStandardDeviation()}};
    }
}

std::vector<std::string> ChargeFluctuations::getPredominantParticleNames() const {
    auto most_frequent_name = [](const AtomHistogram& histogram) {
        const auto [atomid, count] =
            *std::max_element(histogram.begin(), histogram.end(), [](auto& a, auto& b) { return a.second < b.second; });
        return Faunus::atoms[atomid].name;
    }; // in a histogram, find the atom name with most counts

    auto atom_names = atom_histograms | ranges::cpp20::views::transform(most_frequent_name);
    return std::vector<std::string>(ranges::cpp20::begin(atom_names), ranges::cpp20::end(atom_names));
}

std::vector<double> ChargeFluctuations::getChargeStandardDeviation() const {
    auto stdev = atom_mean_charges | ranges::cpp20::views::transform([](auto& i) { return i.stdev(); });
    return std::vector<double>(ranges::cpp20::begin(stdev), ranges::cpp20::end(stdev));
}

std::vector<double> ChargeFluctuations::getMeanCharges() const {
    return std::vector<double>(ranges::cpp20::begin(atom_mean_charges), ranges::cpp20::end(atom_mean_charges));
}

void ChargeFluctuations::_to_disk() {
    if (not filename.empty()) {
        auto molecules = spc.findMolecules(mol_iter->id(), Space::Selection::ALL);
        if (not ranges::cpp20::empty(molecules)) {
            const auto particles_with_avg_charges = averageChargeParticles(*molecules.begin());
            PQRWriter().save(MPI::prefix + filename, particles_with_avg_charges.begin(),
                             particles_with_avg_charges.end(), spc.geometry.getLength());
        }
    }
}

ParticleVector ChargeFluctuations::averageChargeParticles(const Space::GroupType& group) {
    ParticleVector particles;            // temporary particle vector
    particles.reserve(group.capacity()); // allocate required memory already now
    size_t particle_index = 0;

    for (auto it = group.begin(); it < group.trueend(); ++it) {
        // we look for the id that was sampled most often
        auto id_max = std::max_element(std::begin(atom_histograms.at(particle_index)),
                                       std::end(atom_histograms.at(particle_index)),
                                       [](auto& p1, auto& p2) { return p1.second < p2.second; });
        auto& added_particle = particles.emplace_back(Faunus::atoms.at(id_max->first));
        added_particle.charge = atom_mean_charges.at(particle_index).avg();
        added_particle.pos = it->pos - group.mass_center;
        spc.geometry.boundary(added_particle.pos);
        particle_index++;
    }
    return particles;
}

/**
 * @todo replace `mol_iter` with simple molid integer
 */
ChargeFluctuations::ChargeFluctuations(const json& j, const Space& spc) : Analysisbase(spc, "chargefluctuations") {
    from_json(j);
    filename = j.value("pqrfile", ""s);
    verbose = j.value("verbose", true);

    const auto molname = j.at("molecule").get<std::string>();   // molecule name
    const auto& molecule = Faunus::findMoleculeByName(molname); // throws if not found
    if (molecule.atomic) {
        throw ConfigurationError("only molecular groups allowed");
    }

    mol_iter = Faunus::findName(Faunus::molecules, molname); // reference -> iterator
    atom_histograms.resize(molecule.atoms.size());
    atom_mean_charges.resize(molecule.atoms.size());
}
void Multipole::_sample() {
    for (const auto& group : spc.groups) {            // loop over all groups molecules
        if (group.isMolecular() and !group.empty()) { // only active, molecular groups
            const auto particle = Faunus::toMultipole(group, spc.geometry.getBoundaryFunc());
            auto& average = average_moments[group.id];
            average.charge += particle.charge;
            average.dipole_moment += particle.getExt().mulen;
            average.charge_squared += particle.charge * particle.charge;
            average.dipole_moment_squared += particle.getExt().mulen * particle.getExt().mulen;
        }
    }
}
void Multipole::_to_json(json& j) const {
    auto& molecules_json = j["molecules"];

    for (const auto& [molid, average] : average_moments) {
        const auto& molecule_name = Faunus::molecules[molid].name;
        molecules_json[molecule_name] = {{"Z", average.charge.avg()},
                                         {"Z2", average.charge_squared.avg()},
                                         {"C", average.charge_squared.avg() - std::pow(average.charge.avg(), 2)},
                                         {u8::mu, average.dipole_moment.avg()},
                                         {u8::mu + u8::squared, average.dipole_moment_squared.avg()}};
    }
}
Multipole::Multipole(const json& j, const Space& spc) : Analysisbase(spc, "multipole") { from_json(j); }

void ScatteringFunction::_sample() {
    namespace rv = ranges::cpp20::views;
    scatter_positions.clear();
    auto groups = molecule_ids | rv::transform([&](auto id) { return spc.findMolecules(id); }) | rv::join;
    ranges::cpp20::for_each(groups, [&](auto& group) {
        if (mass_center_scattering && group.isMolecular()) {
            scatter_positions.push_back(group.mass_center);
        } else {
            auto positions = group.positions();
            std::copy(positions.begin(), positions.end(), std::back_inserter(scatter_positions));
        }
    });

    // zero-padded suffix to use with `save_after_sample`
    const auto suffix = fmt::format("{:07d}", number_of_samples);

    switch (scheme) {
    case Schemes::DEBYE:
        debye->sample(scatter_positions, spc.geometry.getVolume());
        if (save_after_sample) {
            IO::writeKeyValuePairs(filename + "." + suffix, debye->getIntensity());
        }
        break;
    case Schemes::EXPLICIT_PBC:
        explicit_average_pbc->sample(scatter_positions, spc.geometry.getLength());
        if (save_after_sample) {
            IO::writeKeyValuePairs(filename + "." + suffix, explicit_average_pbc->getSampling());
        }
        break;
    case Schemes::EXPLICIT_IPBC:
        explicit_average_ipbc->sample(scatter_positions, spc.geometry.getLength());
        if (save_after_sample) {
            IO::writeKeyValuePairs(filename + "." + suffix, explicit_average_ipbc->getSampling());
        }
        break;
    }
}

void ScatteringFunction::_to_json(json& j) const {
    j = {{"molecules", molecule_names}, {"com", mass_center_scattering}};
    switch (scheme) {
    case Schemes::DEBYE:
        j["scheme"] = "debye";
        std::tie(j["qmin"], j["qmax"], std::ignore) = debye->getQMeshParameters();
        break;
    case Schemes::EXPLICIT_PBC:
        j["scheme"] = "explicit";
        j["pmax"] = explicit_average_pbc->getQMultiplier();
        j["ipbc"] = false;
        break;
    case Schemes::EXPLICIT_IPBC:
        j["scheme"] = "explicit";
        j["pmax"] = explicit_average_ipbc->getQMultiplier();
        j["ipbc"] = true;
        break;
    }
}

ScatteringFunction::ScatteringFunction(const json& j, const Space& spc) try : Analysisbase(spc, "scatter") {
    from_json(j);

    mass_center_scattering = j.value("com", true);
    save_after_sample = j.value("stepsave", false);
    filename = j.at("file").get<std::string>();
    molecule_names = j.at("molecules").get<decltype(molecule_names)>();
    molecule_ids = Faunus::names2ids(molecules, molecule_names);

    const auto cuboid = std::dynamic_pointer_cast<Geometry::Cuboid>(spc.geometry.asSimpleGeometry());

    if (const auto scheme_str = j.value("scheme", "explicit"s); scheme_str == "debye") {
        scheme = Schemes::DEBYE;
        debye = std::make_unique<Scatter::DebyeFormula<Tformfactor>>(j);
        if (cuboid) {
            faunus_logger->warn("cuboidal cell detected: consider using the `explicit` scheme");
        }
    } else if (scheme_str == "explicit") {
        if (!cuboid) {
            throw ConfigurationError("{} only valid for cuboidal cells", scheme_str);
        }
        const bool ipbc = j.value("ipbc", false);
        const int pmax = j.value("pmax", 15);
        if (ipbc) {
            scheme = Schemes::EXPLICIT_IPBC;
            explicit_average_ipbc = std::make_unique<Scatter::StructureFactorIPBC<>>(pmax);
        } else {
            scheme = Schemes::EXPLICIT_PBC;
            explicit_average_pbc = std::make_unique<Scatter::StructureFactorPBC<>>(pmax);
        }
    } else {
        throw ConfigurationError("unknown scheme");
    }
} catch (std::exception& e) {
    throw ConfigurationError("scatter: "s + e.what());
}

void ScatteringFunction::_to_disk() {
    switch (scheme) {
    case Schemes::DEBYE:
        IO::writeKeyValuePairs(filename, debye->getIntensity());
        break;
    case Schemes::EXPLICIT_PBC:
        IO::writeKeyValuePairs(filename, explicit_average_pbc->getSampling());
        break;
    case Schemes::EXPLICIT_IPBC:
        IO::writeKeyValuePairs(filename, explicit_average_ipbc->getSampling());
        break;
    }
}

void VirtualTranslate::_from_json(const json& j) {
    const auto molname = j.at("molecule").get<std::string>();
    molid = Faunus::findMoleculeByName(molname).id(); // throws if not found
    if (Faunus::molecules[molid].atomic) {
        throw ConfigurationError("atomic molecule {} not allowed", Faunus::molecules[molid].name);
    }

    perturbation_distance = j.at("dL").get<double>() * 1.0_angstrom;
    perturbation_direction = j.value("dir", Point(0.0, 0.0, 1.0));
    perturbation_direction.normalize(); // -> unit vector

    if (filename = j.value("file", ""s); !filename.empty()) {
        filename = MPI::prefix + filename;
        stream = IO::openCompressedOutputStream(filename, true); // throws if error
        *stream << "# steps dL/Å du/kT <force>/kT/Å\n"s;
    }
}
void VirtualTranslate::_sample() {
    if (std::fabs(perturbation_distance) < pc::epsilon_dbl) {
        return;
    }
    if (auto mollist = mutable_space.findMolecules(molid, Space::Selection::ACTIVE); !ranges::cpp20::empty(mollist)) {
        if (ranges::distance(mollist.begin(), mollist.end()) > 1) {
            throw std::runtime_error("exactly ONE active molecule expected");
        }
        if (auto group_it = random.sample(mollist.begin(), mollist.end()); not group_it->empty()) {
            const auto energy_change = momentarilyPerturb(*group_it);
            if (collectWidomAverage(energy_change)) { // collect average (with sanity check)
                writeToFileStream(energy_change);
            }
        }
    }
}

void VirtualTranslate::writeToFileStream(const double energy_change) const {
    if (stream) { // file to disk?
        const auto mean_force = -meanFreeEnergy() / perturbation_distance;
        *stream << fmt::format("{:d} {:.3E} {:.6E} {:.6E}\n", getNumberOfSteps(), perturbation_distance, energy_change,
                               mean_force);
    }
}

/**
 * @param group Group to temporarily displace
 * @return Energy change of perturbation (kT)
 *
 * Calculates the energy change of displacing a group and
 * then restore it to it's original position, leaving Space untouched.
 */
double VirtualTranslate::momentarilyPerturb(Space::GroupType& group) {
    change.groups.at(0).group_index = spc.getGroupIndex(group);
    const auto old_energy = pot.energy(change);
    const Point displacement_vector = perturbation_distance * perturbation_direction;
    group.translate(displacement_vector, spc.geometry.getBoundaryFunc()); // temporarily translate group
    const auto new_energy = pot.energy(change);
    group.translate(-displacement_vector, spc.geometry.getBoundaryFunc()); // restore original position
    return new_energy - old_energy;
}

void VirtualTranslate::_to_json(json& j) const {
    if (number_of_samples > 0 && std::fabs(perturbation_distance) > 0.0) {
        j = {{"dL", perturbation_distance},
             {"force", std::log(mean_exponentiated_energy_change.avg()) / perturbation_distance},
             {"dir", perturbation_direction}};
    }
}
VirtualTranslate::VirtualTranslate(const json& j, Space& spc, Energy::Energybase& pot)
    : PerturbationAnalysisBase("virtualtranslate", pot, spc, j.value("file", ""s)) {
    change.groups.resize(1);
    change.groups.front().internal = false;
    from_json(j);
}

SpaceTrajectory::SpaceTrajectory(const json& j, const Space& spc)
    : Analysisbase(spc, "space trajectory"), groups(spc.groups) {
    from_json(j);
    filename = j.at("file").get<std::string>();

    if (useCompression()) {
        stream = std::make_unique<zstr::ofstream>(MPI::prefix + filename, std::ios::binary);
    } else {
        stream = std::make_unique<std::ofstream>(MPI::prefix + filename, std::ios::binary);
    }

    if (stream != nullptr && *stream) {
        archive = std::make_unique<cereal::BinaryOutputArchive>(*stream);
    }

    if (not archive) {
        throw std::runtime_error("error creating "s + filename);
    }
}

bool SpaceTrajectory::useCompression() const {
    assert(!filename.empty());
    const auto suffix = filename.substr(filename.find_last_of('.') + 1);
    if (suffix == "ztraj") {
        return true;
    }
    if (suffix == "traj") {
        return false;
    }
    throw ConfigurationError("Trajectory file suffix must be `.traj` or `.ztraj`");
}

void SpaceTrajectory::_sample() {
    assert(archive);
    for (auto& group : groups) {
        (*archive)(group);
    }
}

void SpaceTrajectory::_to_json(json& j) const { j = {{"file", filename}}; }

void SpaceTrajectory::_to_disk() { stream->flush(); }

// -----------------------------

ElectricPotential::ElectricPotential(const json& j, const Space& spc)
    : Analysisbase(spc, "electricpotential"), potential_correlation_histogram(histogram_resolution) {
    from_json(j);
    coulomb = std::make_unique<Potential::NewCoulombGalore>();
    coulomb->from_json(j);
    getTargets(j);
    setPolicy(j);
    calculations_per_sample_event = j.value("ncalc", 1);
}

void ElectricPotential::setPolicy(const json& j) {
    output_information.clear();
    policy = j.value("policy", Policies::FIXED);
    auto stride = 0.0;
    switch (policy) {
    case Policies::FIXED:
        applyPolicy = []() {};
        break;
    case Policies::RANDOM_WALK:
        stride = j.at("stride").get<double>();
        output_information["stride"] = stride;
        applyPolicy = [&, stride] {
            auto origin = targets.begin();
            spc.geometry.randompos(origin->position, random);
            std::for_each(std::next(origin), targets.end(), [&](Target& target) {
                target.position = origin->position + stride * randomUnitVector(random);
                std::advance(origin, 1);
            });
        };
        break;
    case Policies::RANDOM_WALK_NO_OVERLAP:
        stride = j.at("stride").get<double>();
        output_information["stride"] = stride;
        applyPolicy = [&, stride] {
            auto origin = targets.begin();
            do {
                spc.geometry.randompos(origin->position, random);
            } while (overlapWithParticles(origin->position));
            std::for_each(std::next(origin), targets.end(), [&](Target& target) {
                do {
                    target.position = origin->position + stride * randomUnitVector(random);
                } while (overlapWithParticles(target.position));
                std::advance(origin, 1);
            });
        };
        break;
    default:
        throw ConfigurationError("unknown policy");
    }
}

void ElectricPotential::getTargets(const json& j) {
    if (const auto& structure = j.find("structure"); structure == j.end()) {
        throw ConfigurationError("missing structure");
    } else {
        PointVector positions;
        if (structure->is_string()) { // load positions from chemical structure file
            auto particles = loadStructure(structure->get<std::string>(), false);
            positions = particles | ranges::cpp20::views::transform(&Particle::pos) | ranges::to_vector;
        } else if (structure->is_array()) { // list of positions
            positions = structure->get<PointVector>();
        }
        std::transform(positions.begin(), positions.end(), std::back_inserter(targets), [&](auto& position) {
            Target target;
            target.position = position;
            target.potential_histogram = std::make_unique<SparseHistogram<double>>(histogram_resolution);
            return target;
        });
        if (targets.empty()) {
            throw ConfigurationError("no targets defined");
        }
    }
}

void ElectricPotential::_sample() {
    for (unsigned int i = 0; i < calculations_per_sample_event; i++) {
        applyPolicy();
        auto potential_correlation = 1.0; // phi1 * phi2 * ...
        for (auto& target : targets) {    // loop over each target point
            auto potential = calcPotentialOnTarget(target);
            target.potential_histogram->add(potential);
            target.mean_potential += potential;
            potential_correlation *= potential;
        }
        mean_potential_correlation += potential_correlation;        // <phi1 * phi2 * ...>
        potential_correlation_histogram.add(potential_correlation); // P(<phi1 * phi2 * ...>)
    }
}

double ElectricPotential::calcPotentialOnTarget(const ElectricPotential::Target& target) {
    auto potential_from_particle = [&](const Particle& particle) {
        const auto distance_to_target = std::sqrt(spc.geometry.sqdist(particle.pos, target.position));
        return coulomb->getCoulombGalore().ion_potential(particle.charge, distance_to_target);
    };
    auto potentials = spc.activeParticles() | ranges::cpp20::views::transform(potential_from_particle);
    return std::accumulate(potentials.begin(), potentials.end(), 0.0);
}

void ElectricPotential::_to_json(json& j) const {
    j = output_information;
    coulomb->to_json(j["coulomb"]);
    j["policy"] = policy;
    j["number of targets"] = targets.size();
    j["calculations per sample"] = calculations_per_sample_event;
    if (number_of_samples > 0) {
        j["correlation (βe)ⁱ⟨ϕ₁ϕ₂...ϕᵢ⟩"] = mean_potential_correlation.avg();
        auto& mean_potential_j = j["mean potentials βe⟨ϕᵢ⟩"] = json::array();
        std::transform(targets.begin(), targets.end(), std::back_inserter(mean_potential_j),
                       [](const auto& target) { return target.mean_potential.avg(); });
    }
}
void ElectricPotential::_to_disk() {
    if (auto stream = std::ofstream(MPI::prefix + "potential_correlation_histogram.dat")) {
        stream << potential_correlation_histogram;
    }
    int filenumber = 1;
    for (const auto& target : targets) {
        if (auto stream = std::ofstream(fmt::format("{}potential_histogram{}.dat", MPI::prefix, filenumber++))) {
            stream << *(target.potential_histogram);
        }
    }
}

/**
 * Checks if position lies within the spheres of diameters `sigma` defined
 * by each active particle in the system. Complexity: N
 *
 * @return True if overlap with any particle
 */
bool ElectricPotential::overlapWithParticles(const Point& position) const {
    auto overlap = [&position, &geometry = spc.geometry](const Particle& particle) {
        const auto radius = 0.5 * particle.traits().sigma;
        return geometry.sqdist(particle.pos, position) < radius * radius;
    };
    auto particles = spc.activeParticles();
    return std::any_of(particles.begin(), particles.end(), overlap);
}

SavePenaltyEnergy::SavePenaltyEnergy(const json& j, const Space& spc, const Energy::Hamiltonian& pot)
    : Analysisbase(spc, "penaltyfunction")
    , filename(j.at("file").get<std::string>())
    , penalty_energy(pot.findFirstOf<Energy::Penalty>()) {
    Analysisbase::from_json(j);
    if (!penalty_energy) {
        faunus_logger->warn("{}: analysis disabled as no penalty function found", name);
    }
}

void SavePenaltyEnergy::_sample() {
    if (penalty_energy) {
        auto name = fmt::format("{}{:06d}.{}", MPI::prefix, filenumber++, filename);
        if (auto stream = IO::openCompressedOutputStream(name, true); stream) {
            penalty_energy->streamPenaltyFunction(*stream);
        }
        name = fmt::format("{}{:06d}-histogram.{}", MPI::prefix, filenumber++, filename);
        if (auto stream = IO::openCompressedOutputStream(name, true); stream) {
            penalty_energy->streamHistogram(*stream);
        }
    }
}

} // namespace Faunus::Analysis
