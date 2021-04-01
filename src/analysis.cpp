#include "analysis.h"
#include "move.h"
#include "energy.h"
#include "reactioncoordinate.h"
#include "multipole.h"
#include "aux/iteratorsupport.h"
#include "aux/eigensupport.h"
#include <spdlog/spdlog.h>
#include <zstr.hpp>
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
                j["relative time"] = _round(timer.result());
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

Analysisbase::Analysisbase(Space& spc, const std::string& name) : spc(spc), name(name) { assert(!name.empty()); }

/**
 * @brief Factory function for generating analysis based on name
 * @param name name of analysis to create
 * @param j configuration for analysis
 * @param spc space to operate on
 * @param pot Hamiltonian
 * @return shared pointer to created analysis base class
 */
std::shared_ptr<Analysisbase> createAnalysis(const std::string& name, const json& j, Space& spc,
                                             Energy::Hamiltonian& pot) {
    try {
        if (name == "atomprofile") {
            return std::make_shared<AtomProfile>(j, spc);
        } else if (name == "atomrdf") {
            return std::make_shared<AtomRDF>(j, spc);
        } else if (name == "atomdipdipcorr") {
            return std::make_shared<AtomDipDipCorr>(j, spc);
        } else if (name == "density") {
            return std::make_shared<Density>(j, spc);
        } else if (name == "chargefluctuations") {
            return std::make_shared<ChargeFluctuations>(j, spc);
        } else if (name == "molrdf") {
            return std::make_shared<MoleculeRDF>(j, spc);
        } else if (name == "multipole") {
            return std::make_shared<Multipole>(j, spc);
        } else if (name == "atominertia") {
            return std::make_shared<AtomInertia>(j, spc);
        } else if (name == "inertia") {
            return std::make_shared<InertiaTensor>(j, spc);
        } else if (name == "moleculeconformation") {
            return std::make_shared<MolecularConformationID>(j, spc);
        } else if (name == "multipolemoments") {
            return std::make_shared<MultipoleMoments>(j, spc);
        } else if (name == "multipoledist") {
            return std::make_shared<MultipoleDistribution>(j, spc);
        } else if (name == "polymershape") {
            return std::make_shared<PolymerShape>(j, spc);
        } else if (name == "qrfile") {
            return std::make_shared<QRtraj>(j, spc);
        } else if (name == "reactioncoordinate") {
            return std::make_shared<FileReactionCoordinate>(j, spc);
        } else if (name == "sanity") {
            return std::make_shared<SanityCheck>(j, spc);
        } else if (name == "savestate") {
            return std::make_shared<SaveState>(j, spc);
        } else if (name == "scatter") {
            return std::make_shared<ScatteringFunction>(j, spc);
        } else if (name == "sliceddensity") {
            return std::make_shared<SlicedDensity>(j, spc);
        } else if (name == "systemenergy") {
            return std::make_shared<SystemEnergy>(j, spc, pot);
        } else if (name == "virtualvolume") {
            return std::make_shared<VirtualVolumeMove>(j, spc, pot);
        } else if (name == "virtualtranslate") {
            return std::make_shared<VirtualTranslate>(j, spc, pot);
        } else if (name == "widom") {
            return std::make_shared<WidomInsertion>(j, spc, pot);
        } else if (name == "xtcfile") {
            return std::make_shared<XTCtraj>(j, spc);
        } else if (name == "spacetraj") {
            return std::make_shared<SpaceTrajectory>(j, spc);
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
            const auto& [key, j_params] = jsonSingleItem(j);
            auto analysis = createAnalysis(key, j_params, spc, pot);
            vec.push_back(analysis);
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
    const auto energies = calculateAllEnergies(); // current energy from all terms in Hamiltonian
    const auto total_energy = std::accumulate(energies.begin(), energies.end(), 0.0);
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
    j = {{"file", file_name}, {"init", initial_energy}, {"final", calculateAllEnergies()}};
    if (!mean_energy.empty()) {
        j["mean"] = mean_energy.avg();
        j["Cv/kB"] = mean_squared_energy.avg() - std::pow(mean_energy.avg(), 2);
    }
    _roundjson(j, 5);
    // normalize();
    // ehist.save( "distofstates.dat" );
}

void SystemEnergy::_from_json(const json& j) {
    file_name = MPI::prefix + j.at("file").get<std::string>();
    output_stream = IO::openCompressedOutputStream(file_name, true);
    if (auto suffix = file_name.substr(file_name.find_last_of('.') + 1); suffix == "csv") {
        separator = ",";
    } else {
        separator = " ";
        *output_stream << "#";
    }
    *output_stream << fmt::format("{:>9}{}{:12}", "step", separator, "total");
    for (const auto& name : names_of_energy_terms) {
        *output_stream << fmt::format("{}{:12}", separator, name);
    }
    *output_stream << "\n";
}

SystemEnergy::SystemEnergy(const json& j, Space& spc, Energy::Hamiltonian& pot) : Analysisbase(spc, "systemenergy") {
    from_json(j);
    std::transform(pot.vec.begin(), pot.vec.end(), std::back_inserter(names_of_energy_terms),
                   [](auto& i) { return i->name; });

    calculateAllEnergies = [&pot]() {
        Change change;
        change.all = true; // trigger full energy calculation
        std::vector<double> energies;
        std::transform(pot.vec.begin(), pot.vec.end(), std::back_inserter(energies),
                       [&](auto i) { return i->energy(change); });
        return energies;
    };
    energy_histogram.setResolution(0.25);
    auto energies = calculateAllEnergies();
    initial_energy = std::accumulate(energies.begin(), energies.end(), 0.0); // initial energy
}

void SystemEnergy::_to_disk() {
    if (*output_stream) {
        output_stream->flush(); // empty buffer
    }
}

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

SaveState::SaveState(json j, Space& spc) : Analysisbase(spc, "savestate") {
    if (j.count("nstep") == 0) { // by default, disable _sample() and
        j["nstep"] = -1;         // store only when _to_disk() is called
    }
    from_json(j);

    save_random_number_generator_state = j.value("saverandom", false);
    filename = MPI::prefix + j.at("file").get<std::string>();
    use_numbered_files = !j.value("overwrite", false);
    convert_hexagonal_prism_to_cuboid = j.value("convert_hexagon", false);

    if (const auto suffix = filename.substr(filename.find_last_of('.') + 1); suffix == "aam") {
        writeFunc = [&](auto& file) { FormatAAM::save(file, spc.p); };
    } else if (suffix == "gro") {
        writeFunc = [&](auto& file) { FormatGRO::save(file, spc); };
    } else if (suffix == "pqr") {
        writeFunc = [&](auto& file) {
            if (convert_hexagonal_prism_to_cuboid) {
                auto hexagonal_prism = std::dynamic_pointer_cast<Geometry::HexagonalPrism>(spc.geo.asSimpleGeometry());
                if (hexagonal_prism) {
                    faunus_logger->debug("creating cuboidal PQR from hexagonal prism");
                    const auto& [cuboid, particles] =
                        Geometry::HexagonalPrismToCuboid(*hexagonal_prism, spc.activeParticles());
                    FormatPQR::save(file, particles, cuboid.getLength());
                } else {
                    throw std::runtime_error("hexagonal prism required for `convert_to_hexagon`");
                }
            } else {
                FormatPQR::save(file, spc.groups, spc.geo.getLength());
            }
        };
    } else if (suffix == "xyz") {
        writeFunc = [&](auto& file) { FormatXYZ::save(file, spc.p, spc.geo.getLength()); };
    } else if (suffix == "json") { // JSON state file
        writeFunc = [&](auto& file) {
            if (std::ofstream f(file); f) {
                json j;
                Faunus::to_json(j, spc);
                if (save_random_number_generator_state) {
                    j["random-move"] = Move::MoveBase::slump;
                    j["random-global"] = Faunus::random;
                }
                f << std::setw(2) << j;
            }
        };
    } else if (suffix == "ubj") { // Universal Binary JSON state file
        writeFunc = [&](auto& file) {
            if (std::ofstream f(file, std::ios::binary); f) {
                json j;
                Faunus::to_json(j, spc);
                if (save_random_number_generator_state) {
                    j["random-move"] = Move::MoveBase::slump;
                    j["random-global"] = Faunus::random;
                }
                auto v = json::to_ubjson(j); // json --> binary
                f.write((const char*)v.data(), v.size() * sizeof(decltype(v)::value_type));
            }
        };
    } else {
        throw ConfigurationError("unknown file extension for '{}'", filename);
    }
}

PairFunctionBase::PairFunctionBase(Space& spc, const json& j, const std::string& name) : Analysisbase(spc, name) {
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
                o << fmt::format("{:.6E} {:.6E}\n", r, N * mean_volume / (volume_at_r * total_number_of_samples));
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
        if (auto hypersphere = std::dynamic_pointer_cast<Geometry::Hypersphere2d>(spc.geo.asSimpleGeometry())) {
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

PairAngleFunctionBase::PairAngleFunctionBase(Space& spc, const json& j, const std::string& name)
    : PairFunctionBase(spc, j, name) {
    from_json(j);
}

void PairAngleFunctionBase::_to_disk() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        hist2.stream_decorator = [&](std::ostream& o, double r, double N) { o << r << " " << N << "\n"; };
        f << hist2;
    }
    file = file + "gofr.dat"; // name file where free g(r) is saved, and make sure that the file is not overwritten by
                              // base-destructor
}

void PairAngleFunctionBase::_from_json(const json&) { hist2.setResolution(dr, 0); }

void PerturbationAnalysisBase::_to_disk() {
    if (stream) {
        stream->flush(); // empty buffer
    }
}
PerturbationAnalysisBase::PerturbationAnalysisBase(const std::string& name, Energy::Energybase& pot, Space& spc,
                                                   const std::string& filename)
    : Analysisbase(spc, name), pot(pot), filename(filename) {
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
    if (std::fabs(volume_displacement) > 0.0) {
        const auto old_volume = spc.geo.getVolume(); // store old volume
        const auto old_energy = pot.energy(change);  // ...and energy
        const auto scale = spc.scaleVolume(old_volume + volume_displacement,
                                           volume_scaling_method); // scale entire system to new volume
        const auto new_energy = pot.energy(change);                // energy after scaling
        spc.scaleVolume(old_volume, volume_scaling_method);        // restore saved system

        const auto energy_change = new_energy - old_energy; // system energy change
        if (collectWidomAverage(energy_change)) {
            writeToFileStream(scale, energy_change);
            sanityCheck(old_energy);
        }
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
    if (stream) { // write to output file if open
        *stream << fmt::format("{:d} {:.3f} {:.10f} {:.6f} {:.6f}", getNumberOfSteps(), volume_displacement,
                               energy_change, exp(-energy_change), -meanFreeEnergy() / volume_displacement);

        // if anisotropic scaling, add an extra column with area or length perturbation
        if (volume_scaling_method == Geometry::XY) {
            const auto box_length = spc.geo.getLength();
            const auto area_change = box_length.x() * box_length.y() * (scale.x() * scale.y() - 1.0);
            *stream << fmt::format(" {:.6f}", area_change);
        } else if (volume_scaling_method == Geometry::Z) {
            const auto box_length = spc.geo.getLength();
            const auto length_change = box_length.z() * (scale.z() - 1.0);
            *stream << fmt::format(" {:.6f}", length_change);
        }
        *stream << "\n"; // trailing newline
    }
}

void VirtualVolumeMove::_from_json(const json& j) {
    volume_displacement = j.at("dV").get<double>();
    volume_scaling_method = j.value("scaling", Geometry::ISOTROPIC);
    if (volume_scaling_method == Geometry::ISOCHORIC) {
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
        _roundjson(j, 5);
    }
}

VirtualVolumeMove::VirtualVolumeMove(const json& j, Space& spc, Energy::Energybase& pot)
    : PerturbationAnalysisBase("virtualvolume", pot, spc, j.value("file", ""s)) {
    cite = "doi:10.1063/1.472721";
    from_json(j);
    change.dV = true;
    change.all = true;
    if (stream) {
        *stream << "# steps dV/" + u8::angstrom + u8::cubed + " du/kT exp(-du/kT) <Pex>/kT/" + u8::angstrom + u8::cubed;
        // if non-isotropic scaling, add another column with dA or dL
        if (volume_scaling_method == Geometry::XY) {
            *stream << " dA/" + u8::angstrom + u8::squared;
        } else if (volume_scaling_method == Geometry::Z) {
            *stream << " dL/" + u8::angstrom;
        }
        *stream << "\n"; // trailing newline
    }
}

void MolecularConformationID::_sample() {
    auto molecules = spc.findMolecules(molid, Space::ACTIVE);
    for (const auto& group : molecules) {
        histogram[group.confid]++;
    }
}
void MolecularConformationID::_to_json(json& j) const { j["histogram"] = histogram; }

MolecularConformationID::MolecularConformationID(const json& j, Space& spc)
    : Analysisbase(spc, "moleculeconformation") {
    from_json(j);
    const auto molname = j.at("molecule").get<std::string>();
    molid = Faunus::findMoleculeByName(molname).id();
}

void QRtraj::_sample() { write_to_file(); }

void QRtraj::_to_json(json& j) const { j = {{"file", filename}}; }

QRtraj::QRtraj(const json& j, Space& spc) : Analysisbase(spc, "qrfile") {
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

void FileReactionCoordinate::_to_json(json& j) const {
    json rcjson = *reaction_coordinate; // invoke to_json(...)
    if (rcjson.count(reaction_coordinate_type) == 0) {
        throw std::runtime_error("error writing json for reaction coordinate");
    }
    j = rcjson[reaction_coordinate_type];
    j["type"] = reaction_coordinate_type;
    j["file"] = filename;
    j.erase("range");      // these are for penalty function
    j.erase("resolution"); // use only, so no need to show
    if (number_of_samples > 0) {
        j["average"] = mean_reaction_coordinate.avg();
    }
}

void FileReactionCoordinate::_sample() {
    if (*stream) {
        const auto reaction_coordinate_value = reaction_coordinate->operator()();
        mean_reaction_coordinate += reaction_coordinate_value;
        (*stream) << fmt::format("{} {:.6f} {:.6f}\n", getNumberOfSteps(), reaction_coordinate_value,
                                 mean_reaction_coordinate.avg());
    }
}

FileReactionCoordinate::FileReactionCoordinate(const json& j, Space& spc) : Analysisbase(spc, "reactioncoordinate") {
    from_json(j);
    filename = MPI::prefix + j.at("file").get<std::string>();
    if (stream = IO::openCompressedOutputStream(filename); not *stream) {
        throw std::runtime_error("could not open create "s + filename);
    }
    reaction_coordinate_type = j.at("type").get<std::string>();
    reaction_coordinate = ReactionCoordinate::createReactionCoordinate({{reaction_coordinate_type, j}}, spc);
}

void FileReactionCoordinate::_to_disk() {
    if (*stream) {
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
    auto inactive_groups = spc.findMolecules(molid, Space::INACTIVE);
    if (!ranges::cpp20::empty(inactive_groups)) {
        const auto& group = *inactive_groups.begin(); // select first group
        if (group.empty() && group.capacity() > 0) {  // must be inactive and have a non-zero capacity
            auto& group_changes = change.groups.emplace_back();
            group_changes.index = spc.getGroupIndex(group); // group index in space
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
        auto& group = spc.groups.at(change.groups.at(0).index); // inactive "ghost" group
        group.resize(group.capacity());                         // activate ghost
        ParticleVector particles;                               // particles to insert
        for (int cnt = 0; cnt < number_of_insertions; ++cnt) {
            particles = inserter->operator()(spc.geo, Faunus::molecules[molid], spc.p); // random pos&orientation
            updateGroup(group, particles);
            const auto energy_change = pot.energy(change); // in kT
            collectWidomAverage(energy_change);
        }
        group.resize(0); // de-activate ghost
    }
}

void WidomInsertion::updateGroup(Space::Tgroup& group, const ParticleVector& particles) {
    assert(particles.size() == group.size());
    std::copy(particles.begin(), particles.end(), group.begin()); // copy to ghost group
    if (absolute_z_coords) {
        std::for_each(group.begin(), group.end(), [](Particle& i) { i.pos.z() = std::fabs(i.pos.z()); });
    }
    if (group.isMolecular()) { // update molecular mass-center for molecular groups
        group.cm =
            Geometry::massCenter(group.begin(), group.end(), this->spc.geo.getBoundaryFunc(), -group.begin()->pos);
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

void Density::_sample() {
    const auto volume = updateVolumeStatistics();
    auto [atom_count, molecular_group_count] = countAtomsAndMolecules();

    for (auto [atomid, number_of_atoms] : atom_count) {
        mean_atom_density[atomid] += number_of_atoms / volume;
    }
    for (auto [molid, number_of_molecules] : molecular_group_count) {
        mean_molecule_density[molid] += number_of_molecules / volume;
        molecular_group_probability_density[molid](number_of_molecules)++;
    }
    for (const auto& reaction : reactions) { // in case of reactions involving atoms (swap moves)
        const auto reactive_atomic_species = reaction.getReactantsAndProducts().first;
        for (auto atomid : reactive_atomic_species) {
            auto atomlist = spc.findAtoms(atomid);
            atomswap_probability_density[atomid](range_size(atomlist))++;
        }
    }
}

/**
 * @brief Counts all molecular groups and all atoms in atomic groups
 * @returns Pair of maps; first = atom count in atomic groups; second = count of molecular groups (keys=molid)
 */
std::pair<std::map<int, int>, std::map<int, int>> Density::countAtomsAndMolecules() {
    std::map<int, int> atom_count;
    std::map<int, int> molecular_group_count;
    // make sure all atom counts are initially zero
    for (const auto& group : spc.groups) {
        if (group.isAtomic()) {
            for (auto particle = group.begin(); particle < group.trueend(); ++particle) {
                atom_count[particle->id] = 0;
            }
        } else {
            molecular_group_count[group.id] = 0;
        }
    }

    for (const auto& group : spc.groups) {
        if (group.isAtomic()) {
            for (const auto& particle : group) {
                atom_count[particle.id]++;
            }
            atomic_group_probability_density[group.id](group.size())++;
        } else if (not group.empty()) {
            molecular_group_count[group.id]++;
        }
    }
    return {atom_count, molecular_group_count};
}

double Density::updateVolumeStatistics() {
    const auto volume = spc.geo.getVolume();
    mean_volume += volume;
    mean_cubic_root_of_volume += std::cbrt(volume);
    mean_inverse_volume += 1.0 / volume;
    return volume;
}

void Density::_to_json(json& j) const {
    j["<V>"] = mean_volume.avg();
    j["<∛V>"] = mean_cubic_root_of_volume.avg();
    j["∛<V>"] = std::cbrt(mean_volume.avg());
    j["<1/V>"] = mean_inverse_volume.avg();

    auto& j_atomic = j["atomic"] = json::object();
    auto& j_molecular = j["molecular"] = json::object();

    for (auto [atomid, density] : mean_atom_density) {
        if (!density.empty()) {
            j_atomic[Faunus::atoms[atomid].name] = json({{"c/M", density.avg() / 1.0_molar}});
        }
    }
    for (auto [molid, density] : mean_molecule_density) {
        if (!density.empty()) {
            j_molecular[Faunus::molecules[molid].name] = json({{"c/M", density.avg() / 1.0_molar}});
        }
    }
    _roundjson(j, 4);
}

Density::Density(const json& j, Space& spc) : Analysisbase(spc, "density") {
    from_json(j);
    for (const auto& molecule : Faunus::molecules) {
        if (molecule.atomic) {
            atomic_group_probability_density[molecule.id()].setResolution(1, 0);
        } else {
            molecular_group_probability_density[molecule.id()].setResolution(1, 0);
        }
    }
    for (const auto& reaction : Faunus::reactions) { // in case of reactions involving atoms (swap moves)
        const auto reactive_atomic_species = reaction.getReactantsAndProducts().first;
        for (auto atomid : reactive_atomic_species) {
            atomswap_probability_density[atomid].setResolution(1, 0);
        }
    }
}

/**
 * Write histograms to disk
 */
void Density::writeTable(const std::string& atom_or_molecule_name, Ttable& table) {
    table.stream_decorator = [&table](std::ostream& stream, int N, double counts) {
        if (counts > 0) {
            stream << fmt::format("{} {} {:.3f}\n", N, counts, counts / table.sumy());
        }
    };
    const auto filename = fmt::format("{}rho-{}.dat", MPI::prefix, atom_or_molecule_name);
    if (std::ofstream file(filename); file) {
        file << "# N counts probability\n" << table;
    } else {
        throw std::runtime_error("could not write "s + filename);
    }
}

void Density::_to_disk() {
    for (auto [molid, table] : atomic_group_probability_density) { // atomic molecules
        writeTable(Faunus::molecules[molid].name, table);
    }
    for (auto [molid, table] : molecular_group_probability_density) { // polyatomic molecules
        writeTable(Faunus::molecules[molid].name, table);
    }
    for (const auto& reaction : Faunus::reactions) {
        const auto reactive_atomic_species = reaction.getReactantsAndProducts().first;
        for (auto atomid : reactive_atomic_species) {
            writeTable(Faunus::atoms[atomid].name, atomswap_probability_density[atomid]);
        }
    }
}

void SanityCheck::_sample() {
    try {
        checkGroupsCoverParticles();
        for (const auto& group : spc.groups) {
            checkWithinContainer(group);
            checkMassCenter(group);
        }
    } catch (std::exception& e) {
        FormatPQR::save(fmt::format("{}step{}-error.pqr", MPI::prefix, getNumberOfSteps()), spc.groups,
                        spc.geo.getLength());
        throw std::runtime_error(e.what());
    }
}
void SanityCheck::checkGroupsCoverParticles() {
    size_t i = 0;
    for (const auto& group : spc.groups) {
        for (auto it = group.begin(); it != group.trueend(); ++it) {
            const auto address_of_particle = &(*it);
            if (address_of_particle != &(spc.p.at(i++))) {
                throw std::runtime_error("group vector out of sync");
            }
        }
    }
    if (i != spc.p.size()) {
        throw std::runtime_error("particle <-> group mismatch");
    }
}
void SanityCheck::checkWithinContainer(const Space::Tgroup& group) {
    for (const auto& particle : group) { // loop over active particles
        if (spc.geo.collision(particle.pos)) {
            const auto atom_index = &particle - &(*group.begin()); // yak!
            const auto group_index = spc.getGroupIndex(group);
            throw std::runtime_error(fmt::format("step {}: particle {}{} of {}{} outside simulation cell",
                                                 getNumberOfSteps(), particle.traits().name, atom_index,
                                                 group.traits().name, group_index));
        }
    }
}
void SanityCheck::checkMassCenter(const Space::Tgroup& group) {
    if (group.isMolecular() && !group.empty()) {
        const auto mass_center = Geometry::massCenter(group.begin(), group.end(), spc.geo.getBoundaryFunc(), -group.cm);
        const auto distance = spc.geo.vdist(group.cm, mass_center).norm();
        if (distance > mass_center_tolerance) {
            throw std::runtime_error(fmt::format("step {}: {}{} mass center out of sync by {:.3f} Å",
                                                 getNumberOfSteps(), group.traits().name, spc.getGroupIndex(group),
                                                 distance));
        }
    }
}
SanityCheck::SanityCheck(const json& j, Space& spc) : Analysisbase(spc, "sanity") {
    from_json(j);
    sample_interval = j.value("nstep", -1);
}

void AtomRDF::sampleDistance(const Point& position1, const Point& position2) {
    const auto distance = spc.geo.vdist(position1, position2);
    if (slicedir.sum() > 0) {
        if (distance.cwiseProduct((Point::Ones() - slicedir.cast<double>()).cwiseAbs()).norm() < thickness) {
            histogram(distance.cwiseProduct(slicedir.cast<double>()).norm())++;
        }
    } else {
        histogram(distance.norm())++;
    }
}

void AtomRDF::_sample() {
    mean_volume += spc.geo.getVolume(dimensions);

    auto active_particles = spc.activeParticles();
    for (auto i = active_particles.begin(); i != active_particles.end(); ++i) {
        for (auto j = i; ++j != active_particles.end();) {
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                sampleDistance(i->pos, j->pos);
            }
        }
    }
}
AtomRDF::AtomRDF(const json& j, Space& spc) : PairFunctionBase(spc, j, "atomrdf") {
    id1 = Faunus::findAtomByName(name1).id();
    id2 = Faunus::findAtomByName(name2).id();
}

void MoleculeRDF::_sample() {
    mean_volume += spc.geo.getVolume(dimensions);
    auto active_molecules = spc.groups | ranges::cpp20::views::filter([](const auto& group) { return !group.empty(); });

    for (auto i = active_molecules.begin(); i != active_molecules.end(); ++i) {
        for (auto j = i; ++j != active_molecules.end();) {
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                const auto distance = std::sqrt(spc.geo.sqdist(i->cm, j->cm));
                histogram(distance)++;
            }
        }
    }
}

MoleculeRDF::MoleculeRDF(const json& j, Space& spc) : PairFunctionBase(spc, j, "molrdf") {
    id1 = findMoleculeByName(name1).id();
    id2 = findMoleculeByName(name2).id();
}

void AtomDipDipCorr::_sample() {
    mean_volume += spc.geo.getVolume(dimensions);
    const auto particles = spc.activeParticles();

    auto sample_distance = [&](const auto& particle1, const auto& particle2, const Point& distance) {
        if (particle1.hasExtension() && particle2.hasExtension()) {
            const auto cosine_angle = particle1.getExt().mu.dot(particle2.getExt().mu);
            const auto r = distance.norm();
            hist2(r) += cosine_angle;
            histogram(r)++; // get g(r) for free
        }
    };

    for (auto i = particles.begin(); i != particles.end(); ++i) {
        for (auto j = i; ++j != particles.end();) {
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                const Point distance = spc.geo.vdist(i->pos, j->pos);
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
AtomDipDipCorr::AtomDipDipCorr(const json& j, Space& spc) : PairAngleFunctionBase(spc, j, "atomdipdipcorr") {
    id1 = findAtomByName(name1).id();
    id2 = findAtomByName(name2).id();
}

// =============== XTCtraj ===============

XTCtraj::XTCtraj(const json& j, Space& spc)
    : Analysisbase(spc, "xtcfile"), atom_filter([](const Particle&) { return true; }) {
    from_json(j);
    assert(atom_filter);
}

void XTCtraj::_to_json(json& j) const {
    j["file"] = writer->filename;
    if (not names.empty()) {
        j["molecules"] = names;
    }
}

void XTCtraj::_from_json(const json& j) {
    const auto file = MPI::prefix + j.at("file").get<std::string>();
    writer = std::make_shared<XTCWriter>(file);

    atom_filter = [](const Particle&) { return true; }; // default: save all particles
    names = j.value("molecules", std::vector<std::string>());
    if (!names.empty()) {
        molecule_ids = Faunus::names2ids(Faunus::molecules, names); // molecule types to save
        if (!molecule_ids.empty()) {
            std::sort(molecule_ids.begin(), molecule_ids.end()); // needed for binary_search
            atom_filter = [&](const Particle& particle) {
                for (const auto& group : spc.groups) {
                    if (group.contains(particle, true)) {
                        return std::binary_search(molecule_ids.begin(), molecule_ids.end(), group.id);
                    }
                }
                return false;
            };
        }
    }
}

void XTCtraj::_sample() {
    assert(atom_filter); // some gcc/clang/ubuntu/macos combinations wronly clear the `filter` function
    auto positions =
        spc.p | ranges::cpp20::views::filter(atom_filter) | ranges::cpp20::views::transform(&Particle::pos);
    writer->writeNext(spc.geo.getLength(), positions.begin(), positions.end());
}

// =============== MultipoleDistribution ===============

double MultipoleDistribution::g2g(const Space::Tgroup& group1, const Space::Tgroup& group2) {
    double energy = 0.0;
    for (const auto& particle_i : group1) {
        for (const auto& particle_j : group2) {
            energy += particle_i.charge * particle_j.charge / spc.geo.vdist(particle_i.pos, particle_j.pos).norm();
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
                                      energy.exact, u_tot, energy.ion_ion, energy.ion_dipole, energy.dipole_dipole,
                                      energy.ion_quadrupole, energy.dipole_dipole_correlation);
            }
        }
    }
}

void MultipoleDistribution::_sample() {
    for (auto& group1 : spc.findMolecules(ids[0])) {     // find active molecules
        for (auto& group2 : spc.findMolecules(ids[1])) { // find active molecules
            if (&group1 != &group2) {
                const auto a = Faunus::toMultipole(group1, spc.geo.getBoundaryFunc());
                const auto b = Faunus::toMultipole(group2, spc.geo.getBoundaryFunc());
                const auto distance = spc.geo.vdist(group1.cm, group2.cm);
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

MultipoleDistribution::MultipoleDistribution(const json& j, Space& spc) : Analysisbase(spc, "Multipole Distribution") {
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
    const auto cm = Geometry::massCenter(slice.begin(), slice.end(), spc.geo.getBoundaryFunc());
    const auto I = Geometry::inertia(slice.begin(), slice.end(), cm, spc.geo.getBoundaryFunc());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(I);
    return esf.eigenvalues();
}
void AtomInertia::_sample() {
    if (output_stream) {
        output_stream << getNumberOfSteps() << " " << compute().transpose() << "\n";
    }
}
AtomInertia::AtomInertia(const json& j, Space& spc) : Analysisbase(spc, "Atomic Inertia Eigenvalues") {
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
}
InertiaTensor::Data InertiaTensor::compute() {
    const auto& group = spc.groups.at(group_index);
    const Space::Tgroup subgroup(group.begin() + particle_range[0], group.begin() + particle_range[1] + 1);
    InertiaTensor::Data d;
    auto I = Geometry::inertia(subgroup.begin(), subgroup.end(), group.cm, spc.geo.getBoundaryFunc());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(I);
    d.eigen_values = esf.eigenvalues();
    std::ptrdiff_t i_eival;
    d.eigen_values.minCoeff(&i_eival);
    d.principle_axis = esf.eigenvectors().col(i_eival).real(); // eigenvector corresponding to the smallest eigenvalue
    return d;
}
void InertiaTensor::_sample() {
    InertiaTensor::Data d = compute();
    if (output_stream) {
        output_stream << getNumberOfSteps() << " " << d.eigen_values.transpose() << " " << d.principle_axis.transpose()
                      << "\n";
    }
}
InertiaTensor::InertiaTensor(const json& j, Space& spc) : Analysisbase(spc, "Inertia Tensor") {
    from_json(j);
    filename = MPI::prefix + j.at("file").get<std::string>();
    output_stream.open(filename);
    group_index = j.at("index").get<size_t>();
    particle_range =
        j.value("indexes", std::vector<size_t>({0, spc.groups[group_index].size()})); // whole molecule by default
}
void InertiaTensor::_to_disk() {
    if (output_stream) {
        output_stream.flush(); // empty buffer
    }
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
    Space::Tgroup subgroup(group.begin() + particle_range[0], group.begin() + particle_range[1] + 1);
    const auto mass_center = use_molecular_mass_center
                                 ? group.cm
                                 : Geometry::massCenter(subgroup.begin(), subgroup.end(), spc.geo.getBoundaryFunc());

    MultipoleMoments::Data multipole;
    Tensor quadrupole; // quadrupole tensor
    quadrupole.setZero();
    for (const auto& particle : subgroup) {
        Point position = particle.pos - mass_center;
        spc.geo.boundary(position);
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
MultipoleMoments::MultipoleMoments(const json& j, Space& spc) : Analysisbase(spc, "Multipole Moments") {
    from_json(j);
    try {
        use_molecular_mass_center = j.value("mol_cm", true); // use the mass center of the whole molecule
        filename = MPI::prefix + j.at("file").get<std::string>();
        output_stream.open(filename); // output file

        group_index = j.at("index").get<int>();
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
    if (data.gyration_radius.cnt > 0) {
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
    auto molecules = spc.findMolecules(molid, Space::ACTIVE);
    const auto num_molecules = std::distance(molecules.begin(), molecules.end());

    if (num_molecules > 1 && tensor_output_stream) {
        throw std::runtime_error("tensor output `file` cannot be used with multiple molecules");
    }

    for (const auto& group : molecules) {
        if (group.size() >= 2) { // two or more particles required to form a polymer
            const auto gyration_tensor =
                Geometry::gyration(group.begin(), group.end(), group.cm, spc.geo.getBoundaryFunc());
            const auto principal_moment = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(gyration_tensor).eigenvalues();
            const auto gyration_radius_squared = gyration_tensor.trace();
            const auto end_to_end_squared = spc.geo.sqdist(group.begin()->pos, std::prev(group.end())->pos);

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

PolymerShape::PolymerShape(const json& j, Space& spc) : Analysisbase(spc, "Polymer Shape") {
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
        tensor_output_stream = IO::openCompressedOutputStream(MPI::prefix + filename);
        *tensor_output_stream << "# step Rg xx xy xz xy yy yz xz yz zz\n";
    }
}

void AtomProfile::_from_json(const json& j) {
    origin = j.value("origo", Point(0, 0, 0));
    dir = j.value("dir", dir);
    file = j.at("file").get<std::string>();
    names = j.at("atoms").get<decltype(names)>();                            // atom names
    const auto vec_of_ids = Faunus::names2ids(Faunus::atoms, names);         // names --> molids
    atom_id_selection = std::set<int>(vec_of_ids.begin(), vec_of_ids.end()); // copy vector to set
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
        origin =
            Geometry::massCenter(mass_center_particles.begin(), mass_center_particles.end(), spc.geo.getBoundaryFunc());
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
    const Point distance = spc.geo.vdist(position, origin);
    return distance.cwiseProduct(dir.cast<double>()).norm();
}

AtomProfile::AtomProfile(const json& j, Space& spc) : Analysisbase(spc, "atomprofile") { from_json(j); }

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
        z_offset =
            Geometry::massCenter(mass_center_particles.begin(), mass_center_particles.end(), spc.geo.getBoundaryFunc())
                .z();
    }

    auto filtered_particles = spc.activeParticles() | ranges::cpp20::views::filter([&](const auto& particle) {
                                  return std::find(atom_ids.begin(), atom_ids.end(), particle.id) != atom_ids.end();
                              });

    for (const auto& particle : filtered_particles) {
        histogram(particle.pos.z() - z_offset)++;
    }
}
SlicedDensity::SlicedDensity(const json& j, Space& spc) : Analysisbase(spc, "sliceddensity") { from_json(j); }

void SlicedDensity::_to_disk() {
    if (std::ofstream f(MPI::prefix + file); f and number_of_samples > 0) {
        f << "# z rho/M\n";
        const Point box_length = spc.geo.getLength();
        const auto half_z_length = 0.5 * box_length.z();
        const auto volume = box_length.x() * box_length.y() * dz;
        for (double z = -half_z_length; z <= half_z_length; z += dz) {
            f << fmt::format("{:.6E} {:.6E}\n", z, histogram(z) / volume / number_of_samples * 1e27 / pc::Nav);
        }
    }
}
void ChargeFluctuations::_sample() {
    auto filtered_molecules = spc.findMolecules(mol_iter->id(), Space::ACTIVE);
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
        auto molecules = spc.findMolecules(mol_iter->id(), Space::ALL);
        if (not ranges::cpp20::empty(molecules)) {
            const auto particles_with_avg_charges = averageChargeParticles(*molecules.begin());
            FormatPQR::save(MPI::prefix + filename, particles_with_avg_charges, spc.geo.getLength());
        }
    }
}

ParticleVector ChargeFluctuations::averageChargeParticles(const Space::Tgroup& group) {
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
        added_particle.pos = it->pos - group.cm;
        spc.geo.boundary(added_particle.pos);
        particle_index++;
    }
    return particles;
}

/**
 * @todo replace `mol_iter` with simple molid integer
 */
ChargeFluctuations::ChargeFluctuations(const json& j, Space& spc) : Analysisbase(spc, "chargefluctuations") {
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
            const auto particle = Faunus::toMultipole(group, spc.geo.getBoundaryFunc());
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
Multipole::Multipole(const json& j, Space& spc) : Analysisbase(spc, "multipole") { from_json(j); }
void ScatteringFunction::_sample() {
    scatter_positions.clear();
    for (int id : molecule_ids) {                         // loop over molecule names
        for (const auto& group : spc.findMolecules(id)) { // loop over groups
            if (mass_center_scattering && group.isMolecular()) {
                scatter_positions.push_back(group.cm);
            } else {
                for (const auto& particle : group) { // loop over particle index in group
                    scatter_positions.push_back(particle.pos);
                }
            }
        }
    }

    // zero-padded suffix to use with `save_after_sample`
    const auto suffix = fmt::format("{:07d}", number_of_samples);

    switch (scheme) {
    case DEBYE:
        debye->sample(scatter_positions, spc.geo.getVolume());
        if (save_after_sample) {
            IO::write(filename + "." + suffix, debye->getIntensity());
        }
        break;
    case EXPLICIT_PBC:
        explicit_average_pbc->sample(scatter_positions, spc.geo.getLength());
        if (save_after_sample) {
            IO::write(filename + "." + suffix, explicit_average_pbc->getSampling());
        }
        break;
    case EXPLICIT_IPBC:
        explicit_average_ipbc->sample(scatter_positions, spc.geo.getLength());
        if (save_after_sample) {
            IO::write(filename + "." + suffix, explicit_average_ipbc->getSampling());
        }
        break;
    }
}
void ScatteringFunction::_to_json(json& j) const {
    j = {{"molecules", molecule_names}, {"com", mass_center_scattering}};
    switch (scheme) {
    case DEBYE:
        j["scheme"] = "debye";
        std::tie(j["qmin"], j["qmax"], std::ignore) = debye->getQMeshParameters();
        break;
    case EXPLICIT_PBC:
        j["scheme"] = "explicit";
        j["pmax"] = explicit_average_pbc->getQMultiplier();
        j["ipbc"] = false;
        break;
    case EXPLICIT_IPBC:
        j["scheme"] = "explicit";
        j["pmax"] = explicit_average_ipbc->getQMultiplier();
        j["ipbc"] = true;
        break;
    }
}

ScatteringFunction::ScatteringFunction(const json& j, Space& spc) try : Analysisbase(spc, "scatter") {
    from_json(j);

    mass_center_scattering = j.value("com", true);
    save_after_sample = j.value("stepsave", false);
    filename = j.at("file").get<std::string>();
    molecule_names = j.at("molecules").get<decltype(molecule_names)>();
    molecule_ids = Faunus::names2ids(molecules, molecule_names);

    const auto cuboid = std::dynamic_pointer_cast<Geometry::Cuboid>(spc.geo.asSimpleGeometry());

    if (const auto scheme_str = j.value("scheme", "explicit"s); scheme_str == "debye") {
        scheme = DEBYE;
        debye = std::make_shared<Scatter::DebyeFormula<Tformfactor>>(j);
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
            scheme = EXPLICIT_IPBC;
            explicit_average_ipbc = std::make_shared<Scatter::StructureFactorIPBC<>>(pmax);
        } else {
            scheme = EXPLICIT_PBC;
            explicit_average_pbc = std::make_shared<Scatter::StructureFactorPBC<>>(pmax);
        }
    } else {
        throw ConfigurationError("unknown scheme");
    }
} catch (std::exception& e) {
    throw ConfigurationError("scatter: "s + e.what());
}

void ScatteringFunction::_to_disk() {
    switch (scheme) {
    case DEBYE:
        IO::write(filename, debye->getIntensity());
        break;
    case EXPLICIT_PBC:
        IO::write(filename, explicit_average_pbc->getSampling());
        break;
    case EXPLICIT_IPBC:
        IO::write(filename, explicit_average_ipbc->getSampling());
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
    if (std::fabs(perturbation_distance) > 0.0) {
        if (auto mollist = spc.findMolecules(molid, Space::ACTIVE); !ranges::cpp20::empty(mollist)) {
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
}
void VirtualTranslate::writeToFileStream(const double energy_change) const {
    if (stream) { // file to disk?
        const auto mean_force = -meanFreeEnergy() / perturbation_distance;
        *stream << fmt::format("{:d} {:.3f} {:.10f} {:.6f}\n", getNumberOfSteps(), perturbation_distance, energy_change,
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
double VirtualTranslate::momentarilyPerturb(Space::Tgroup& group) {
    change.groups.at(0).index = spc.getGroupIndex(group);
    const auto old_energy = pot.energy(change);
    const Point displacement_vector = perturbation_distance * perturbation_direction;
    group.translate(displacement_vector, spc.geo.getBoundaryFunc()); // temporarily translate group
    const auto new_energy = pot.energy(change);
    group.translate(-displacement_vector, spc.geo.getBoundaryFunc()); // restore original position
    return new_energy - old_energy;
}

void VirtualTranslate::_to_json(json& j) const {
    if (number_of_samples > 0 && std::fabs(perturbation_distance) > 0.0) {
        j = {{"dL", perturbation_distance},
             {"force", std::log(mean_exponentiated_energy_change) / perturbation_distance},
             {"dir", perturbation_direction}};
    }
}
VirtualTranslate::VirtualTranslate(const json& j, Space& spc, Energy::Energybase& pot)
    : PerturbationAnalysisBase("virtualtranslate", pot, spc, j.value("file", ""s)) {
    change.groups.resize(1);
    change.groups.front().internal = false;
    from_json(j);
}

SpaceTrajectory::SpaceTrajectory(const json& j, Space& spc)
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
} // namespace Faunus::Analysis
