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

void to_json(json &j, const Analysisbase &base) { base.to_json(j); }

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
    number_of_steps++;
    if (sample_interval > 0 && number_of_steps > number_of_skipped_steps) {
        if ((number_of_steps % sample_interval) == 0) {
            number_of_samples++;
            timer.start();
            _sample();
            timer.stop();
        }
    }
}

void Analysisbase::from_json(const json &j) {
    number_of_skipped_steps = j.value("nskip", 0);
    sample_interval = j.value("nstep", 0);
    _from_json(j);
}

void Analysisbase::to_json(json &json_output) const {
    assert(not name.empty());
    try {
        auto &j = json_output[name];
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
        if (not cite.empty()) {
            j["reference"] = cite;
        }
    } catch (std::exception &e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}

void Analysisbase::_to_json(json &) const {}

void Analysisbase::_from_json(const json &) {}

int Analysisbase::getNumberOfSteps() const { return number_of_steps; }

void SystemEnergy::normalize() {
    double sum = energy_histogram.sumy();
    for (auto &i : energy_histogram.getMap()) {
        i.second = i.second / sum;
    }
}

void SystemEnergy::_sample() {
    auto energies = energyFunc(); // current energy from all terms in Hamiltonian
    double total_energy = std::accumulate(energies.begin(), energies.end(), 0.0);
    if (std::isfinite(total_energy)) {
        mean_energy += total_energy;
        mean_squared_energy += total_energy * total_energy;
    }
    *output_stream << getNumberOfSteps() << separator << total_energy;
    for (auto energy : energies) {
        *output_stream << separator << energy;
    }
    *output_stream << "\n";
    // ehist(tot)++;
}

void SystemEnergy::_to_json(json &j) const {
    j = {{"file", file_name}, {"init", initial_energy}, {"final", energyFunc()}};
    if (number_of_samples > 0) {
        j["mean"] = mean_energy.avg();
        j["Cv/kB"] = mean_squared_energy.avg() - std::pow(mean_energy.avg(), 2);
    }
    _roundjson(j, 5);
    // normalize();
    // ehist.save( "distofstates.dat" );
}

void SystemEnergy::_from_json(const json &j) {
    file_name = MPI::prefix + j.at("file").get<std::string>();
    if (output_stream = IO::openCompressedOutputStream(file_name); output_stream == nullptr) {
        throw std::runtime_error(name + ": cannot open output file " + file_name);
    } else {
        if (auto suffix = file_name.substr(file_name.find_last_of(".") + 1); suffix == "csv") {
            separator = ",";
        } else {
            separator = " ";
            *output_stream << "#";
        }
        *output_stream << "total";
        for (auto &name : names_of_energy_terms) {
            *output_stream << separator << name;
        }
        *output_stream << "\n";
        output_stream->precision(16);
    }
}

SystemEnergy::SystemEnergy(const json &j, Energy::Hamiltonian &pot) {
    assert(!pot.vec.empty());
    name = "systemenergy";
    for (auto i : pot.vec) {
        names_of_energy_terms.push_back(i->name);
    }
    from_json(j);
    energyFunc = [&pot]() {
        Change change;
        change.all = true;
        std::vector<double> energies;
        std::transform(pot.vec.begin(), pot.vec.end(), std::back_inserter(energies),
                       [&](auto i) { return i->energy(change); });
        return energies;
    };
    energy_histogram.setResolution(0.25);
    auto energies = energyFunc();
    initial_energy = std::accumulate(energies.begin(), energies.end(), 0.0); // initial energy
}

void SystemEnergy::_to_disk() {
    if (*output_stream) {
        output_stream->flush(); // empty buffer
    }
}

void SaveState::_to_json(json &j) const { j["file"] = filename; }

void SaveState::_sample() {
    assert(sample_interval >= 0);
    // tag filename with step number:
    if (use_numbered_files) {
        auto numbered_filename = filename;
        numbered_filename.insert(filename.find_last_of("."), "_"s + std::to_string(getNumberOfSteps()));
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

SaveState::SaveState(json j, Space &spc) {
    name = "savestate";

    if (j.count("nstep") == 0) { // by default, disable _sample() and
        j["nstep"] = -1;         // store only when _to_disk() is called
    }
    from_json(j);

    save_random_number_generator_state = j.value("saverandom", false);
    filename = MPI::prefix + j.at("file").get<std::string>();
    use_numbered_files = !j.value("overwrite", false);
    convert_hexagonal_prism_to_cuboid = j.value("convert_hexagon", false);

    if (auto suffix = filename.substr(filename.find_last_of(".") + 1); suffix == "aam") {
        writeFunc = [&](auto &file) { FormatAAM::save(file, spc.p); };
    } else if (suffix == "gro") {
        writeFunc = [&](auto &file) { FormatGRO::save(file, spc); };
    } else if (suffix == "pqr") {
        writeFunc = [&](auto &file) {
            if (convert_hexagonal_prism_to_cuboid) {
                auto hexagonal_prism =
                    std::dynamic_pointer_cast<Geometry::HexagonalPrism>(spc.geo.asSimpleGeometry());
                if (hexagonal_prism) {
                    faunus_logger->debug("creating cuboidal PQR from hexagonal prism");
                    const auto &[cuboid, particles] =
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
        writeFunc = [&](auto &file) { FormatXYZ::save(file, spc.p, spc.geo.getLength()); };
    } else if (suffix == "json") { // JSON state file
        writeFunc = [&](auto &file) {
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
        writeFunc = [&](auto &file) {
            if (std::ofstream f(file, std::ios::binary); f) {
                json j;
                Faunus::to_json(j, spc);
                if (save_random_number_generator_state) {
                    j["random-move"] = Move::MoveBase::slump;
                    j["random-global"] = Faunus::random;
                }
                auto v = json::to_ubjson(j); // json --> binary
                f.write((const char *)v.data(), v.size() * sizeof(decltype(v)::value_type));
            }
        };
    } else {
        throw ConfigurationError("unknown file extension for '{}'", filename);
    }
}

PairFunctionBase::PairFunctionBase(const json &j) { from_json(j); }

void PairFunctionBase::_to_json(json &j) const {
    j = {{"dr", dr / 1.0_angstrom}, {"name1", name1},        {"name2", name2}, {"file", file}, {"dim", dim},
         {"slicedir", slicedir},    {"thickness", thickness}};
    if (Rhypersphere > 0)
        j["Rhyper"] = Rhypersphere;
}

void PairFunctionBase::_from_json(const json &j) {
    file = j.at("file");
    name1 = j.at("name1");
    name2 = j.at("name2");
    dim = j.value("dim", 3);
    dr = j.value("dr", 0.1) * 1.0_angstrom;
    slicedir = j.value("slicedir", slicedir);
    thickness = j.value("thickness", 0);
    hist.setResolution(dr, 0);
    Rhypersphere = j.value("Rhyper", -1.0);
}
void PairFunctionBase::_to_disk() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        double Vr = 1, sum = hist.sumy();
        hist.stream_decorator = [&](std::ostream &o, double r, double N) {
            if (dim == 3)
                Vr = 4 * pc::pi * std::pow(r, 2) * dr;
            else if (dim == 2) {
                Vr = 2 * pc::pi * r * dr;
                if (Rhypersphere > 0)
                    Vr = 2.0 * pc::pi * Rhypersphere * std::sin(r / Rhypersphere) * dr;
            } else if (dim == 1)
                Vr = dr;
            if (Vr > 0)
                o << r << " " << N * V / (Vr * sum) << "\n";
        };
        f << hist;
    }
}

PairAngleFunctionBase::PairAngleFunctionBase(const json &j) : PairFunctionBase(j) { from_json(j); }

void PairAngleFunctionBase::_to_disk() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        hist2.stream_decorator = [&](std::ostream &o, double r, double N) {
            o << r << " " << N << "\n";
        };
        f << hist2;
    }
    file = file + "gofr.dat"; // name file where free g(r) is saved, and make sure that the file is not overwritten by
                              // base-destructor
}

void PairAngleFunctionBase::_from_json(const json &) { hist2.setResolution(dr, 0); }

void PerturbationAnalysisBase::_to_disk() {
    if (output_stream) {
        output_stream->flush(); // empty buffer
    }
}
PerturbationAnalysisBase::PerturbationAnalysisBase(Energy::Energybase& pot, Space& spc, const std::string& filename)
    : spc(spc), pot(pot), filename(filename) {
    if (!filename.empty()) {
        this->filename = MPI::prefix + filename;
        output_stream = IO::openCompressedOutputStream(this->filename, true); // throws if error
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
    if (output_stream) { // write to output file if open
        *output_stream << fmt::format("{:d} {:.3f} {:.10f} {:.6f} {:.6f}", getNumberOfSteps(), volume_displacement,
                                      energy_change, exp(-energy_change),
                                      -meanFreeEnergy() / volume_displacement);

        // if anisotropic scaling, add an extra column with area or length perturbation
        if (volume_scaling_method == Geometry::XY) {
            const auto box_length = spc.geo.getLength();
            const auto area_change = box_length.x() * box_length.y() * (scale.x() * scale.y() - 1.0);
            *output_stream << fmt::format(" {:.6f}", area_change);
        } else if (volume_scaling_method == Geometry::Z) {
            const auto box_length = spc.geo.getLength();
            const auto length_change = box_length.z() * (scale.z() - 1.0);
            *output_stream << fmt::format(" {:.6f}", length_change);
        }
        *output_stream << "\n"; // trailing newline
    }
}

void VirtualVolumeMove::_from_json(const json& j) {
    volume_displacement = j.at("dV").get<double>();
    volume_scaling_method = j.value("scaling", Geometry::ISOTROPIC);
    if (volume_scaling_method == Geometry::ISOCHORIC) {
        throw std::runtime_error(name + ": isochoric volume scaling not allowed");
    }
}

void VirtualVolumeMove::_to_json(json &j) const {
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
    : PerturbationAnalysisBase(pot, spc, j.value("file", ""s)) {
    name = "virtualvolume";
    cite = "doi:10.1063/1.472721";
    from_json(j);
    change.dV = true;
    change.all = true;
    if (output_stream) {
        *output_stream << "# steps dV/" + u8::angstrom + u8::cubed + " du/kT exp(-du/kT) <Pex>/kT/" + u8::angstrom +
                          u8::cubed;
        // if non-isotropic scaling, add another column with dA or dL
        if (volume_scaling_method == Geometry::XY) {
            *output_stream << " dA/" + u8::angstrom + u8::squared;
        } else if (volume_scaling_method == Geometry::Z) {
            *output_stream << " dL/" + u8::angstrom;
        }
        *output_stream << "\n"; // trailing newline
    }
}

void MolecularConformationID::_sample() {
    auto molecules = spc.findMolecules(molid, Space::ACTIVE);
    for (auto &group : molecules) {
        histogram[group.confid]++;
    }
}
void MolecularConformationID::_to_json(json &j) const { j["histogram"] = histogram; }

MolecularConformationID::MolecularConformationID(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "moleculeconformation";
    const std::string molname = j.at("molecule");
    molid = findMoleculeByName(molname).id();
}

void QRtraj::_sample() { write_to_file(); }

void QRtraj::_to_json(json &j) const { j = {{"file", filename}}; }

QRtraj::QRtraj(const json &j, Space &spc) {
    from_json(j);
    name = "qrfile";
    filename = MPI::prefix + j.value("file", "qrtraj.dat"s);
    if (stream = IO::openCompressedOutputStream(filename); !*stream) { // may be gzip compressed
        throw std::runtime_error("could not open create "s + filename);
    }
    write_to_file = [&groups = spc.groups, &stream = stream]() {
        for (auto &g : groups) {
            for (auto it = g.begin(); it != g.trueend(); ++it) { // loop over *all* particles
                if (it < g.end()) {                              // active particles...
                    *stream << fmt::format("{} {} ", it->charge, atoms[it->id].sigma * 0.5);
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

void CombinedAnalysis::sample() {
    for (auto &ptr : this->vec)
        ptr->sample();
}

void CombinedAnalysis::to_disk() {
    for (auto &ptr : this->vec)
        ptr->to_disk();
}

CombinedAnalysis::CombinedAnalysis(const json& j, Space& spc, Energy::Hamiltonian& pot) {
    if (!j.is_array()) {
        throw ConfigurationError("analysis: json array expected");
    }
    for (auto& j_analysis : j) {
        try {
            const auto& [key, j_params] = jsonSingleItem(j_analysis);
            try {
                if (key == "atomprofile") {
                    emplace_back<AtomProfile>(j_params, spc);
                } else if (key == "atomrdf") {
                    emplace_back<AtomRDF>(j_params, spc);
                } else if (key == "atomdipdipcorr") {
                    emplace_back<AtomDipDipCorr>(j_params, spc);
                } else if (key == "density") {
                    emplace_back<Density>(j_params, spc);
                } else if (key == "chargefluctuations") {
                    emplace_back<ChargeFluctuations>(j_params, spc);
                } else if (key == "molrdf") {
                    emplace_back<MoleculeRDF>(j_params, spc);
                } else if (key == "multipole") {
                    emplace_back<Multipole>(j_params, spc);
                } else if (key == "atominertia") {
                    emplace_back<AtomInertia>(j_params, spc);
                } else if (key == "inertia") {
                    emplace_back<InertiaTensor>(j_params, spc);
                } else if (key == "moleculeconformation") {
                    emplace_back<MolecularConformationID>(j_params, spc);
                } else if (key == "multipolemoments") {
                    emplace_back<MultipoleMoments>(j_params, spc);
                } else if (key == "multipoledist") {
                    emplace_back<MultipoleDistribution>(j_params, spc);
                } else if (key == "polymershape") {
                    emplace_back<PolymerShape>(j_params, spc);
                } else if (key == "qrfile") {
                    emplace_back<QRtraj>(j_params, spc);
                } else if (key == "reactioncoordinate") {
                    emplace_back<FileReactionCoordinate>(j_params, spc);
                } else if (key == "sanity") {
                    emplace_back<SanityCheck>(j_params, spc);
                } else if (key == "savestate") {
                    emplace_back<SaveState>(j_params, spc);
                } else if (key == "scatter") {
                    emplace_back<ScatteringFunction>(j_params, spc);
                } else if (key == "sliceddensity") {
                    emplace_back<SlicedDensity>(j_params, spc);
                } else if (key == "systemenergy") {
                    emplace_back<SystemEnergy>(j_params, pot);
                } else if (key == "virtualvolume") {
                    emplace_back<VirtualVolumeMove>(j_params, spc, pot);
                } else if (key == "virtualtranslate") {
                    emplace_back<VirtualTranslate>(j_params, spc, pot);
                } else if (key == "widom") {
                    emplace_back<WidomInsertion>(j_params, spc, pot);
                } else if (key == "xtcfile") {
                    emplace_back<XTCtraj>(j_params, spc);
                } else if (key == "spacetraj") {
                    emplace_back<SpaceTrajectory>(j_params, spc.groups);
                } else {
                    throw ConfigurationError("unknown analysis");
                }
                // append additional analysis to the if-chain
            } catch (std::exception& e) {
                usageTip.pick(key);
                throw ConfigurationError("'{}': {}", key, e.what());
            }
        } catch (std::exception& e) {
            throw ConfigurationError("analysis: {}", e.what()).attachJson(j_analysis);
        }
    }
}

void FileReactionCoordinate::_to_json(json &j) const {
    json rcjson = *rc; // invoke to_json(...)
    if (rcjson.count(type) == 0)
        throw std::runtime_error("error writing json for reaction coordinate");
    j = rcjson[type];
    j["type"] = type;
    j["file"] = filename;
    j.erase("range");      // these are for penalty function
    j.erase("resolution"); // use only, so no need to show
    if (number_of_samples > 0)
        j["average"] = avg.avg();
}

void FileReactionCoordinate::_sample() {
    if (*stream) {
        double val = (*rc)();
        avg += val;
        (*stream) << fmt::format("{} {:.6f} {:.6f}\n", getNumberOfSteps(), val, avg.avg());
    }
}

FileReactionCoordinate::FileReactionCoordinate(const json &j, Space &spc) {
    from_json(j);
    name = "reactioncoordinate";
    filename = MPI::prefix + j.at("file").get<std::string>();
    if (stream = IO::openCompressedOutputStream(filename); not*stream) {
        throw std::runtime_error("could not open create "s + filename);
    }
    type = j.at("type").get<std::string>();
    rc = ReactionCoordinate::createReactionCoordinate({{type, j}}, spc);
}

void FileReactionCoordinate::_to_disk() {
    if (*stream) {
        stream->flush(); // empty buffer
    }
}

/**
 * This searches for an inactive group of type `molid`
 * and prepares the `change` object for energy evaluation
 */
void WidomInsertion::selectGhostGroup() {
    change.clear();
    auto inactive_groups = spc.findMolecules(molid, Space::INACTIVE); // list of inactive molecules
    if (!ranges::cpp20::empty(inactive_groups)) {                     // did we find any?
        const auto& group = *inactive_groups.begin();                 // select first group
        if (group.empty() && group.capacity() > 0) {                  // must be inactive and have a non-zero capacity
            auto& group_changes = change.groups.emplace_back();
            group_changes.index = spc.getGroupIndex(group); // group index in space
            group_changes.all = true;
            group_changes.internal = group.isAtomic(); // internal energy only if non-molecular
            return;
        }
    }
    std::runtime_error(fmt::format("{}: no inactive {} groups available", name, Faunus::molecules.at(molid).name));
}

void WidomInsertion::_sample() {
    selectGhostGroup();
    auto& group = spc.groups.at(change.groups.at(0).index); // inactive "ghost" group
    group.resize(group.capacity());                         // activate ghost
    ParticleVector particles;                               // particles to insert
    for (int cnt = 0; cnt < number_of_insertions; ++cnt) {
        particles = inserter->operator()(spc.geo, Faunus::molecules[molid], spc.p); // random pos&orientation
        updateGroup(group, particles);
        const auto energy_change = pot.energy(change); // in kT
        collectWidomAverage(energy_change);
    }
    group.resize(0); // de-activate group
}

void WidomInsertion::updateGroup(Space::Tgroup& group, const ParticleVector& particles) {
    assert(particles.size() == group.size());
    std::copy(particles.begin(), particles.end(), group.begin()); // copy to ghost group
    if (absolute_z_coords) {
        std::for_each(group.begin(), group.end(), [](Particle& i) { i.pos.z() = std::fabs(i.pos.z()); });
    }
    if (group.isMolecular()) { // update molecular mass-center for molecular groups
        group.cm = Geometry::massCenter(group.begin(), group.end(), this->spc.geo.getBoundaryFunc(), -group.begin()->pos);
    }
}

void WidomInsertion::_to_json(json &j) const {
    const auto number_of_insertions = mean_exponentiated_energy_change.cnt;
    if (number_of_insertions > 0) {
        const double excess_chemical_potential = meanFreeEnergy();
        j = {{"molecule", Faunus::molecules[molid].name},
             {"insertions", number_of_insertions},
             {"absz", absolute_z_coords},
             {"insertscheme", *inserter},
             {u8::mu + "/kT", {{"excess", excess_chemical_potential}}}};
    }
}

void WidomInsertion::_from_json(const json &j) {
    number_of_insertions = j.at("ninsert").get<int>();
    absolute_z_coords = j.value("absz", false);
    if (auto ptr = std::dynamic_pointer_cast<RandomInserter>(inserter); ptr) {
        ptr->dir = j.value("dir", Point({1, 1, 1}));
    } // set insert directions for RandomInserter

    const auto molecule_name = j.at("molecule").get<std::string>();
    molid = findMoleculeByName(molecule_name).id();
}

WidomInsertion::WidomInsertion(const json& j, Space& spc, Energy::Hamiltonian& pot)
    : PerturbationAnalysisBase(pot, spc) {
    name = "widom";
    cite = "doi:10/dkv4s6";
    inserter = std::make_shared<RandomInserter>();
    from_json(j);
}

void Density::_sample() {
    auto volume = updateVolumeStatistics();
    auto [atom_count, molecular_group_count] = countAtomsAndMolecules();

    for (auto [atomid, N] : atom_count) {
        mean_atom_density[atomid] += N / volume;
    }

    for (auto [molid, N] : molecular_group_count) {
        mean_molecule_density[molid] += N / volume;
        molecular_group_probability_density[molid](N)++;
    }

    for (const auto &reaction : reactions) { // in case of reactions involving atoms (swap moves)
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
    for (const auto &group : spc.groups) {
        if (group.atomic) {
            for (auto particle = group.begin(); particle < group.trueend(); ++particle) {
                atom_count[particle->id] = 0;
            }
        } else {
            molecular_group_count[group.id] = 0;
        }
    }

    for (const auto &group : spc.groups) {
        if (group.atomic) {
            for (const auto &particle : group) {
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
    mean_cubic_root_of_volume += cbrt(volume);
    mean_inverse_volume += 1.0 / volume;
    return volume;
}

void Density::_to_json(json &j) const {
    j["<V>"] = mean_volume.avg();
    j["<∛V>"] = mean_cubic_root_of_volume.avg();
    j["∛<V>"] = std::cbrt(mean_volume.avg());
    j["<1/V>"] = mean_inverse_volume.avg();

    auto &j_atomic = j["atomic"] = json::object();
    auto &j_molecular = j["molecular"] = json::object();

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

Density::Density(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "density";
    for (const auto &molecule : Faunus::molecules) {
        if (molecule.atomic) {
            atomic_group_probability_density[molecule.id()].setResolution(1, 0);
        } else {
            molecular_group_probability_density[molecule.id()].setResolution(1, 0);
        }
    }
    for (const auto &reaction : Faunus::reactions) { // in case of reactions involving atoms (swap moves)
        const auto reactive_atomic_species = reaction.getReactantsAndProducts().first;
        for (auto atomid : reactive_atomic_species) {
            atomswap_probability_density[atomid].setResolution(1, 0);
        }
    }
}

/**
 * Write histograms to disk
 */
void Density::writeTable(const std::string &atom_or_molecule_name, Ttable &table) {
    table.stream_decorator = [&table](std::ostream &stream, int N, double counts) {
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
    for (auto &reaction : Faunus::reactions) {
        const auto reactive_atomic_species = reaction.getReactantsAndProducts().first;
        for (auto atomid : reactive_atomic_species) {
            writeTable(Faunus::atoms[atomid].name, atomswap_probability_density[atomid]);
        }
    }
}

void SanityCheck::_sample() {
    // loop over all groups
    for (auto &g : spc.groups) {
        // check if particles are inside container
        for (auto &i : g) // loop over active particles
            if (spc.geo.collision(i.pos))
                throw std::runtime_error("step "s + std::to_string(number_of_samples) + ": index " +
                                         std::to_string(&i - &(*g.begin())) + " of group " +
                                         std::to_string(std::distance(spc.groups.begin(), spc.findGroupContaining(i))) +
                                         " outside container");

        // The groups must exactly contain all particles in `p`
        size_t i = 0;
        for (auto &g : spc.groups)
            for (auto it = g.begin(); it != g.trueend(); ++it)
                if (&*it != &spc.p.at(i++))
                    throw std::runtime_error("group vector out of sync");
        assert(i == spc.p.size());

        // check if molecular mass centers are correct
        if (not g.atomic)
            if (not g.empty()) {
                Point cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.cm);
                double sqd = spc.geo.sqdist(g.cm, cm);
                if (sqd > 1e-6) {
                    std::cerr << "step:      " << number_of_samples << std::endl
                              << "molecule:  " << &g - &*spc.groups.begin() << std::endl
                              << "dist:      " << sqrt(sqd) << std::endl
                              << "g.cm:      " << g.cm.transpose() << std::endl
                              << "actual cm: " << cm.transpose() << std::endl;
                    FormatPQR::save(MPI::prefix + "sanity-" + std::to_string(number_of_samples) + ".pqr", spc.p,
                                    spc.geo.getLength());
                    throw std::runtime_error("mass center-out-of-sync");
                }
            }
    }
}
SanityCheck::SanityCheck(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "sanity";
    sample_interval = j.value("nstep", -1);
}
void AtomRDF::_sample() {
    V += spc.geo.getVolume(dim);
    auto active = spc.activeParticles();
    for (auto i = active.begin(); i != active.end(); ++i)
        for (auto j = i; ++j != active.end();)
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                Point rvec = spc.geo.vdist(i->pos, j->pos);
                if (slicedir.sum() > 0) {
                    if (rvec.cwiseProduct(slicedir.cast<double>()).norm() < thickness) {
                        // rvec = rvec.cwiseProduct( Point(1.,1.,1.) - slice.cast<double>() );
                        hist(rvec.norm())++;
                    }
                } else {
                    hist(rvec.norm())++;
                }
            }
}
AtomRDF::AtomRDF(const json &j, Space &spc) : PairFunctionBase(j), spc(spc) {
    name = "atomrdf";
    id1 = findAtomByName(name1).id();
    id2 = findAtomByName(name2).id();
}
void MoleculeRDF::_sample() {
    V += spc.geo.getVolume(dim);
    auto mollist1 = spc.findMolecules(id1, Space::ACTIVE);
    auto mollist2 = spc.findMolecules(id2, Space::ACTIVE);
    // auto mollist = ranges::views::concat(mollist1, mollist2);

    std::vector<std::reference_wrapper<Tspace::Tgroup>> mollist;
    mollist.insert(mollist.end(), mollist1.begin(), mollist1.end());
    mollist.insert(mollist.end(), mollist2.begin(), mollist2.end());

    for (auto i = mollist.begin(); i != mollist.end(); ++i) {
        for (auto j = i; ++j != mollist.end();) {
            if ((i->get().id == id1 && j->get().id == id2) || (i->get().id == id2 && j->get().id == id1)) {
                double r = std::sqrt(spc.geo.sqdist(i->get().cm, j->get().cm));
                hist(r)++;
            }
        }
    }
}
MoleculeRDF::MoleculeRDF(const json &j, Space &spc) : PairFunctionBase(j), spc(spc) {
    name = "molrdf";
    id1 = findMoleculeByName(name1).id();
    id2 = findMoleculeByName(name2).id();
    assert(id1 >= 0 && id2 >= 0);
}

void AtomDipDipCorr::_sample() {
    V += spc.geo.getVolume(dim);
    auto active = spc.activeParticles();
    for (auto i = active.begin(); i != active.end(); ++i) {
        for (auto j = i; ++j != active.end();) {
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                Point rvec = spc.geo.vdist(i->pos, j->pos);
                if (slicedir.sum() > 0) {
                    if (rvec.cwiseProduct(slicedir.cast<double>()).norm() < thickness) {
                        if(i->hasExtension() && j->hasExtension()) {
                            double dipdip = i->getExt().mu.dot(j->getExt().mu);
                            double r1 = rvec.norm();
                            hist2(r1) += dipdip;
                            hist(r1)++; // get g(r) for free
                        }
                    }
                } else {
                    if(i->hasExtension() && j->hasExtension()) {
                        double dipdip = i->getExt().mu.dot(j->getExt().mu);
                        double r1 = rvec.norm();
                        hist2(r1) += dipdip;
                        hist(r1)++; // get g(r) for free
                    }
                }
            }
        }
    }
}
AtomDipDipCorr::AtomDipDipCorr(const json &j, Space &spc) : PairAngleFunctionBase(j), spc(spc) {
    name = "atomdipdipcorr";
    id1 = findAtomByName(name1).id();
    id2 = findAtomByName(name2).id();
}

// =============== XTCtraj ===============

XTCtraj::XTCtraj(const json &j, Space &s) : filter([](Particle &) { return true; }), spc(s) {
    from_json(j);
    name = "xtcfile";
    assert(filter); // filter must be callable
}

void XTCtraj::_to_json(json &j) const {
    j["file"] = writer->filename;
    if (not names.empty())
        j["molecules"] = names;
}

void XTCtraj::_from_json(const json &j) {
    auto file = MPI::prefix + j.at("file").get<std::string>();
    writer = std::make_shared<XTCWriter>(file);

    // By default, *all* active and inactive groups are saved,
    // but here allow for a user defined list of molecule ids
    names = j.value("molecules", std::vector<std::string>());
    if (not names.empty()) {
        molids = Faunus::names2ids(Faunus::molecules, names); // molecule types to save
        if (not molids.empty())
            filter = [&](Particle &i) {
                for (auto &g : spc.groups)     // loop over all active and inactive groups
                    if (g.contains(i, true)) { // does group contain particle?
                        if (std::find(molids.begin(), molids.end(), g.id) != molids.end())
                            return true;
                        else
                            return false;
                    }
                return false;
            };
    };
}

void XTCtraj::_sample() {
    // On some gcc/clang and certain ubuntu/macos combinations,
    // the ranges::view::filter(rng,unaryp) clears the `filter` function.
    // Using the ranges piping seem to solve the issue.
    assert(filter);
    auto coordinates = spc.p | ranges::cpp20::views::filter(filter) |
                       ranges::cpp20::views::transform([](auto &particle) -> Point & { return particle.pos; });
    assert(filter);
    writer->writeNext(spc.geo.getLength(), coordinates.begin(), coordinates.end());
}

// =============== MultipoleDistribution ===============

double MultipoleDistribution::g2g(const MultipoleDistribution::Tgroup &g1, const MultipoleDistribution::Tgroup &g2) {
    double u = 0;
    for (auto &i : g1)
        for (auto &j : g2)
            u += i.charge * j.charge / spc.geo.vdist(i.pos, j.pos).norm();
    return u;
}

/**
 * @note `fmt` is currently included w. spdlog but has been accepted into c++20.
 */
void MultipoleDistribution::save() const {
    if (number_of_samples > 0) {
        if (std::ofstream file(MPI::prefix + filename.c_str()); file) {
            file << "# Multipolar energies (kT/lB)\n"
                 << fmt::format("# {:>8}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n", "R", "exact", "tot", "ii", "id",
                                "dd", "iq", "mucorr");
            for (auto [r, u] : m) {
                double u_tot = u.ii.avg() + u.id.avg() + u.dd.avg() + u.iq.avg();
                file << fmt::format("{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}\n", r * dr,
                                    u.exact, u_tot, u.ii, u.id, u.dd, u.iq, u.mucorr);
            }
        }
    }
}
void MultipoleDistribution::_sample() {
    for (auto &gi : spc.findMolecules(ids[0]))     // find active molecules
        for (auto &gj : spc.findMolecules(ids[1])) // find active molecules
            if (gi != gj) {
                auto a = Faunus::toMultipole(gi, spc.geo.getBoundaryFunc());
                auto b = Faunus::toMultipole(gj, spc.geo.getBoundaryFunc());
                Point R = spc.geo.vdist(gi.cm, gj.cm);
                auto &d = m[to_bin(R.norm(), dr)];
                d.exact += g2g(gi, gj);
                d.ii += a.charge * b.charge / R.norm();
                d.id += q2mu(a.charge * b.getExt().mulen, b.getExt().mu, b.charge * a.getExt().mulen, a.getExt().mu, R);
                d.dd += mu2mu(a.getExt().mu, b.getExt().mu, a.getExt().mulen * b.getExt().mulen, R);
                d.iq += q2quad(a.charge, b.getExt().Q, b.charge, a.getExt().Q, R);
                d.mucorr += a.getExt().mu.dot(b.getExt().mu);
            }
}

void MultipoleDistribution::_to_json(json &j) const { j = {{"molecules", names}, {"file", filename}, {"dr", dr}}; }

MultipoleDistribution::MultipoleDistribution(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "Multipole Distribution";
    dr = j.at("dr").get<double>();
    filename = j.at("file").get<std::string>();
    names = j.at("molecules").get<decltype(names)>(); // molecule names
    ids = names2ids(molecules, names);                // names --> molids
    if (ids.size() != 2)
        throw std::runtime_error("specify exactly two molecules");
}

void MultipoleDistribution::_to_disk() { save(); }

// =============== AtomInertia ===============

void AtomInertia::_to_json(json &j) const {
    j["index"] = index; // atom id
}
Point AtomInertia::compute() {
    auto slice = spc.findAtoms(index);
    auto cm = Geometry::massCenter(slice.begin(), slice.end(), spc.geo.getBoundaryFunc());
    auto I = Geometry::inertia(slice.begin(), slice.end(), cm, spc.geo.getBoundaryFunc());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(I);
    return esf.eigenvalues();
}
void AtomInertia::_sample() {
    if (file)
        file << getNumberOfSteps() << " " << compute().transpose() << "\n";
}
AtomInertia::AtomInertia(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "Atomic Inertia Eigenvalues";
    filename = MPI::prefix + j.at("file").get<std::string>();
    file.open(filename); // output file
    index = j.at("index").get<size_t>(); // atom id
}
void AtomInertia::_to_disk() {
    if (file)
        file.flush(); // empty buffer
}

// =============== InertiaTensor ===============

void InertiaTensor::_to_json(json &j) const {
    j["indexes"] = indexes; // range of indexes within the group
    j["index"] = index; // group index
}
InertiaTensor::Data InertiaTensor::compute() {
    Space::Tgroup g(spc.groups[index].begin()+indexes[0], spc.groups[index].begin()+indexes[1]+1);
    InertiaTensor::Data d;
    auto I = Geometry::inertia(g.begin(), g.end(), spc.groups[index].cm, spc.geo.getBoundaryFunc());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(I);
    d.eivals = esf.eigenvalues();
    std::ptrdiff_t i_eival;
    d.eivals.minCoeff(&i_eival);
    d.eivec = esf.eigenvectors().col(i_eival).real(); // eigenvector corresponding to the smallest eigenvalue
    return d;
}
void InertiaTensor::_sample() {
    InertiaTensor::Data d = compute();
    if (file)
        file << getNumberOfSteps() << " " << d.eivals.transpose() << " " << d.eivec.transpose() << "\n";
}
InertiaTensor::InertiaTensor(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "Inertia Tensor";
    filename = MPI::prefix + j.at("file").get<std::string>();
    file.open(filename); // output file
    index = j.at("index").get<size_t>(); // group index
    indexes = j.value("indexes", std::vector<size_t>({0, spc.groups[index].size()})); // whole molecule by default
}
void InertiaTensor::_to_disk() {
    if (file)
        file.flush(); // empty buffer
}

// =============== MultipoleMoments ===============

void MultipoleMoments::_to_json(json &j) const {
    std::ostringstream o;
    std::string atom1 = atoms[(spc.groups[index].begin()+indexes[0])->id].name;
    std::string atom2 = atoms[(spc.groups[index].begin()+indexes[1])->id].name;
    o << atom1 << " " << indexes[0] << " – " << atom2 << " " << indexes[1];
    j["particles"] = o.str(); // range of particles within the group
    j["molecule"] = molecules[spc.groups[index].id].name; // group name
}
MultipoleMoments::Data MultipoleMoments::compute() {
    Space::Tgroup g(spc.groups[index].begin()+indexes[0], spc.groups[index].begin()+indexes[1]+1);
    auto cm = spc.groups[index].cm;
    if (not mol_cm)
        cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc());
    MultipoleMoments::Data d;
    Tensor S; // quadrupole tensor
    S.setZero();
    for (auto &i : g) {
        Point t = i.pos - cm;
        spc.geo.boundary(t);
        d.q += i.charge;
        d.mu += i.charge * t;
        S += i.charge * ( 3 * t * t.transpose() - Eigen::Matrix<double, 3, 3>::Identity() * t.squaredNorm() );
    }
    S = 0.5 * S;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(S);
    d.eivals = esf.eigenvalues();
    std::ptrdiff_t i_eival;
    d.eivals.minCoeff(&i_eival);
    d.eivec = esf.eigenvectors().col(i_eival).real(); // eigenvector corresponding to the smallest eigenvalue
    d.center = cm;
    return d;
}
void MultipoleMoments::_sample() {
    MultipoleMoments::Data d = compute();
    if (file)
        file << getNumberOfSteps() << " " << d.q << " " << d.mu.transpose() << " " << d.center.transpose() << " "
             << d.eivals.transpose() << " " << d.eivec.transpose() << "\n";
}
MultipoleMoments::MultipoleMoments(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "Multipole Moments";
    filename = MPI::prefix + j.at("file").get<std::string>();
    file.open(filename); // output file
    index = j.at("index").get<size_t>(); // group index
    indexes = j.value("indexes", std::vector<size_t>({0, spc.groups[index].size()-1})); // whole molecule by default
    mol_cm = j.value("mol_cm", true); // use the mass center of the whole molecule
}
void MultipoleMoments::_to_disk() {
    if (file)
        file.flush(); // empty buffer
}

// =============== PolymerShape ===============

void PolymerShape::_to_json(json &j) const {
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
        throw std::runtime_error(name + ": tensor output `file` cannot be used with multiple molecules");
    }
    for (const auto &group : molecules) {
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
                const auto &t = gyration_tensor;
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
        gyration_radius_histogram.stream_decorator = [](auto &stream, auto Rg, auto observations) {
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

PolymerShape::PolymerShape(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "Polymer Shape";
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

void AtomProfile::_from_json(const json &j) {
    ref = j.value("origo", Point(0, 0, 0));
    dir = j.value("dir", dir);
    file = j.at("file").get<std::string>();
    names = j.at("atoms").get<decltype(names)>();              // atom names
    auto vec_of_ids = names2ids(Faunus::atoms, names);         // names --> molids
    ids = std::set<int>(vec_of_ids.begin(), vec_of_ids.end()); // copy vector to set
    dr = j.value("dr", 0.1);
    tbl.setResolution(dr, 0);
    count_charge = j.value("charge", false);
    if (j.count("atomcom") == 1) {
        atom_com = j.at("atomcom");
        id_com = findAtomByName(atom_com).id();
    }
}
void AtomProfile::_to_json(json &j) const {
    j = {{"origo", ref}, {"dir", dir}, {"atoms", names}, {"file", file}, {"dr", dr}, {"charge", count_charge}};
}
void AtomProfile::_sample() {
    Group<Particle> all(spc.p.begin(), spc.p.end());
    if (id_com >= 0) { // calc. mass center of selected atoms
        auto slice = all.find_id(id_com);
        ref = Geometry::massCenter(slice.begin(), slice.end(), spc.geo.getBoundaryFunc());
    }
    for (auto &g : spc.groups)
        for (auto &p : g)
            if (ids.count(p.id) > 0) {
                Point rvec = spc.geo.vdist(p.pos, ref);
                double r = rvec.cwiseProduct(dir.cast<double>()).norm();
                if (count_charge)
                    tbl(r) += p.charge; // count charges
                else
                    tbl(r) += 1; // count atoms
            }
}
AtomProfile::AtomProfile(const json &j, Space &spc) : spc(spc) {
    name = "atomprofile";
    from_json(j);
}
void AtomProfile::_to_disk() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        tbl.stream_decorator = [&](std::ostream &o, double r, double N) {
            double Vr = 1;
            int dim = dir.sum();
            switch (dim) {
            case 3:
                Vr = 4 * pc::pi * std::pow(r, 2) * dr;
                break;
            case 2:
                Vr = 2 * pc::pi * r * dr;
                break;
            case 1:
                Vr = dr;
                break;
            default:
                throw std::runtime_error("bad dimension");
            }
            if (Vr < dr) // take care of the case where Vr=0
                Vr = dr; // then the volume element is simply dr

            N = N / double(number_of_samples);                            // average number of particles/charges
            o << r << " " << N << " " << N / Vr * 1e27 / pc::Nav << "\n"; // ... and molar concentration
        };
        f << "# r N rho/M\n" << tbl;
    }
}
void SlicedDensity::_from_json(const json &j) {
    file = j.at("file").get<std::string>();
    names = j.at("atoms").get<decltype(names)>(); // molecule names
    ids = names2ids(atoms, names);                // names --> molids
    dz = j.value("dz", 0.1);
    if (j.count("atomcom") == 1) {
        atom_com = j.at("atomcom");
        id_com = findAtomByName(atom_com).id();
    }
    N.setResolution(dz);
}
void SlicedDensity::_to_json(json &j) const {
    j = {{"atoms", names}, {"file", file}, {"dz", dz}, {"atomcom", atom_com}};
}

void SlicedDensity::_sample() {
    Group<Particle> all(spc.p.begin(), spc.p.end());
    double zcm = 0;
    if (id_com >= 0) { // calc. mass center of selected atoms
        auto slice = all.find_id(id_com);
        zcm = Geometry::massCenter(slice.begin(), slice.end(), spc.geo.getBoundaryFunc()).z();
    }
    // count atoms in slices
    for (auto &g : spc.groups) // loop over all groups
        for (auto &i : g)      // loop over active particles
            if (std::find(ids.begin(), ids.end(), i.id) not_eq ids.end())
                N(i.pos.z() - zcm)++;
}
SlicedDensity::SlicedDensity(const json &j, Space &spc) : spc(spc) {
    name = "sliceddensity";
    from_json(j);
}
void SlicedDensity::_to_disk() {
    if (std::ofstream f(MPI::prefix + file); f and number_of_samples > 0) {
        f << "# z rho/M\n";
        Point L = spc.geo.getLength();
        double halfz = 0.5 * L.z();
        double volume = L.x() * L.y() * dz;
        for (double z = -halfz; z <= halfz; z += dz)
            f << z << " " << N(z) / volume / number_of_samples * 1e27 / pc::Nav << "\n";
    }
}
void ChargeFluctuations::_sample() {
    for (auto &g : spc.findMolecules(mol_iter->id(), Space::ACTIVE)) {
        size_t cnt = 0;
        for (auto &p : g) {
            idcnt[cnt][p.id]++;
            charge[cnt] += p.charge;
            cnt++;
        }
    }
}
void ChargeFluctuations::_to_json(json &j) const {
    std::vector<std::string> mainname; // main name of atom with fluctuating charge
    std::vector<double> qavg;          // average charge
    std::vector<double> qstdev;        // standard deviation of the charge
    for (size_t i = 0; i < idcnt.size(); ++i) {
        qavg.push_back(charge.at(i).avg());
        qstdev.push_back(charge.at(i).stdev());
        // we look for the id that was sampled most often
        auto id_max = std::max_element(std::begin(idcnt.at(i)), std::end(idcnt.at(i)),
                                       [](auto &p1, auto &p2) { return p1.second < p2.second; });
        mainname.push_back(atoms.at(id_max->first).name);
    }
    if (verbose)
        j = {{"dominant atoms", mainname}, {"<q>", qavg}, {"std", qstdev}};
    j["molecule"] = mol_iter->name;
    if (not file.empty())
        j["pqrfile"] = file;
}
void ChargeFluctuations::_to_disk() {
    if (not file.empty()) {
        auto molecules = spc.findMolecules(mol_iter->id(), Space::ALL);
        if (not ranges::cpp20::empty(molecules)) {
            auto &g = *molecules.begin();
            ParticleVector pvec;        // temporary particle vector
            pvec.reserve(g.capacity()); // allocate required memory already now
            size_t cnt = 0;
            for (auto p = g.begin(); p < g.trueend(); ++p) {
                // we look for the id that was sampled most often
                auto id_max = std::max_element(std::begin(idcnt.at(cnt)), std::end(idcnt.at(cnt)),
                                               [](auto &p1, auto &p2) { return p1.second < p2.second; });
                pvec.push_back(Faunus::atoms.at(id_max->first));
                pvec.back().charge = charge.at(cnt).avg();
                pvec.back().pos = p->pos - g.cm;
                spc.geo.boundary(pvec.back().pos);
                cnt++;
            }
            FormatPQR::save(MPI::prefix + file, pvec, spc.geo.getLength());
        }
    }
}

/**
 * @todo replace `mol_iter` with simple molid integer
 */
ChargeFluctuations::ChargeFluctuations(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "chargefluctuations";
    file = j.value("pqrfile", ""s);
    verbose = j.value("verbose", true);
    const std::string molname = j.at("molecule"); // molecule name
    auto molecule = findMoleculeByName(molname); // throws if not found
    if (molecule.atomic) {
        throw ConfigurationError("only molecular groups allowed");
    }
    mol_iter = Faunus::findName(Faunus::molecules, molname);
    idcnt.resize(molecule.atoms.size());
    charge.resize(molecule.atoms.size());
}
void Multipole::_sample() {
    for (const auto &group : spc.groups) {            // loop over all groups molecules
        if (group.isMolecular() and !group.empty()) { // only active, molecular groups
            auto particle = Faunus::toMultipole(group, spc.geo.getBoundaryFunc());
            auto &average = average_moments[group.id];
            average.charge += particle.charge;
            average.dipole_moment += particle.getExt().mulen;
            average.charge_squared += particle.charge * particle.charge;
            average.dipole_moment_squared += particle.getExt().mulen * particle.getExt().mulen;
        }
    }
}
void Multipole::_to_json(json &j) const {
    auto &molecules_json = j["molecules"];
    for (auto &[molid, average] : average_moments) {
        auto &molecule_name = Faunus::molecules[molid].name;
        molecules_json[molecule_name] = {{"Z", average.charge.avg()},
                                         {"Z2", average.charge_squared.avg()},
                                         {"C", average.charge_squared.avg() - std::pow(average.charge.avg(), 2)},
                                         {u8::mu, average.dipole_moment.avg()},
                                         {u8::mu + u8::squared, average.dipole_moment_squared.avg()}};
    }
}
Multipole::Multipole(const json &j, const Space &spc) : spc(spc) {
    from_json(j);
    name = "multipole";
}
void ScatteringFunction::_sample() {
    p.clear();
    for (int id : ids) { // loop over molecule names
        auto groups = spc.findMolecules(id);
        for (auto &g : groups) // loop over groups
            if (use_com && !g.atomic)
                p.push_back(g.cm);
            else
                for (auto &i : g) // loop over particle index in group
                    p.push_back(i.pos);
    }

    // zero-padded suffix to use with `save_after_sample`
    std::string suffix = fmt::format("{:07d}", number_of_samples);
    switch (scheme) {
    case DEBYE:
        debye->sample(p, spc.geo.getVolume());
        if (save_after_sample)
            IO::write(filename + "." + suffix, debye->getIntensity());
        break;
    case EXPLICIT_PBC:
        explicit_average_pbc->sample(p, spc.geo.getLength());
        if (save_after_sample)
            IO::write(filename + "." + suffix, explicit_average_pbc->getSampling());
        break;
    case EXPLICIT_IPBC:
        explicit_average_ipbc->sample(p, spc.geo.getLength());
        if (save_after_sample)
            IO::write(filename + "." + suffix, explicit_average_ipbc->getSampling());
        break;
    }
}
void ScatteringFunction::_to_json(json &j) const {
    j = {{"molecules", names}, {"com", use_com}};
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

ScatteringFunction::ScatteringFunction(const json &j, Space &spc) try : spc(spc) {
    from_json(j);
    name = "scatter";
    use_com = j.value("com", true);
    save_after_sample = j.value("stepsave", false);   // save to disk for each sample
    filename = j.at("file").get<std::string>();
    names = j.at("molecules").get<decltype(names)>(); // molecule names
    ids = names2ids(molecules, names);                // names --> molids

    if (std::string scheme_str = j.value("scheme", "explicit"); scheme_str == "debye") {
        scheme = DEBYE;
        debye = std::make_shared<Scatter::DebyeFormula<Tformfactor>>(j);
        // todo: add warning if used on PBC system
        // if (spc.geo.geometry->boundary_conditions.x() == Geometry::PERIODIC)
        //    faunus_logger->warn("{0}: '{1}' scheme used on PBC system; consider 'explicit' instead", name,
        //    scheme_str);
    } else if (scheme_str == "explicit") {
        bool ipbc = j.value("ipbc", false);
        int pmax = j.value("pmax", 15);
        if (ipbc) {
            scheme = EXPLICIT_IPBC;
            explicit_average_ipbc = std::make_shared<Scatter::StructureFactorIPBC<>>(pmax);
        } else {
            scheme = EXPLICIT_PBC;
            explicit_average_pbc = std::make_shared<Scatter::StructureFactorPBC<>>(pmax);
        }
        // todo: add warning if used a non-cubic system
    } else
        throw std::runtime_error("unknown scheme");
} catch (std::exception &e) {
    throw std::runtime_error("scatter: "s + e.what());
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
    try {
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
            output_stream = IO::openCompressedOutputStream(filename, true); // throws if error
            *output_stream << "# steps dL/Å du/kT <force>/kT/Å\n"s;
        }
    } catch (std::exception& e) {
        throw ConfigurationError("{}", e.what());
    }
}
void VirtualTranslate::_sample() {
    if (std::fabs(perturbation_distance) > 0.0) {
        if (auto mollist = spc.findMolecules(molid, Space::ACTIVE); !ranges::cpp20::empty(mollist)) {
            if (ranges::distance(mollist.begin(), mollist.end()) > 1) {
                throw std::runtime_error(name + ": exactly ONE active molecule expected");
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
    if (output_stream) { // file to disk?
        const auto mean_force = -meanFreeEnergy() / perturbation_distance;
        *output_stream << fmt::format("{:d} {:.3f} {:.10f} {:.6f}\n", getNumberOfSteps(), perturbation_distance,
                                      energy_change, mean_force);
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
    : PerturbationAnalysisBase(pot, spc, j.value("file", ""s)) {
    name = "virtualtranslate";
    change.groups.resize(1);
    change.groups.front().internal = false;
    from_json(j);
}

SpaceTrajectory::SpaceTrajectory(const json &j, Space::Tgvec &groups) : groups(groups) {
    from_json(j);
    name = "space trajectory";
    filename = j.at("file");
    if (useCompression())
        stream = std::make_unique<zstr::ofstream>(MPI::prefix + filename, std::ios::binary);
    else
        stream = std::make_unique<std::ofstream>(MPI::prefix + filename, std::ios::binary);

    if (stream != nullptr && *stream)
        archive = std::make_unique<cereal::BinaryOutputArchive>(*stream);

    if (not archive)
        throw std::runtime_error("error creating "s + filename);
}

bool SpaceTrajectory::useCompression() const {
    assert(!filename.empty());
    std::string suffix = filename.substr(filename.find_last_of(".") + 1);
    if (suffix == "ztraj")
        return true;
    else if (suffix == "traj")
        return false;
    else
        throw std::runtime_error("Trajectory file suffix must be `.traj` or `.ztraj`");
}

void SpaceTrajectory::_sample() {
    assert(archive);
    for (auto &group : groups) {
        (*archive)(group);
    }
}

void SpaceTrajectory::_to_json(json &j) const { j = {{"file", filename}}; }

void SpaceTrajectory::_to_disk() {
    assert(*stream);
    stream->flush();
}
} // namespace Faunus::Analysis
