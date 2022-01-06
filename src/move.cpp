#include <doctest/doctest.h>
#include "core.h"
#include "move.h"
#include "speciation.h"
#include "clustermove.h"
#include "chainmove.h"
#include "forcemove.h"
#include "montecarlo.h"
#include "regions.h"
#include "aux/iteratorsupport.h"
#include "aux/eigensupport.h"
#include "spdlog/spdlog.h"
#include <range/v3/view/counted.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/algorithm/count.hpp>
#include <range/v3/view/transform.hpp>

namespace Faunus::Move {

Random MoveBase::slump; // static instance of Random (shared for all moves)

void MoveBase::from_json(const json& j) {
    if (const auto it = j.find("repeat"); it != j.end()) {
        if (it->is_number()) {
            repeat = it->get<int>();
        } else if (it->is_string() && it->get<std::string>() == "N") {
            repeat = -1;
        } else {
            throw std::runtime_error("invalid 'repeat'");
        }
    }
    sweep_interval = j.value("nstep", 1); // Non-stochastic moves are defined with `repeat=0`...
    if (sweep_interval > 1) {             // ...the move is then instead run at a fixed sweep interval
        repeat = 0;
    }
    _from_json(j);
    if (repeat < 0) {
        repeat = 0;
    }
}

void MoveBase::to_json(json& j) const {
    _to_json(j);
    if (timer_move.result() > 0.01) { // only print if more than 1% of the time
        j["relative time (without energy calc)"] = timer_move.result();
    }
    if (timer.result() > 0.01) { // only print if more than 1% of the time
        j["relative time"] = timer.result();
    }
    j["acceptance"] = double(number_of_accepted_moves) / number_of_attempted_moves;
    j["repeat"] = repeat;
    j["stochastic"] = isStochastic();
    j["moves"] = number_of_attempted_moves;
    if (!cite.empty()) {
        j["cite"] = cite;
    }
    roundJSON(j, 3);
}

void MoveBase::move(Change& change) {
    timer.start();
    timer_move.start();
    number_of_attempted_moves++;
    change.clear();
    _move(change);
    if (change.empty()) {
        timer.stop();
    }
    timer_move.stop();
}

void MoveBase::accept(Change& change) {
    number_of_accepted_moves++;
    _accept(change);
    timer.stop();
}

void MoveBase::reject(Change& change) {
    number_of_rejected_moves++;
    _reject(change);
    timer.stop();
}

/**
 * This can be used to introduce an extra energy change that will be added to the
 * total trial energy in the Metropolis acceptance criterion. Typically used
 * to include bias not contained in Energy::Hamiltonian.
 *
 * @param change Change made to the system
 * @param old_energy Energy from hamiltonian before the change (kT)
 * @param new_energy Energy from hamiltonian after the change (kT)
 * @return Energy due to custom bias from the particular move (kT)
 */
double MoveBase::bias([[maybe_unused]] Change& change, [[maybe_unused]] double old_energy,
                      [[maybe_unused]] double new_energy) {
    return 0.0;
}

void MoveBase::_accept([[maybe_unused]] Change& change) {}

void MoveBase::_reject([[maybe_unused]] Change& change) {}

MoveBase::MoveBase(Space& spc, std::string_view name, std::string_view cite) : cite(cite), name(name), spc(spc) {}

void MoveBase::setRepeat(const int new_repeat) { repeat = new_repeat; }
bool MoveBase::isStochastic() const { return repeat != 0; }

void from_json(const json &j, MoveBase &move) { move.from_json(j); }

void to_json(json& j, const MoveBase& move) { move.to_json(j[move.getName()]); }

const std::string& MoveBase::getName() const {
    assert(!name.empty());
    return name;
}

// -----------------------------------

ReplayMove::ReplayMove(Space& spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {}

ReplayMove::ReplayMove(Space& spc) : ReplayMove(spc, "replay", "") {}

void ReplayMove::_to_json(json& j) const { j["file"] = reader->filename; }

void ReplayMove::_from_json(const json& j) { reader = std::make_unique<XTCReader>(j.at("file").get<std::string>()); }

void ReplayMove::_move(Change &change) {
    assert(reader);
    if (!end_of_trajectory) {
        if (reader->read(frame.step, frame.timestamp, frame.box, spc.positions().begin(), spc.positions().end())) {
            spc.geometry.setLength(frame.box);
            change.everything = true;
        } else {
            // nothing to do, simulation shall stop
            end_of_trajectory = true;
            mcloop_logger->warn("No more frames to read from {}. Running on empty.", reader->filename);
        }
    }
}

double ReplayMove::bias(Change&, double, double) {
    return force_accept; // always accept
}

void AtomicTranslateRotate::_to_json(json& j) const {
    j = {{"dir", directions},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(mean_square_displacement.avg())},
         {"molecule", molecule_name}};
    roundJSON(j, 3);
}

void AtomicTranslateRotate::_from_json(const json& j) {
    assert(!molecules.empty());
    molecule_name = j.at("molecule");
    const auto molecule = findMoleculeByName(molecule_name);
    molid = molecule.id();
    if (molecule.rigid) {
        faunus_logger->warn("structure of rigid molecule {} may be disturbed by {}", molecule_name, name);
    }
    directions = j.value("dir", Point(1, 1, 1));
    if (repeat < 0) {
        auto mollist = spc.findMolecules(molid, Space::Selection::ALL);
        repeat = std::distance(mollist.begin(), mollist.end()); // repeat for each molecule...
        if (repeat > 0) {
            repeat = repeat * mollist.front().size(); // ...and for each atom
        }
    }
    energy_resolution = j.value("energy_resolution", 0.0);
}

void AtomicTranslateRotate::translateParticle(ParticleVector::iterator particle, double displacement) {
    const auto old_position = particle->pos; // backup old position
    particle->pos += randomUnitVector(slump, directions) * displacement * slump();
    spc.geometry.boundary(particle->pos);
    latest_displacement_squared = spc.geometry.sqdist(old_position, particle->pos); // square displacement

    auto& group = spc.groups.at(cdata.group_index);
    if (group.isMolecular()) {
        group.mass_center =
            Geometry::massCenter(group.begin(), group.end(), spc.geometry.getBoundaryFunc(), -group.mass_center);
        checkMassCenter(group);
    }
}

void AtomicTranslateRotate::checkMassCenter(Space::GroupType& group) const {
    const auto allowed_threshold = 1e-6;
    const auto old_mass_center = group.mass_center;
    group.translate(-old_mass_center, spc.geometry.getBoundaryFunc()); // translate to origin
    const auto should_be_zero = spc.geometry.sqdist({0, 0, 0}, Geometry::massCenter(group.begin(), group.end()));
    if (should_be_zero > allowed_threshold) {
        faunus_logger->error("{}: error calculating mass center for {}", name, group.traits().name);
        groupToDisk(group);
        throw std::runtime_error("molecule likely too large for periodic boundaries; increase box size?");
    }
    group.translate(old_mass_center, spc.geometry.getBoundaryFunc());
}

void AtomicTranslateRotate::groupToDisk(const Space::GroupType& group) const {
    if (auto stream = std::ofstream("mass-center-failure.pqr"); stream) {
        const auto group_iter = spc.groups.cbegin() + spc.getGroupIndex(group);
        auto groups = ranges::cpp20::views::counted(group_iter, 1); // slice out single group
        PQRWriter().save(stream, groups, spc.geometry.getLength());
    }
}

void AtomicTranslateRotate::_move(Change& change) {
    if (auto particle = randomAtom(); particle != spc.particles.end()) {
        latest_particle = particle;
        const auto translational_displacement = particle->traits().dp;
        const auto rotational_displacement = particle->traits().dprot;

        if (translational_displacement > 0.0) { // translate
            translateParticle(particle, translational_displacement);
        }

        if (rotational_displacement > 0.0) { // rotate
            const auto random_unit_vector = Faunus::randomUnitVector(slump);
            const auto angle = rotational_displacement * (slump() - 0.5);
            Eigen::Quaterniond quaternion(Eigen::AngleAxisd(angle, random_unit_vector));
            particle->rotate(quaternion, quaternion.toRotationMatrix());
        }

        if (translational_displacement > 0.0 or rotational_displacement > 0.0) {
            change.groups.push_back(cdata); // add to list of moved groups
        }
    } else {
        latest_particle = spc.particles.end();
        latest_displacement_squared = 0.0; // no particle found --> no movement
    }
}

void AtomicTranslateRotate::_accept(Change&) {
    mean_square_displacement += latest_displacement_squared;
    sampleEnergyHistogram();
}

void AtomicTranslateRotate::_reject(Change&) { mean_square_displacement += 0; }

AtomicTranslateRotate::AtomicTranslateRotate(Space& spc, const Energy::Hamiltonian& hamiltonian, std::string name,
                                             std::string cite)
    : MoveBase(spc, name, cite), hamiltonian(hamiltonian) {
    repeat = -1; // meaning repeat N times
    cdata.relative_atom_indices.resize(1);
    cdata.internal = true;
}

AtomicTranslateRotate::AtomicTranslateRotate(Space& spc, const Energy::Hamiltonian& hamiltonian)
    : AtomicTranslateRotate(spc, hamiltonian, "transrot", "") {}

/**
 * For atomic groups, select `ALL` since these may be partially filled and thereby
 * appear inactive. Note also that only one instance of atomic molecules can exist.
 * For molecular groups, select only active ones.
 *
 * @return Iterator to particle to move; `end()` if nothing selected
 */
ParticleVector::iterator AtomicTranslateRotate::randomAtom() {
    assert(molid >= 0);
    auto particle = spc.particles.end(); // particle iterator
    auto selection = (Faunus::molecules[molid].atomic) ? Space::Selection::ALL : Space::Selection::ACTIVE;
    auto mollist = spc.findMolecules(molid, selection);
    if (auto group = slump.sample(mollist.begin(), mollist.end()); group != mollist.end()) { // random molecule
        if (not group->empty()) {
            particle = slump.sample(group->begin(), group->end());           // random particle
            cdata.group_index = Faunus::distance(spc.groups.begin(), group); // index of touched group
            cdata.relative_atom_indices[0] =
                std::distance(group->begin(), particle); // index of moved particle relative to group
        }
    }
    return particle;
}

/**
 * Here we access the Hamiltonian and sum all energy terms from the just performed MC move.
 * This is used to updated the histogram of energies, sampled for each individual particle type.
 */
void AtomicTranslateRotate::sampleEnergyHistogram() {
    if (energy_resolution > 0.0) {
        assert(latest_particle != spc.particles.end());
        const auto particle_energy =
            std::accumulate(hamiltonian.latestEnergies().begin(), hamiltonian.latestEnergies().end(), 0.0);
        auto& particle_histogram =
            energy_histogram.try_emplace(latest_particle->id, SparseHistogram(energy_resolution)).first->second;
        particle_histogram.add(particle_energy);
    }
}

void AtomicTranslateRotate::saveHistograms() {
    if (energy_resolution) {
        for (const auto& [atom_id, histogram] : energy_histogram) {
            const auto filename = fmt::format("energy-histogram-{}.dat", Faunus::atoms[atom_id].name);
            if (auto stream = std::ofstream(filename); stream) {
                stream << "# energy/kT observations\n" << histogram;
            }
        }
    }
}

AtomicTranslateRotate::~AtomicTranslateRotate() { saveHistograms(); }

std::unique_ptr<MoveBase> createMove(const std::string& name, const json& properties, Space& spc,
                                     Energy::Hamiltonian& hamiltonian) {
    try {
        std::unique_ptr<MoveBase> move;
        if (name == "moltransrot") {
            move = std::make_unique<TranslateRotate>(spc);
        } else if (name == "smartmoltransrot") {
            move = std::make_unique<SmartTranslateRotate>(spc);
        } else if (name == "smartmoltransrot2") {
            return std::make_unique<SmartTranslateRotate2>(spc, properties);
        } else if (name == "conformationswap") {
            move = std::make_unique<ConformationSwap>(spc);
        } else if (name == "transrot") {
            move = std::make_unique<AtomicTranslateRotate>(spc, hamiltonian);
        } else if (name == "pivot") {
            move = std::make_unique<PivotMove>(spc);
        } else if (name == "crankshaft") {
            move = std::make_unique<CrankshaftMove>(spc);
        } else if (name == "volume") {
            move = std::make_unique<VolumeMove>(spc);
        } else if (name == "charge") {
            move = std::make_unique<ChargeMove>(spc);
        } else if (name == "chargetransfer") {
            move = std::make_unique<ChargeTransfer>(spc);
        } else if (name == "rcmc") {
            move = std::make_unique<SpeciationMove>(spc);
        } else if (name == "quadrantjump") {
            move = std::make_unique<QuadrantJump>(spc);
        } else if (name == "cluster") {
            move = std::make_unique<Cluster>(spc);
        } else if (name == "replay") {
            move = std::make_unique<ReplayMove>(spc);
        } else if (name == "langevin_dynamics") {
            move = std::make_unique<LangevinDynamics>(spc, hamiltonian);
        } else if (name == "temper") {
#ifdef ENABLE_MPI
            move = std::make_unique<ParallelTempering>(spc, MPI::mpi);
            move->setRepeat(0); // zero weight moves are run at the end of each sweep
#else
            throw ConfigurationError("{} requires that Faunus is compiled with MPI", name);
#endif
        }
        if (!move) {
            throw ConfigurationError("unknown move '{}'", name);
        }
        move->from_json(properties);
        return move;
    } catch (std::exception& e) { throw ConfigurationError("error creating move -> {}", e.what()); }
}

void MoveCollection::addMove(std::shared_ptr<MoveBase>&& move) {
    if (!move) {
        throw std::runtime_error("invalid move");
    }
    moves.vec.emplace_back(move);
    repeats.push_back(static_cast<double>(move->repeat));
    distribution = std::discrete_distribution<unsigned int>(repeats.begin(), repeats.end());
    number_of_moves_per_sweep = static_cast<unsigned int>(std::accumulate(repeats.begin(), repeats.end(), 0.0));
}

MoveCollection::MoveCollection(const json& list_of_moves, Space& spc, Energy::Hamiltonian& hamiltonian) {
    assert(list_of_moves.is_array());
    for (const auto& j : list_of_moves) { // loop over move list
        const auto& [name, parameters] = jsonSingleItem(j);
        try {
            addMove(createMove(name, parameters, spc, hamiltonian));
        } catch (std::exception& e) {
            usageTip.pick(name);
            throw ConfigurationError("{}", e.what()).attachJson(j);
        }
    }
}

void to_json(json& j, const MoveCollection& propagator) { j = propagator.moves; }

void to_json(json& j, const SmartMonteCarlo& smc) {
    j["region"] = static_cast<json>(*smc.region);
    j["p"] = smc.outside_rejection_probability;
}

const BasePointerVector<MoveBase>& MoveCollection::getMoves() const { return moves; }

MoveCollection::move_iterator MoveCollection::sample() {
#ifdef ENABLE_MPI
    auto& random_engine = MPI::mpi.random.engine; // parallel processes (tempering) must be in sync
#else
    auto& random_engine = Move::MoveBase::slump.engine;
#endif
    if (!moves.empty()) {
        return moves.begin() + distribution(random_engine);
    }
    return moves.end();
}

#ifdef ENABLE_MPI

void ParallelTempering::_to_json(json& j) const {
    j = {{"replicas", mpi.world.size()},
         {"format", exchange_particles.getFormat()},
         {"partner_policy", partner->policy},
         {"volume_scale", volume_scaling_method}};
    auto& exchange_json = j["exchange"] = json::object();
    for (const auto& [pair, acceptance] : acceptance_map) {
        auto id = fmt::format("{} <-> {}", pair.first, pair.second);
        exchange_json[id] = {{"attempts", acceptance.size()}, {"acceptance", acceptance.avg()}};
    }
}

/**
 * Exchange groups sizes with partner MPI rank and Resize all groups to the exchanged values
 */
void ParallelTempering::exchangeGroupSizes(Space::GroupVector& groups, int partner_rank) {
    std::vector<size_t> sizes = groups | ranges::cpp20::views::transform(&Group::size) | ranges::to_vector;
    mpi.world.sendrecv_replace(sizes.begin(), sizes.end(), partner_rank, mpl::tag_t(0), partner_rank, mpl::tag_t(0));
    auto it = sizes.begin();
    ranges::cpp20::for_each(groups, [&it](Group& group) { group.resize(*it++); });
}

/**
 * This will exchange the states between two partner replicas and set the change object accordingy
 */
void ParallelTempering::exchangeState(Change& change) {
    if (MPI::exchangeVolume(mpi, *partner->rank, spc.geometry, volume_scaling_method)) {
        change.volume_change = true;
    }
    exchangeGroupSizes(spc.groups, *partner->rank);
    exchange_particles.replace(mpi.world, *partner->rank, spc.particles);
    spc.updateParticles(spc.particles.begin(), spc.particles.end(), spc.particles.begin());
    change.everything = true;
}

void ParallelTempering::_move(Change& change) {
    mpi.world.barrier(); // wait until all ranks reach here
    if (!MPI::checkRandomEngineState(mpi.world, MoveBase::slump)) {
        faunus_logger->error("Random numbers out of sync across MPI nodes. Do not use 'hardware' seed.");
        mpi.world.abort(1); // neighbor search *requires* that random engines are in sync
    }
    partner->generate(mpi.world, MoveBase::slump);
    if (partner->rank.has_value()) {
        exchangeState(change);
    }
}

/**
 * @param energy_change Energy change of current replica
 * @return Energy change in partner replica
 *
 * In the MC move, the energy of the current replica and
 * the bias (== energy change of the partner) are added together
 * to form the final trial energy for the tempering move.
 */
double ParallelTempering::exchangeEnergy(const double energy_change) {
    double partner_energy_change = 0.0;
    const auto tag = mpl::tag_t(0);
    mpi.world.sendrecv(energy_change, *partner->rank, tag, partner_energy_change, *partner->rank, tag);
    return partner_energy_change; // return partner energy change
}

/**
 * The bias() function takes the current old and new energy and exchanges
 * the resulting energy change with the partner replica. The change in the
 * replica is returned and will be added to the total trial energy in
 * the MetropolicMonteCarlo class.
 */
double ParallelTempering::bias([[maybe_unused]] Change& change, double uold, double unew) {
    return exchangeEnergy(unew - uold); // energy change in partner replica
}

void ParallelTempering::_accept([[maybe_unused]] Change& change) {
    acceptance_map[partner->getPair(mpi.world)] += 1.0;
}
void ParallelTempering::_reject([[maybe_unused]] Change& change) {
    acceptance_map[partner->getPair(mpi.world)] += 0.0;
}

void ParallelTempering::_from_json(const json& j) {
    exchange_particles.setFormat(j.value("format", MPI::ParticleBuffer::Format::XYZQI));
    partner = createMPIPartnerPolicy(j.value("partner_policy", MPI::PartnerPolicy::ODDEVEN));
    volume_scaling_method = j.value("volume_scale", Geometry::VolumeMethod::ISOTROPIC);
}

ParallelTempering::ParallelTempering(Space& spc, const MPI::Controller& mpi)
    : MoveBase(spc, "temper", "doi:10/b3vcw7"), mpi(mpi) {
    if (mpi.world.size() < 2) {
        throw std::runtime_error(name + " requires two or more MPI processes");
    }
    partner = MPI::createMPIPartnerPolicy(MPI::PartnerPolicy::ODDEVEN);
}

#endif

void VolumeMove::_to_json(json& j) const {
    if (number_of_attempted_moves > 0) {
        j = {{"dV", logarithmic_volume_displacement_factor},
             {"method", volume_scaling_method},
             {"⟨V⟩", mean_volume.avg()},
             {"√⟨ΔV²⟩", std::sqrt(mean_square_volume_change.avg())},
             {"∛√⟨ΔV²⟩", std::cbrt(std::sqrt(mean_square_volume_change.avg()))}};
        roundJSON(j, 3);
    }
}
void VolumeMove::_from_json(const json& j) {
    logarithmic_volume_displacement_factor = j.at("dV").get<double>();
    volume_scaling_method = j.value("method", Geometry::VolumeMethod::ISOTROPIC);
    if (volume_scaling_method == Geometry::VolumeMethod::INVALID) {
        throw ConfigurationError("invalid volume scaling method");
    }
}
void VolumeMove::_move(Change& change) {
    if (logarithmic_volume_displacement_factor > 0.0) {
        change.volume_change = true;
        change.everything = true;
        old_volume = spc.geometry.getVolume();
        new_volume = std::exp(std::log(old_volume) + (slump() - 0.5) * logarithmic_volume_displacement_factor);
        spc.scaleVolume(new_volume, volume_scaling_method);
    }
}
void VolumeMove::_accept([[maybe_unused]] Change& change) {
    mean_square_volume_change += std::pow(new_volume - old_volume, 2);
    mean_volume += new_volume;
    assert(std::fabs(spc.geometry.getVolume() - new_volume) < 1.0e-9);
}

VolumeMove::VolumeMove(Space& spc) : MoveBase(spc, "volume"s, ""s) { repeat = 1; }

void VolumeMove::_reject([[maybe_unused]] Change& change) {
    mean_square_volume_change += 0.0;
    mean_volume += old_volume;
    assert(std::fabs(spc.geometry.getVolume() - old_volume) < 1.0e-9);
}

// ------------------------------------------------

void ChargeMove::_to_json(json& j) const {
    using namespace u8;
    j = {{"index", atomIndex},
         {"dq", dq},
         {rootof + bracket(Delta + "q" + squared), std::sqrt(msqd.avg())},
         {cuberoot + rootof + bracket(Delta + "q" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
    roundJSON(j, 3);
}
void ChargeMove::_from_json(const json& j) {
    dq = j.at("dq").get<double>();
    atomIndex = j.at("index").get<int>();
    auto git = spc.findGroupContaining(spc.particles.at(atomIndex)); // group containing atomIndex
    cdata.group_index = std::distance(spc.groups.begin(), git);      // integer *index* of moved group
    cdata.relative_atom_indices[0] =
        std::distance(git->begin(), spc.particles.begin() + atomIndex); // index of particle rel. to group
}
void ChargeMove::_move(Change& change) {
    if (dq > 0) {
        auto& p = spc.particles.at(atomIndex); // refence to particle
        double qold = p.charge;
        p.charge += dq * (slump() - 0.5);
        deltaq = p.charge - qold;
        change.groups.push_back(cdata); // add to list of moved groups
    } else
        deltaq = 0;
}
void ChargeMove::_accept(Change&) { msqd += deltaq * deltaq; }
void ChargeMove::_reject(Change&) { msqd += 0; }

ChargeMove::ChargeMove(Space& spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = 1;
    cdata.internal = true;                 // the group is internally changed
    cdata.relative_atom_indices.resize(1); // we change exactly one atom
}

ChargeMove::ChargeMove(Space& spc) : ChargeMove(spc, "charge", "") {}

void ChargeTransfer::_to_json(json& j) const {
    using namespace u8;
    j = {{"dq", dq},
         {rootof + bracket(Delta + "q" + squared), std::sqrt(msqd.avg())},
         {cuberoot + rootof + bracket(Delta + "q" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
    roundJSON(j, 3);
}
void ChargeTransfer::_from_json(const json& j) {
    dq = j.at("dq").get<double>();
    mol1.molname = j.at("mol1");                                  // string containing name of molecule 1
    mol2.molname = j.at("mol2");                                  // string containing name of molecule 2
    mol1.id = findMoleculeByName(mol1.molname).id();              // group containing mol1.molname
    mol2.id = findMoleculeByName(mol2.molname).id();              // group containing mol2.molname
    mol1.molrange = j.at("molrange1").get<std::vector<double>>(); // vector containing lower and upper limit of
                                                                  // total charge of molecule 1
    mol2.molrange = j.at("molrange2").get<std::vector<double>>(); // vector containing lower and upper limit of
                                                                  // total charge of molecule 2
    mol1.min = j.at("min1").get<std::vector<double>>();           // vector containing lower limits of atomic charges in
                                                                  // molecule 1
    mol1.max = j.at("max1").get<std::vector<double>>();           // vector containing upper limits of atomic charges in
                                                                  // molecule 1
    mol2.min = j.at("min2").get<std::vector<double>>();           // vector containing lower limits of atomic charges in
                                                                  // molecule 2
    mol2.max = j.at("max2").get<std::vector<double>>();           // vector containing upper limits of atomic charges in
                                                                  // molecule 2

    if (repeat < 0) {
        auto v = spc.findMolecules(mol1.id);
        repeat = std::distance(v.begin(), v.end());
    }
    if (repeat < 0) {
        auto v = spc.findMolecules(mol2.id);
        repeat = std::distance(v.begin(), v.end());
    }

    if (mol1.min.size() != mol1.max.size()) {
        // checking so that mol1.min and mol1.max contains equal number of entries
        throw ConfigurationError("mol1.min and mol1.max need to have the same number of entries. "
                                 "mol1.min has {} and mol1.max has {} entries.",
                                 mol1.min.size(), mol1.max.size());
    }

    if (mol1.min.size() == 0 || mol1.max.size() == 0) {
        // checking so that mol1.min and mol1.max are not empty
        throw ConfigurationError("mol1.min and mol1.max both need to have nonzero number of entries. "
                                 "mol1.min has {} and mol1.max has {}  entries.",
                                 mol1.min.size(), mol1.max.size());
    }

    if (mol2.min.size() != mol2.max.size()) {
        // checking so that mol2.min and mol2.max contains equal number of entries
        throw ConfigurationError("mol2.min and mol2.max need to have the same number of entries. "
                                 "mol2.min has {} and mol2.max has {} entries.",
                                 mol2.min.size(), mol2.max.size());
    }

    if (mol2.min.size() == 0 || mol2.max.size() == 0) {
        // checking so that mol2.min and mol2.max are not empty
        throw ConfigurationError("mol2.min and mol2.max both need to have nonzero number of entries. "
                                 "mol2.min has {} and mol2.max has {} entries.",
                                 mol2.min.size(), mol2.max.size());
    }
}

void ChargeTransfer::_move(Change& change) {
    auto mollist1 = spc.findMolecules(mol1.id, Space::Selection::ACTIVE);
    auto mollist2 = spc.findMolecules(mol2.id, Space::Selection::ACTIVE);
    if ((not ranges::cpp20::empty(mollist1)) and (not ranges::cpp20::empty(mollist2))) {
        auto git1 = slump.sample(mollist1.begin(), mollist1.end()); // selecting a random molecule of type molecule1
        auto git2 = slump.sample(mollist2.begin(), mollist2.end()); // selecting a random molecule of type molecule2

        if (!git1->empty() && !git2->empty()) { // check that both molecule1 and molecule 2 exist

            if (dq > 0) {
                // change.chargeMove =
                //    true; // setting to true makes the self-energy being computed and added to the total energy
                mol1.numOfAtoms = Faunus::distance(git1->begin(), git1->end());
                mol2.numOfAtoms = Faunus::distance(git2->begin(), git2->end());

                mol1.ratio.clear(); // clearing vector containing ratio of atomic charge ranges and the charge range of
                                    // the whole molecule1
                mol2.ratio.clear(); // clearing vector containing ratio of atomic charge ranges and the charge range of
                                    // the whole molecule2

                for (i = 0; i < mol1.numOfAtoms; i++) {
                    mol1.ratio.push_back(
                        (mol1.max[i] - mol1.min[i]) /
                        (mol1.molrange[1] - mol1.molrange[0])); // calculating ratio of atom i in molecule 1
                }

                for (i = 0; i < mol2.numOfAtoms; i++) {
                    mol2.ratio.push_back(
                        (mol2.max[i] - mol2.min[i]) /
                        (mol2.molrange[1] - mol2.molrange[0])); // calculating ratio of atom i in molecule 2
                }

                mol1.charges = 0; // setting sum of all atomic charges in molecule1 to zero
                mol2.charges = 0; // setting sum of all atomic charges in molecule2 to zero
                deltaq = dq * (slump() - 0.5);
                mol1.changeQ.clear(); // clearing vector containing attempted charge moves on all atoms in molecule1
                mol2.changeQ.clear(); // clearing vector containing attempted charge moves on all atoms in molecule2
                mol1.cdata.group_index = Faunus::distance(spc.groups.begin(), git1);
                mol2.cdata.group_index = Faunus::distance(spc.groups.begin(), git2);

                for (i = 0; i < mol1.numOfAtoms; i++) {
                    auto p = git1->begin() + i; // object containing atom i in molecule1
                    mol1.changeQ.push_back(
                        deltaq * mol1.ratio[i]); // assigning attempted charge move of atom i in molecule1 to vector
                    // sumChanges1 += changeQ1[i];
                    mol1.charges +=
                        p->charge + mol1.changeQ[i]; // adding new attempted charge of atom i in molecule1 to sum
                }
                for (i = 0; i < mol2.numOfAtoms; i++) { // Doing the same as above loop but for molecule2
                    auto p = git2->begin() + i;
                    mol2.changeQ.push_back(-deltaq * mol2.ratio[i]);
                    // sumMoves2 += changeQ2[i];
                    mol2.charges += p->charge + mol2.changeQ[i];
                }

                // Torodial boundary conditions
                if (mol1.charges < mol1.molrange[0]) { // Checking if sum of new attempted atomic charges in molecule1
                                                       // will fall below lower limit in molrange1
                    sumTemp = 0;                       // resetting temporary sum of atomic charges
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        // temporary sum of charge moves attempted on all atoms in molecule1
                        sumTemp += p->charge - (2 * mol1.min[i] - (p->charge + mol1.changeQ[i]));
                        // new attempted charge of atom i in molecule1, obeying torodial boundary conditions
                        p->charge = 2 * mol1.min[i] - (p->charge + mol1.changeQ[i]);
                    }
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        // new attempted charge of atom i in molecule2, obeying torodial boundary conditions
                        p->charge += sumTemp * mol2.ratio[i];
                    }
                }

                else if (mol1.charges > mol1.molrange[1]) {
                    // same procedure as above if statement, but if sum of new atempted charges in molecule1 falls above
                    // upper limit in molrange1
                    sumTemp = 0;
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        sumTemp += p->charge - (2 * mol1.max[i] - (p->charge + mol1.changeQ[i]));
                        p->charge = 2 * mol1.max[i] - (p->charge + mol1.changeQ[i]);
                    }
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        p->charge += sumTemp * mol2.ratio[i];
                    }
                }

                else if (mol2.charges < mol2.molrange[0]) {
                    // same as first if statement, but with respect to molecule2 and its molrange
                    sumTemp = 0;
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        sumTemp += p->charge - (2 * mol2.min[i] - (p->charge + mol2.changeQ[i]));
                        p->charge = 2 * mol2.min[i] - (p->charge + mol2.changeQ[i]);
                    }
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        p->charge += sumTemp * mol1.ratio[i];
                    }
                }

                else if (mol2.charges > mol2.molrange[1]) {
                    // same as previous if statement, but if sum of new attempted charges in molecule2 falls above upper
                    // limit in molrange2
                    sumTemp = 0;
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        sumTemp += p->charge - (2 * mol2.max[i] - (p->charge + mol2.changeQ[i]));
                        p->charge = 2 * mol2.max[i] - (p->charge + mol2.changeQ[i]);
                    }
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        p->charge += sumTemp * mol1.ratio[i];
                    }
                }

                else {
                    // in case no boundaries were crossed, i.e. all new charges lies within their respective ranges
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        p->charge += mol1.changeQ[i];
                    }
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        p->charge += mol2.changeQ[i];
                    }
                }
                mol1.cdata.all = true;               // change all atoms in molecule1
                mol2.cdata.all = true;               // change all atoms in molecule2
                change.groups.push_back(mol1.cdata); // add to list of moved groups
                change.groups.push_back(mol2.cdata); // add to list of moved groups

            } else
                deltaq = 0;
        }
    }
}

void ChargeTransfer::_accept(Change&) { msqd += deltaq * deltaq; }
void ChargeTransfer::_reject(Change&) { msqd += 0; }

ChargeTransfer::ChargeTransfer(Space& spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
    mol1.cdata.internal = true;
    mol2.cdata.internal = true;
    // cdata1.atoms.resize(numOfAtoms1);
    // cdata2.atoms.resize(numOfAtoms2);
}

ChargeTransfer::ChargeTransfer(Space& spc) : ChargeTransfer(spc, "chargetransfer", "") {}

void QuadrantJump::_to_json(json& j) const {
    j = {{"dir", dir},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molecules[molid].name}};
    roundJSON(j, 3);
}
void QuadrantJump::_from_json(const json& j) {
    assert(!molecules.empty());
    const std::string molname = j.at("molecule");
    molid = findMoleculeByName(molname).id();
    dir = j.value("dir", Point(1, 1, 1));
    index = j.value("index", decltype(index)());
    if (repeat < 0) {
        auto v = spc.findMolecules(molid);
        repeat = std::distance(v.begin(), v.end());
    }
}
void QuadrantJump::_move(Change& change) {
    assert(molid >= 0);
    assert(!spc.groups.empty());
    assert(spc.geometry.getVolume() > 0);

    _sqd = 0.0;

    // pick random group from the system matching molecule type
    // TODO: This can be slow -- implement look-up-table in Space
    auto mollist = spc.findMolecules(molid, Space::Selection::ACTIVE); // list of molecules w. 'molid'
    if (not ranges::cpp20::empty(mollist)) {
        auto it = slump.sample(mollist.begin(), mollist.end());
        if (not it->empty()) {
            assert(it->id == molid);
            Point oldcm = it->mass_center;
            if (index.size() == 2) {
                auto cm_O = Geometry::massCenter(spc.particles.begin() + index[0], spc.particles.begin() + index[1] + 1,
                                                 spc.geometry.getBoundaryFunc());
                it->translate(-2 * spc.geometry.vdist(oldcm, cm_O).cwiseProduct(dir.cast<double>()),
                              spc.geometry.getBoundaryFunc());
            } else {
                it->translate(-2 * oldcm.cwiseProduct(dir.cast<double>()), spc.geometry.getBoundaryFunc());
            }
            _sqd = spc.geometry.sqdist(oldcm, it->mass_center); // squared displacement
            Change::GroupChange d;
            d.group_index = Faunus::distance(spc.groups.begin(), it); // integer *index* of moved group
            d.all = true;                                             // *all* atoms in group were moved
            change.groups.push_back(d);                               // add to list of moved groups

            assert(spc.geometry.sqdist(it->mass_center,
                                       Geometry::massCenter(it->begin(), it->end(), spc.geometry.getBoundaryFunc(),
                                                            -it->mass_center)) < 1e-9);
        }
    } else
        faunus_logger->warn("{0}: no molecules found", name);
}

QuadrantJump::QuadrantJump(Space& spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
}

QuadrantJump::QuadrantJump(Space& spc) : QuadrantJump(spc, "quadrantjump", "") {}

void AtomicSwapCharge::_to_json(json& j) const {
    j = {{"pH", pH},
         {"pka", pKa},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molname}};
    roundJSON(j, 3);
}
void AtomicSwapCharge::_from_json(const json& j) {
    assert(!molecules.empty());
    molname = j.at("molecule");
    molid = findMoleculeByName(molname).id();
    pH = j.at("pH").get<double>();
    pKa = j.at("pKa").get<double>();
    if (repeat < 0) {
        auto v = spc.findMolecules(molid);
        repeat = std::distance(v.begin(), v.end()); // repeat for each molecule...
        if (repeat > 0) {
            repeat = repeat * v.front().size(); // ...and for each atom
        }
    }
}
ParticleVector::iterator AtomicSwapCharge::randomAtom() {
    assert(molid >= 0);
    auto mollist = spc.findMolecules(molid); // all `molid` groups
    if (not ranges::cpp20::empty(mollist)) {
        auto git = slump.sample(mollist.begin(), mollist.end()); // random molecule iterator
        if (!git->empty()) {
            auto p = slump.sample(git->begin(), git->end());                 // random particle iterator
            cdata.group_index = Faunus::distance(spc.groups.begin(), git);   // integer *index* of moved group
            cdata.relative_atom_indices[0] = std::distance(git->begin(), p); // index of particle rel. to group
            return p;
        }
    }
    return spc.particles.end();
}
void AtomicSwapCharge::_move(Change& change) {
    _sqd = 0.0;
    auto p = randomAtom();
    if (p != spc.particles.end()) {
        // auto &g = spc.groups[cdata.index];
        double oldcharge = p->charge;
        p->charge = fabs(oldcharge - 1);
        _sqd = fabs(oldcharge - 1) - oldcharge;
        change.groups.push_back(cdata);   // add to list of moved groups
        _bias = _sqd * (pH - pKa) * ln10; // one may add bias here...
    }
}
double AtomicSwapCharge::bias(Change&, double, double) { return _bias; }
void AtomicSwapCharge::_accept(Change&) { msqd += _sqd; }
void AtomicSwapCharge::_reject(Change&) { msqd += 0; }

AtomicSwapCharge::AtomicSwapCharge(Space& spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
    cdata.relative_atom_indices.resize(1);
    cdata.internal = true;
}

AtomicSwapCharge::AtomicSwapCharge(Space& spc) : AtomicSwapCharge(spc, "swapcharge", "") {}

void TranslateRotate::_to_json(json& j) const {
    j = {{"dir", translational_direction},
         {"dp", translational_displacement},
         {"dprot", rotational_displacement},
         {"dirrot", fixed_rotation_axis},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(mean_squared_displacement.avg())},
         {"√⟨θ²⟩/°", std::sqrt(mean_squared_rotation_angle.avg()) / 1.0_deg},
         {"molecule", Faunus::molecules[molid].name}};
    roundJSON(j, 3);
}
void TranslateRotate::_from_json(const json& j) {
    const std::string molname = j.at("molecule");
    const auto molecule = findMoleculeByName(molname);
    if (molecule.atomic) {
        throw ConfigurationError("molecule '{}' cannot be atomic", molname);
    }
    molid = molecule.id();
    translational_direction = j.value("dir", Point(1, 1, 1));
    translational_displacement = j.at("dp").get<double>();
    rotational_displacement = j.at("dprot").get<double>();
    fixed_rotation_axis = j.value("dirrot", Point(0.0, 0.0, 0.0)); // predefined axis of rotation
    if (fixed_rotation_axis.count() > 0.0) {
        fixed_rotation_axis.normalize();
        faunus_logger->debug("{}: fixed rotation axis [{:.3f},{:.3f},{:.3f}]", name, fixed_rotation_axis.x(),
                             fixed_rotation_axis.y(), fixed_rotation_axis.z());
    }
    if (repeat < 0) {
        repeat = spc.numMolecules<Space::GroupType::ACTIVE>(molid);
        if (repeat == 0) {
            faunus_logger->warn("no initial '{}' molecules found; setting repeat to 1", molname);
            repeat = 1;
        } else {
            faunus_logger->debug("repeat = {} for molecule '{}'", repeat, molname);
        }
    }
}

/**
 * @todo `mollist` scales linearly w. system size -- implement look-up-table in Space?
 */
std::optional<std::reference_wrapper<Space::GroupType>> TranslateRotate::findRandomMolecule() {
    if (auto mollist = spc.findMolecules(molid, Space::Selection::ACTIVE); not ranges::cpp20::empty(mollist)) {
        if (auto group_it = slump.sample(mollist.begin(), mollist.end()); not group_it->empty()) {
            return *group_it;
        }
    }
    return std::nullopt;
}

/**
 * @param group Group to translate
 * @return Squared translation distance of mass center
 */
double TranslateRotate::translateMolecule(Space::GroupType& group) {
    if (translational_displacement > 0.0) { // translate
        const auto old_mass_center = group.mass_center;
        const auto displacement_vector =
            Faunus::randomUnitVector(slump, translational_direction) * translational_displacement * slump();

        group.translate(displacement_vector, spc.geometry.getBoundaryFunc());
        return spc.geometry.sqdist(old_mass_center, group.mass_center);
    } else {
        return 0.0;
    }
}

/**
 * @param group Group to rotate
 * @return Squared rotation angle around mass-center
 */
double TranslateRotate::rotateMolecule(Space::GroupType& group) {
    if (rotational_displacement <= pc::epsilon_dbl) {
        return 0.0;
    }
    Point rotation_axis;
    if (fixed_rotation_axis.count() > 0.0) { // fixed user-defined axis
        rotation_axis = fixed_rotation_axis;
    } else {
        rotation_axis = Faunus::randomUnitVector(slump);
    }
    const auto angle = rotational_displacement * (slump() - 0.5);
    const Eigen::Quaterniond quaternion(Eigen::AngleAxisd(angle, rotation_axis));
    group.rotate(quaternion, spc.geometry.getBoundaryFunc());
    return angle * angle;
}

void TranslateRotate::_move(Change& change) {
    if (auto group = findRandomMolecule()) { // note that group is of type std::optional
        latest_displacement_squared = translateMolecule(group->get());
        latest_rotation_angle_squared = rotateMolecule(group->get());
        if (latest_displacement_squared > 0.0 || latest_rotation_angle_squared > 0.0) { // report changes
            auto& change_data = change.groups.emplace_back();
            change_data.group_index = spc.getGroupIndex(group->get()); // integer *index* of moved group
            change_data.all = true;                                    // *all* atoms in group were moved
            change_data.internal = false;                              // internal energy is unchanged
        }
        checkMassCenter(group->get());
    } else {
        latest_displacement_squared = 0.0;   // these are used to track mean squared
        latest_rotation_angle_squared = 0.0; // translational and rotational displacements
    }
}

void TranslateRotate::checkMassCenter(const Space::GroupType& group) const {
    const auto allowed_threshold = 1e-6;
    const auto cm_recalculated =
        Geometry::massCenter(group.begin(), group.end(), spc.geometry.getBoundaryFunc(), -group.mass_center);
    const auto should_be_small = spc.geometry.sqdist(group.mass_center, cm_recalculated);
    if (should_be_small > allowed_threshold) {
        faunus_logger->error("{}: error calculating mass center for {}", name, group.traits().name);
        PQRWriter().save("mass-center-failure.pqr", spc.groups, spc.geometry.getLength());
        throw std::runtime_error("molecule likely too large for periodic boundaries; increase box size?");
    }
}

void TranslateRotate::_accept(Change&) {
    mean_squared_displacement += latest_displacement_squared;
    mean_squared_rotation_angle += latest_rotation_angle_squared;
}

void TranslateRotate::_reject(Change&) {
    mean_squared_displacement += 0.0;
    mean_squared_rotation_angle += 0.0;
}

TranslateRotate::TranslateRotate(Space& spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
}

TranslateRotate::TranslateRotate(Space& spc) : TranslateRotate(spc, "moltransrot", "") {}

/**
 * @param outside_rejection_probability Probability that any element outside the region will be rejected
 * @param region Region to preferentially pick from
 */
SmartMonteCarlo::SmartMonteCarlo(double outside_rejection_probability, std::unique_ptr<Region::RegionBase> region)
    : outside_rejection_probability(outside_rejection_probability), region(std::move(region)) {}

/**
 * @param number_total Total number of elements
 * @param number_inside Total number of elements inside region
 * @param direction Did element exit or enter the region? Returns zero if no boundary crossing
 */
double SmartMonteCarlo::bias(const int number_total, const int number_inside,
                             SmartMonteCarlo::BiasDirection direction) const {
    const auto p = outside_rejection_probability;
    switch (direction) {
    case BiasDirection::EXIT_REGION:
        return -std::log(p / (1.0 - (1.0 - p) / (p * number_total + (1.0 - p) * number_inside)));
    case BiasDirection::ENTER_REGION:
        return -std::log(1.0 / (1.0 + (1.0 - p) / (p * number_total + (1.0 - p) * number_inside)));
    case BiasDirection::NO_CROSSING:
        return 0.0;
    }
}

/**
 * This is called *after* the move and the `bias()` function will determine if the move
 * resulted in a flux over the region boundary and return the appropriate bias.
 */
double SmartTranslateRotate2::bias([[maybe_unused]] Change& change, [[maybe_unused]] double old_energy,
                                   [[maybe_unused]] double new_energy) {
    if (selection) {
        updateAverageCountInside(selection->number_inside);
        return smart_monte_carlo->bias(*selection);
    }
    return 0.0;
}

/**
 * This is used to sample the average number of groups inside the region. Once converged, the
 * expensive region search can be disabled by giving `fixed_count_inside` a mean value which will
 * be used for all further bias evaluations.
 */
void SmartTranslateRotate2::updateAverageCountInside(int count_inside) {
    mean_count_inside += static_cast<double>(count_inside);
    if (mean_count_inside.size() % bias_update_interval == 0) {
        if (mean_count_inside.stdev() / mean_count_inside.avg() < 0.05) {
            smart_monte_carlo->fixed_count_inside = static_cast<int>(mean_count_inside.avg());
            faunus_logger->info("Stopping bias update since threshold reached.");
        }
    }
}

/**
 * Upon calling `select()`, the `outside_rejection_probability` is used to exclude particles
 * outside and may thus often return `std::nullopt`.
 */
std::optional<std::reference_wrapper<Space::GroupType>> SmartTranslateRotate2::findRandomMolecule() {
    if (auto mollist = spc.findMolecules(molid, Space::Selection::ACTIVE); ranges::cpp20::empty(mollist)) {
        selection = std::nullopt;
    } else {
        selection = smart_monte_carlo->select<Group>(mollist, slump);
        if (selection) {
            return *(selection->item);
        }
    }
    return std::nullopt;
}
void SmartTranslateRotate2::_to_json(json& j) const {
    TranslateRotate::_to_json(j);
    j["smart monte carlo"] = static_cast<json>(*smart_monte_carlo);
}
} // namespace Faunus::Move

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] TranslateRotate") {
    using namespace Faunus;
    CHECK(!atoms.empty());     // set in a previous test
    CHECK(!molecules.empty()); // set in a previous test

    Space spc;
    Move::TranslateRotate mv(spc);
    json j = R"( {"molecule":"A", "dp":1.0, "dprot":0.5, "dir":[0,1,0], "repeat":2 })"_json;
    mv.from_json(j);

    j = json(mv).at(mv.getName());
    CHECK(j.at("molecule") == "A");
    // CHECK(j.at("dir") == Point(0, 1, 0));
    CHECK(j.at("dp") == 1.0);
    CHECK(j.at("repeat") == 2);
    CHECK(j.at("dprot") == 0.5);
}
#endif

namespace Faunus::Move {

void SmartTranslateRotate::_to_json(json& j) const {
    j = {//{"Number of counts inside geometry", cntInner},
         //{"Number of counts outside geometry", cnt - cntInner},
         {"dir", dir},
         {"dp", dptrans},
         {"dprot", dprot},
         {"p", p},
         {"rx", r_x},
         {"ry", r_y},
         {"molid", molid},
         {"refid", refid1},
         {"refid", refid2},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molecules[molid].name},
         {"ref1", atoms[refid1].name},
         {"ref2", atoms[refid2].name}};
    roundJSON(j, 3);
}
void SmartTranslateRotate::_from_json(const json& j) {
    assert(!molecules.empty());
    const std::string molname = j.at("molecule");
    const std::string refname1 = j.at("ref1");
    const std::string refname2 = j.at("ref2");
    molid = findMoleculeByName(molname).id();
    refid1 = findAtomByName(refname1).id();
    refid2 = findAtomByName(refname2).id();
    dir = j.value("dir", Point(1, 1, 1));
    dprot = j.at("dprot");
    dptrans = j.at("dp");
    p = j.at("p");
    r_x = j.at("rx");  // length of ellipsoidal radius along axis connecting reference atoms (in Å)
    r_y = j.at("ry");  // length of ellipsoidal radius perpendicular to axis connecting reference atoms (in Å)
    rsd = j.at("rsd"); // threshold for relative standard deviation of molecules inside geometry. When it goes below
    // this value, a constant bias is used

    if (repeat < 0) {
        auto v = spc.findMolecules(molid);
        repeat = std::distance(v.begin(), v.end());
    }
}

void SmartTranslateRotate::_move(Change& change) {
    squared_displacement = 0.0;
    bias_energy = 0.0;

    auto mollist = spc.findMolecules(molid, Space::Selection::ACTIVE); // list of molecules w. 'molid'
    if (ranges::cpp20::empty(mollist)) {
        return;
    }
    auto selected_molecule =
        slump.sample(mollist.begin(), mollist.end()); // chosing random molecule in group of type molname
    if (selected_molecule != mollist.end() || selected_molecule->empty()) {
        return;
    }

    const auto is_inside = createInsideLambdaFunc();
    const auto molecule_is_inside = is_inside(selected_molecule->mass_center);
    cnt++; // total number of counts

    if (!molecule_is_inside && slump() < p) {
        return; // reject outside molecule with probability p
    }

    if (update_bias) { // continuing to adjust bias according to number of molecules inside and outside geometry
        auto mass_centers = mollist | ranges::cpp20::views::transform(&Group::mass_center);
        num_molecules_total = ranges::cpp20::distance(mollist.begin(), mollist.end());
        num_molecules_inside = ranges::cpp20::count_if(mass_centers, is_inside);
        num_molecules_outside = num_molecules_total - num_molecules_inside;
        updateAveragesAndBias();
    }

    if (dptrans > 0) { // translate
        Point oldcm = selected_molecule->mass_center;
        Point dp = randomUnitVector(slump, dir) * dptrans * slump();

        selected_molecule->translate(dp, spc.geometry.getBoundaryFunc());
        squared_displacement = spc.geometry.sqdist(oldcm, selected_molecule->mass_center); // squared displacement
    }

    if (dprot > 0) { // rotate
        Point u = randomUnitVector(slump);
        double angle = dprot * (slump() - 0.5);
        Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, u));
        selected_molecule->rotate(Q, spc.geometry.getBoundaryFunc());
    }

    if (dptrans > 0 || dprot > 0) { // define changes
        Change::GroupChange d;
        d.group_index = Faunus::distance(spc.groups.begin(), selected_molecule); // integer *index* of moved group
        d.all = true;                                             // *all* atoms in group were moved
        change.groups.push_back(d);                               // add to list of moved groups
    } else {
        return;
    }

    assert(spc.geometry.sqdist(selected_molecule->mass_center,
                               Geometry::massCenter(selected_molecule->begin(), selected_molecule->end(),
                                                    spc.geometry.getBoundaryFunc(), -selected_molecule->mass_center)) <
           1e-6);

    const auto moved_molecule_is_inside = is_inside(selected_molecule->mass_center);
    const double N = update_bias ? static_cast<double>(num_molecules_inside) : average_num_molecules_inside;
    if (molecule_is_inside && !moved_molecule_is_inside) {
        // inside --> outside
        bias_energy = -std::log(p / (1.0 - (1.0 - p) / (p * num_molecules_total + (1.0 - p) * N)));
    } else if (!molecule_is_inside && moved_molecule_is_inside) {
        // outside --> inside
        bias_energy = -std::log(1.0 / (1.0 + (1.0 - p) / (p * num_molecules_total + (1.0 - p) * N)));
    } else {
        // boundary not crossed
        bias_energy = 0.0;
    }
}
void SmartTranslateRotate::updateAveragesAndBias() {
    mean_num_molecules_inside += num_molecules_inside;
    mean_num_molecules_outside += num_molecules_outside;

    if (cnt % 100 == 0) {
        countNin_avgBlocks += mean_num_molecules_inside.avg();
        countNout_avgBlocks += mean_num_molecules_outside.avg();
    }

    if (cnt % 100000 == 0) {
        average_num_molecules_inside = countNin_avgBlocks.avg(); // block average number of molecules inside geometry
        if (countNin_avgBlocks.stdev() / average_num_molecules_inside <
            rsd) { // if block standard deviation is below specified threshold
            faunus_logger->info(
                "Stopping bias update since rsd ({}) lower than threshold. Avg. molecules inside region = {}",
                countNin_avgBlocks.stdev() / average_num_molecules_inside, average_num_molecules_inside);
            update_bias = false; // stop updating bias, use constant value
        }
    }
}

/**
 * This creates a lambda function to determine if a given position is inside
 * the region or not. Note that this randomly generates a region that is then
 * fixed in the returned lambda.
 */
std::function<bool(const Point&)> SmartTranslateRotate::createInsideLambdaFunc() {
    auto reflist1 = spc.findAtoms(refid1); // list of atoms w. 'refid1'
    auto reflist2 = spc.findAtoms(refid2); // list of atoms w. 'refid2'
    auto ref1 = slump.sample(reflist1.begin(), reflist1.end());
    auto ref2 = slump.sample(reflist2.begin(), reflist2.end());
    Point cylAxis = spc.geometry.vdist(ref2->pos, ref1->pos) * 0.5; // half vector between reference atoms
    Point origo = ref2->pos - cylAxis; // coordinates of middle point between reference atoms: new origo
    if (r_x < cylAxis.norm()) {        // checking so that a is larger than length of cylAxis
        throw std::runtime_error(
            "specified radius of ellipsoid along the axis connecting reference atoms (rx) must be larger or equal "
            "to half the distance between reference atoms. Specified radius is " +
            std::to_string(r_x) + " Å whereas half the distance between reference atoms is " +
            std::to_string(cylAxis.norm()) + "Å");
    }

    return [=](const Point& position) {
        auto molV = spc.geometry.vdist(position, origo); // vector between selected molecule and center of geometry
        auto cosTheta = molV.dot(cylAxis) / molV.norm() / cylAxis.norm(); // cosinus of angle between coordinate
                                                                          // vector of selected molecule and axis
                                                                          // connecting reference atoms
        auto theta = acos(cosTheta); // angle between coordinate vector of sel. molecule and axis connecting ref. atoms
        auto x = cosTheta * molV.norm();   // x coordinate of selected molecule with respect to center of geometry
                                           // (in plane including vectors molV and cylAxis)
        auto y = sin(theta) * molV.norm(); // y coordinate of selected molecule with respect to center of geometry
                                           // (in plane including vectors molV and cylAxis)
        auto coord =
            x * x / (r_x * r_x) + y * y / (r_y * r_y); // calculating normalized coordinate with respect to
                                                       // dimensions of geometry (>1.0 → outside, <1.0 → inside)
        return coord <= 1.0;
    };
}

double SmartTranslateRotate::bias([[maybe_unused]] Change& change, double, double) { return bias_energy; }

SmartTranslateRotate::SmartTranslateRotate(Space& spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
}

SmartTranslateRotate::SmartTranslateRotate(Space& spc) : SmartTranslateRotate(spc, "smartmoltransrot", "") {}
void SmartTranslateRotate::_accept([[maybe_unused]] Change& change) { msqd += squared_displacement; }
void SmartTranslateRotate::_reject([[maybe_unused]] Change& change) { msqd += 0.0; }

void ConformationSwap::_to_json(json& j) const {
    j = {{"molid", molid},
         {"molecule", Faunus::molecules.at(molid).name},
         {"keeppos", inserter.keep_positions},
         {"copy_policy", copy_policy}};
    roundJSON(j, 3);
}

void ConformationSwap::_from_json(const json& j) {
    const auto molecule_name = j.at("molecule").get<std::string>();
    const auto molecule = Faunus::findMoleculeByName(molecule_name);
    molid = molecule.id();
    if (molecule.conformations.size() < 2) {
        throw ConfigurationError("minimum two conformations required for {}", molecule_name);
    }
    checkConformationSize(); // do conformations fit periodic boundaries?
    inserter.keep_positions = j.value("keeppos", false);
    copy_policy = j.value("copy_policy", CopyPolicy::ALL);
    if (copy_policy == CopyPolicy::INVALID) {
        throw ConfigurationError("invalid copy policy");
    }
    setRepeat();
}

void ConformationSwap::setRepeat() {
    assert(molid >= 0);
    if (repeat < 0) { // negative value signals repeat = N number of molecules
        auto groups = spc.findMolecules(molid, Space::Selection::ALL);
        repeat = std::distance(groups.begin(), groups.end());
        if (repeat == 0) {
            faunus_logger->warn("{}: no molecules found; repeat set to ZERO", name, repeat);
        }
    }
}
void ConformationSwap::_move(Change& change) {
    auto groups = spc.findMolecules(molid, Space::Selection::ACTIVE);
    if (auto group = slump.sample(groups.begin(), groups.end()); group != groups.end()) {
        inserter.offset = group->mass_center; // insert on top of mass center
        auto particles = inserter(spc.geometry, Faunus::molecules[molid], spc.particles); // new conformation
        if (particles.size() == group->size()) {
            checkMassCenterDrift(group->mass_center, particles); // throws if not OK
            copyConformation(particles, group->begin());
            group->conformation_id = Faunus::molecules[molid].conformations.getLastIndex(); // store conformation id
            registerChanges(change, *group);                                                // update change object
        } else {
            throw std::out_of_range(name + ": conformation atom count mismatch");
        }
    }
}

/**
 * This will copy the new conformation onto the destination group. By default
 * all information is copied, but can be limited to positions, only
 *
 * @param particles Source particles
 * @param destination Iterator to first particle in destination
 */
void ConformationSwap::copyConformation(ParticleVector& particles, ParticleVector::iterator destination) const {
    std::function<void(const Particle&, Particle&)> copy_function; // how to copy particle information
    switch (copy_policy) {
    case CopyPolicy::ALL:
        copy_function = [](const Particle& src, Particle& dst) { dst = src; };
        break;
    case CopyPolicy::POSITIONS:
        copy_function = [](const Particle& src, Particle& dst) { dst.pos = src.pos; };
        break;
    case CopyPolicy::CHARGES:
        copy_function = [](const Particle& src, Particle& dst) { dst.charge = src.charge; };
        break;
    default:
        throw std::runtime_error("invalid copy policy");
    }

    std::for_each(particles.cbegin(), particles.cend(), [&](const Particle& source) {
        copy_function(source, *destination);
        destination++;
    });
}

void ConformationSwap::registerChanges(Change& change, const Space::GroupType& group) const {
    auto& group_change = change.groups.emplace_back();
    group_change.group_index = spc.getGroupIndex(group); // index of moved group
    group_change.all = true;                             // all atoms in group were moved
    group_change.internal = false;                       // skip internal energy calculation
}
/**
 * @throw if there's a mass-center drift
 *
 * Move shouldn't move mass centers, so let's check if this is true
 */
void ConformationSwap::checkMassCenterDrift(const Point& old_mass_center, const ParticleVector& particles) {
    switch (copy_policy) {
    case CopyPolicy::CHARGES: // positions untouched; no check needed
        return;
    default:
        const auto max_allowed_distance = 1.0e-6;
        const auto new_mass_center =
            Geometry::massCenter(particles.begin(), particles.end(), spc.geometry.getBoundaryFunc(), -old_mass_center);
        if ((new_mass_center - old_mass_center).norm() > max_allowed_distance) {
            throw std::runtime_error(name + ": unexpected mass center movement");
        }
    }
}

ConformationSwap::ConformationSwap(Space& spc, const std::string& name, const std::string& cite)
    : MoveBase(spc, name, cite) {}

ConformationSwap::ConformationSwap(Space& spc) : ConformationSwap(spc, "conformationswap", "doi:10/dmc3") {
    repeat = -1; // meaning repeat n times
    inserter.dir = Point::Zero();
    inserter.rotate = true;
    inserter.allow_overlap = true;
}

/**
 * Checks if any two particles in any conformation is father away than half the shortest cell length
 */
void ConformationSwap::checkConformationSize() const {
    assert(molid >= 0);

    const auto is_periodic = spc.geometry.boundaryConditions().isPeriodic();
    if (is_periodic.cast<int>().sum() == 0) { // if cell has no periodicity ...
        return;                               // ... then no need to check further
    }

    // find smallest periodic side-length
    const auto infinity = Point::Constant(pc::infty);
    const auto max_allowed_separation =
        (is_periodic.array() == true).select(spc.geometry.getLength(), infinity).minCoeff() * 0.5;

    auto find_max_distance = [&geometry = spc.geometry](const auto& positions) {
        double max_squared_distance = 0.0;
        for (auto i = positions.begin(); i != positions.end(); ++i) {
            for (auto j = i; ++j != positions.end();) {
                max_squared_distance = std::max(max_squared_distance, geometry.sqdist(*i, *j));
            }
        }
        return std::sqrt(max_squared_distance);
    }; // find internal maximum distance in a set of positions

    size_t conformation_id = 0;
    for (const auto& conformation : Faunus::molecules.at(molid).conformations.data) {
        const auto positions = conformation | ranges::cpp20::views::transform(&Particle::pos);
        const auto max_separation = find_max_distance(positions);
        if (max_separation > max_allowed_separation) {
            faunus_logger->warn("particles in conformation {} separated by {:.3f} Å which *may* break periodic "
                                "boundaries. If so, you'll know.",
                                conformation_id, max_separation);
        }
        conformation_id++;
    }
}

} // namespace Faunus::Move
