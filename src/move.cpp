#include <doctest/doctest.h>
#include "core.h"
#include "move.h"
#include "speciation.h"
#include "clustermove.h"
#include "chainmove.h"
#include "forcemove.h"
#include "montecarlo.h"
#include "aux/iteratorsupport.h"
#include "aux/eigensupport.h"
#include "spdlog/spdlog.h"

namespace Faunus::Move {

Random MoveBase::slump; // static instance of Random (shared for all moves)

void MoveBase::from_json(const json &j) {
    if (auto it = j.find("repeat"); it != j.end()) {
        if (it->is_number()) {
            repeat = it->get<int>();
        } else if (it->is_string()) {
            if (it->get<std::string>() == "N") {
                repeat = -1;
            } else {
                throw std::runtime_error("Unknown string value for 'repeat'");
            }
        }
    }
    _from_json(j);
    if (repeat < 0) {
        repeat = 0;
    }
}

void MoveBase::to_json(json &j) const {
    _to_json(j);
    if (timer_move.result() > 0.01) // only print if more than 1% of the time
        j["relative time (without energy calc)"] = timer_move.result();
    if (timer.result() > 0.01) // only print if more than 1% of the time
        j["relative time"] = timer.result();
    j["acceptance"] = double(accepted) / cnt;
    j["repeat"] = repeat;
    j["moves"] = cnt;
    if (!cite.empty())
        j["cite"] = cite;
    _roundjson(j, 3);
}

void MoveBase::move(Change &change) {
    timer.start();
    timer_move.start();
    cnt++;
    change.clear();
    _move(change);
    if (change.empty())
        timer.stop();
    timer_move.stop();
}

void MoveBase::accept(Change &c) {
    accepted++;
    _accept(c);
    timer.stop();
}

void MoveBase::reject(Change &c) {
    rejected++;
    _reject(c);
    timer.stop();
}

double MoveBase::bias(Change &, double, double) {
    return 0; // du
}

void MoveBase::_accept(Change &) {}

void MoveBase::_reject(Change &) {}

void from_json(const json &j, MoveBase &m) { m.from_json(j); }

void to_json(json &j, const MoveBase &m) {
    assert(!m.name.empty());
    m.to_json(j[m.name]);
}

ReplayMove::ReplayMove(Space &spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {}

ReplayMove::ReplayMove(Space &spc) : ReplayMove(spc, "replay", "") {}

void ReplayMove::_to_json(json &j) const { j["file"] = reader->filename; }

void ReplayMove::_from_json(const json &j) { reader = std::make_shared<XTCReader>(j.at("file")); }

void ReplayMove::_move(Change &change) {
    assert(reader != nullptr);
    if (!end_of_trajectory) {
        if (reader->read(frame.step, frame.timestamp, frame.box, spc.positions().begin(), spc.positions().end())) {
            spc.geo.setLength(frame.box);
            change.all = true;
        } else {
            // nothing to do, simulation shall stop
            end_of_trajectory = true;
            mcloop_logger->warn("No more frames to read from {}. Running on empty.", reader->filename);
        }
    }
}

double ReplayMove::bias(Change &, double, double) {
    return force_accept; // always accept
}

void AtomicTranslateRotate::_to_json(json &j) const {
    j = {{"dir", directions},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(mean_square_displacement.avg())},
         {"molecule", molecule_name}};
    _roundjson(j, 3);
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
        auto mollist = spc.findMolecules(molid, Space::ALL);
        repeat = std::distance(mollist.begin(), mollist.end()); // repeat for each molecule...
        if (repeat > 0) {
            repeat = repeat * mollist.front().size(); // ...and for each atom
        }
    }
}

void AtomicTranslateRotate::translateParticle(ParticleVector::iterator particle, double displacement) {
    auto old_position = particle->pos; // backup old position
    particle->pos += ranunit(slump, directions) * displacement * slump();
    spc.geo.boundary(particle->pos);
    _sqd = spc.geo.sqdist(old_position, particle->pos); // square displacement

    if (auto &group = spc.groups[cdata.index]; group.atomic == false) { // update COM for molecular groups
        group.cm = Geometry::massCenter(group.begin(), group.end(), spc.geo.getBoundaryFunc(), -group.cm);
#ifndef NDEBUG
        // Perform test to see if the move violates PBC
        auto old_mass_center = group.cm;                              // backup mass center
        group.translate(-old_mass_center, spc.geo.getBoundaryFunc()); // translate to {0,0,0}
        double should_be_zero = spc.geo.sqdist({0, 0, 0}, Geometry::massCenter(group.begin(), group.end()));
        if (should_be_zero > 1e-6) {
            throw std::runtime_error("atomic move too large");
        } else {
            group.translate(old_mass_center, spc.geo.getBoundaryFunc());
        }
#endif
    }
}

void AtomicTranslateRotate::_move(Change &change) {
    if (auto particle = randomAtom(); particle != spc.p.end()) {
        double translational_displacement = atoms.at(particle->id).dp;
        double rotational_displacement = atoms.at(particle->id).dprot;

        assert(translational_displacement >= 0.0);

        if (translational_displacement > 0.0) { // translate
            translateParticle(particle, translational_displacement);
        }

        if (rotational_displacement > 0.0) { // rotate
            Point u = ranunit(slump);
            double angle = rotational_displacement * (slump() - 0.5);
            Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, u));
            particle->rotate(Q, Q.toRotationMatrix());
        }

        if (translational_displacement > 0.0 or rotational_displacement > 0.0) {
            change.groups.push_back(cdata); // add to list of moved groups
        }
    } else {
        _sqd = 0.0; // no particle found --> no movement
    }
}

void AtomicTranslateRotate::_accept(Change &) { mean_square_displacement += _sqd; }
void AtomicTranslateRotate::_reject(Change &) { mean_square_displacement += 0; }

AtomicTranslateRotate::AtomicTranslateRotate(Space &spc, std::string name, std::string cite)
    : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
    cdata.atoms.resize(1);
    cdata.internal = true;
}

AtomicTranslateRotate::AtomicTranslateRotate(Space &spc) : AtomicTranslateRotate(spc, "transrot", "") {}

/**
 * For atomic groups, select `ALL` since these may be partially filled and thereby
 * appear inactive. Note also that only one instance of atomic molecules can exist.
 * For molecular groups, select only active ones.
 *
 * @return Iterator to particle to move; `end()` if nothing selected
 */
ParticleVector::iterator AtomicTranslateRotate::randomAtom() {
    assert(molid >= 0);
    auto particle = spc.p.end(); // particle iterator
    auto selection = (Faunus::molecules[molid].atomic) ? Space::ALL : Space::ACTIVE;
    auto mollist = spc.findMolecules(molid, selection);
    if (auto group = slump.sample(mollist.begin(), mollist.end()); group != mollist.end()) { // random molecule
        if (not group->empty()) {
            particle = slump.sample(group->begin(), group->end());     // random particle
            cdata.index = Faunus::distance(spc.groups.begin(), group); // index of touched group
            cdata.atoms[0] = std::distance(group->begin(), particle);  // index of moved particle relative to group
        }
    }
    return particle;
}

Propagator::Propagator(const json& j, Space& spc, Energy::Hamiltonian& pot, [[maybe_unused]] MPI::MPIController& mpi) {
    if (j.count("random") == 1) {
        MoveBase::slump = j["random"]; // slump is static --> shared for all moves
        Faunus::random = j["random"];
    }

    for (auto& j_move : j.at("moves")) { // loop over move list
        try {
            const auto& [key, j_params] = jsonSingleItem(j_move);
            if (key == "moltransrot") {
                _moves.emplace_back<Move::TranslateRotate>(spc);
            } else if (key == "smartmoltransrot") {
                _moves.emplace_back<Move::SmartTranslateRotate>(spc);
            } else if (key == "conformationswap") {
                _moves.emplace_back<Move::ConformationSwap>(spc);
            } else if (key == "transrot") {
                _moves.emplace_back<Move::AtomicTranslateRotate>(spc);
            } else if (key == "pivot") {
                _moves.emplace_back<Move::PivotMove>(spc);
            } else if (key == "crankshaft") {
                _moves.emplace_back<Move::CrankshaftMove>(spc);
            } else if (key == "volume") {
                _moves.emplace_back<Move::VolumeMove>(spc);
            } else if (key == "charge") {
                _moves.emplace_back<Move::ChargeMove>(spc);
            } else if (key == "chargetransfer") {
                _moves.emplace_back<Move::ChargeTransfer>(spc);
            } else if (key == "rcmc") {
                _moves.emplace_back<Move::SpeciationMove>(spc);
            } else if (key == "quadrantjump") {
                _moves.emplace_back<Move::QuadrantJump>(spc);
            } else if (key == "cluster") {
                _moves.emplace_back<Move::Cluster>(spc);
            } else if (key == "replay") {
                _moves.emplace_back<Move::ReplayMove>(spc);
            } else if (key == "langevin_dynamics") {
                auto move_ptr = std::make_shared<Move::LangevinDynamics>(spc, pot);
                _moves.push_back<Move::LangevinDynamics>(move_ptr);
            } else if (key == "temper") {
#ifdef ENABLE_MPI
                _moves.emplace_back<Move::ParallelTempering>(spc, mpi);
                _moves.back()->repeat = 0; // this gives the move ZERO weight
#else
                throw ConfigurationError("'{}': not included - recompile Faunus with MPI support", key);
#endif
            } else {
                throw ConfigurationError("'{}': unknown move", key);
            }
            // append additional moves to the if-chain

            try {
                _moves.back()->from_json(j_params);
                addWeight(_moves.back()->repeat);
            } catch (std::exception& e) {
                usageTip.pick(key);
                throw ConfigurationError("'{}': {}", key, e.what());
            }
        } catch (std::exception& e) {
            throw ConfigurationError("move: {}", e.what()).attachJson(j_move);
        }
    }
}

void Propagator::addWeight(double weight) {
    _weights.push_back(weight);
    distribution = std::discrete_distribution<>(_weights.begin(), _weights.end());
    _repeat = int(std::accumulate(_weights.begin(), _weights.end(), 0.0));
}

void to_json(json &j, const Propagator &propagator) { j = propagator._moves; }

#ifdef ENABLE_MPI

void ParallelTempering::_to_json(json &j) const {
    j = {{"replicas", mpi.nproc()}, {"datasize", particle_transmitter.getFormat()}};
    json &_j = j["exchange"];
    _j = json::object();
    for (const auto &[id, acceptance] : acceptance_map) {
        _j[id] = {{"attempts", acceptance.cnt}, {"acceptance", acceptance.avg()}};
    }
}

void ParallelTempering::findPartner() {
    auto true_or_false = static_cast<bool>(mpi.random.range(0, 1));
    int rank_increment = true_or_false ? 1 : -1;
    if (mpi.rank() % 2 == 0) { // even replica
        partner = mpi.rank() + rank_increment;
    } else { // odd replica
        partner = mpi.rank() - rank_increment;
    }
}

bool ParallelTempering::goodPartner() {
    if (partner >= 0 && partner < mpi.nproc() && partner != mpi.rank()) {
        return true;
    } else {
        partner = -1;
        return false;
    }
}

/**
 * This will exchange the states between two partner replicas and set the change object accordingy
 */
void ParallelTempering::exchangeState(Change &change) {
    assert(partner != -1);
    auto old_volume = spc.geo.getVolume();
    particle_transmitter.sendExtra.at(VOLUME) = old_volume;      // copy current volume for sending
    partner_particles->resize(spc.p.size());                     // temparary storage
    particle_transmitter.recv(mpi, partner, *partner_particles); // receive particles
    particle_transmitter.send(mpi, spc.p, partner);              // send everything
    particle_transmitter.waitrecv();
    particle_transmitter.waitsend();

    auto new_volume = particle_transmitter.recvExtra.at(VOLUME);
    if (new_volume < very_small_volume || spc.p.size() != partner_particles->size()) {
        MPI_Abort(mpi.comm, 1);
    } else {
        change.all = true;
        if (std::fabs(new_volume - old_volume) > pc::epsilon_dbl) {
            change.dV = true;
            spc.geo.setVolume(new_volume);
        }
        spc.updateParticles(partner_particles->begin(), partner_particles->end(), spc.p.begin());
    }
}

void ParallelTempering::_move(Change &change) {
    mpi.barrier(); // wait until all ranks reach here
    findPartner();
    if (goodPartner()) {
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
double ParallelTempering::exchangeEnergy(double energy_change) {
    assert(partner >= 0);
    std::vector<MPI::FloatTransmitter::floatp> energy_change_vector = {energy_change};
    auto energy_change_partner = float_transmitter.swapf(mpi, energy_change_vector, partner);
    return energy_change_partner.at(0); // return partner energy change
}

/**
 * The bias() function takes the current old and new energy and exchanges
 * the resulting energy change with the partner replica. The change in the
 * replica is returned.
 *
 * @todo Here we could run a custom Metropolis criterion and return +/- infinity
 * to trigger accept/reject. This would remedy the dangerous expectation that
 * the states of the random number generators (Movebase) are aligned on all nodes.
 * This would however ignore other bias contributions, particularly it would prove
 * problematic with grand canonical moves.
 */
double ParallelTempering::bias(Change &, double uold, double unew) {
    assert(partner != -1);
    if constexpr (false) {
        // todo: add sanity check for random number generator state in partnering replicas.
        return exchangeEnergy(unew - uold); // exchange change with partner (MPI)
    } else {
        double energy_change = unew - uold;
        double partner_energy_change = exchangeEnergy(energy_change);
        if (MetropolisMonteCarlo::metropolis(energy_change + partner_energy_change)) {
            return pc::neg_infty; // accept!
        } else {
            return pc::infty; // reject!
        }
    }
}

std::string ParallelTempering::id() const {
    assert(partner != -1);
    // note `std::minmax(a,b)` takes _references_; the initializer list (used here) takes a _copy_
    const auto pair = std::minmax({mpi.rank(), partner});
    return fmt::format("{} <-> {}", pair.first, pair.second);
}

void ParallelTempering::_accept(Change &) {
    acceptance_map[id()] += 1;
}
void ParallelTempering::_reject(Change &) {
    acceptance_map[id()] += 0;
}

void ParallelTempering::_from_json(const json &j) {
    particle_transmitter.setFormat(j.value("format", std::string("XYZQI")));
}

ParallelTempering::ParallelTempering(Space &spc, MPI::MPIController &mpi)
    : MoveBase(spc, "temper", "doi:10/b3vcw7"), mpi(mpi) {
    if (mpi.nproc() < 2) {
        throw std::runtime_error(name + " requires two or more MPI processes");
    }
    partner_particles = std::make_shared<ParticleVector>();
    partner_particles->reserve(spc.p.size());
    particle_transmitter.recvExtra.resize(1);
    particle_transmitter.sendExtra.resize(1);
}

/**
 * At the end of the simulation, the state of the random number generators
 * must be the same on all ranks. Run with verbose logging (trace) and observe output!
 */
ParallelTempering::~ParallelTempering() {
#ifndef NDEBUG
    faunus_logger->trace("mpi{}: last random number (Movebase) = {}", mpi.rank(), slump());
    faunus_logger->trace("mpi{}: last random number (Temper) = {}", mpi.rank(), random());
    faunus_logger->trace("mpi{}: last random number (MPI) = {}", mpi.rank(), mpi.random());
#endif
}

#endif

void VolumeMove::_to_json(json &j) const {
    using namespace u8;
    if (cnt > 0) {
        j = {{"dV", volume_displacement_factor},
             {"method", method->first},
             {bracket("V"), mean_volume.avg()},
             {rootof + bracket(Delta + "V" + squared), std::sqrt(mean_square_volume_change.avg())},
             {cuberoot + rootof + bracket(Delta + "V" + squared),
              std::cbrt(std::sqrt(mean_square_volume_change.avg()))}};
        _roundjson(j, 3);
    }
}
void VolumeMove::_from_json(const json &j) {
    method = methods.find(j.value("method", "isotropic"));
    if (method == methods.end()) {
        throw ConfigurationError("unknown volume change method");
    }
    volume_displacement_factor = j.at("dV").get<double>();
}
void VolumeMove::_move(Change &change) {
    if (volume_displacement_factor > 0.0) {
        change.dV = true;
        change.all = true;
        old_volume = spc.geo.getVolume();
        new_volume = std::exp(std::log(old_volume) + (slump() - 0.5) * volume_displacement_factor);
        volume_change = new_volume - old_volume;
        spc.scaleVolume(new_volume, method->second);
    } else
        volume_change = 0.0;
}
void VolumeMove::_accept(Change &) {
    mean_square_volume_change += volume_change * volume_change;
    mean_volume += new_volume;
    assert(std::fabs(spc.geo.getVolume() - new_volume) < 1.0e-9);
}

VolumeMove::VolumeMove(Space &spc, std::string name, std::string cite) : MoveBase(spc, name, cite) { repeat = 1; }

VolumeMove::VolumeMove(Space &spc) : VolumeMove(spc, "volume", "") {}

void VolumeMove::_reject(Change &) {
    mean_square_volume_change += 0.0;
    mean_volume += old_volume;
    assert(std::fabs(spc.geo.getVolume() - old_volume) < 1.0e-9);
}

void ChargeMove::_to_json(json &j) const {
    using namespace u8;
    j = {{"index", atomIndex},
         {"dq", dq},
         {rootof + bracket(Delta + "q" + squared), std::sqrt(msqd.avg())},
         {cuberoot + rootof + bracket(Delta + "q" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
    _roundjson(j, 3);
}
void ChargeMove::_from_json(const json &j) {
    dq = j.at("dq").get<double>();
    atomIndex = j.at("index").get<int>();
    auto git = spc.findGroupContaining(spc.p[atomIndex]);                    // group containing atomIndex
    cdata.index = std::distance(spc.groups.begin(), git);                    // integer *index* of moved group
    cdata.atoms[0] = std::distance(git->begin(), spc.p.begin() + atomIndex); // index of particle rel. to group
}
void ChargeMove::_move(Change &change) {
    if (dq > 0) {
        auto &p = spc.p[atomIndex]; // refence to particle
        double qold = p.charge;
        p.charge += dq * (slump() - 0.5);
        deltaq = p.charge - qold;
        change.groups.push_back(cdata); // add to list of moved groups
    } else
        deltaq = 0;
}
void ChargeMove::_accept(Change &) { msqd += deltaq * deltaq; }
void ChargeMove::_reject(Change &) { msqd += 0; }

ChargeMove::ChargeMove(Space &spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = 1;
    cdata.internal = true; // the group is internally changed
    cdata.atoms.resize(1); // we change exactly one atom
}

ChargeMove::ChargeMove(Space &spc) : ChargeMove(spc, "charge", "") {}

void ChargeTransfer::_to_json(json &j) const {
    using namespace u8;
    j = {{"dq", dq},
         {rootof + bracket(Delta + "q" + squared), std::sqrt(msqd.avg())},
         {cuberoot + rootof + bracket(Delta + "q" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
    _roundjson(j, 3);
}
void ChargeTransfer::_from_json(const json &j) {
    dq = j.at("dq").get<double>();
    mol1.molname = j.at("mol1"); // string containing name of molecule 1
    mol2.molname = j.at("mol2"); // string containing name of molecule 2
    mol1.id = findMoleculeByName(mol1.molname).id(); // group containing mol1.molname
    mol2.id = findMoleculeByName(mol2.molname).id(); // group containing mol2.molname
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

void ChargeTransfer::_move(Change &change) {
    auto mollist1 = spc.findMolecules(mol1.id, Space::ACTIVE);
    auto mollist2 = spc.findMolecules(mol2.id, Space::ACTIVE);
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
                mol1.cdata.index = Faunus::distance(spc.groups.begin(), git1);
                mol2.cdata.index = Faunus::distance(spc.groups.begin(), git2);

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

void ChargeTransfer::_accept(Change &) { msqd += deltaq * deltaq; }
void ChargeTransfer::_reject(Change &) { msqd += 0; }

ChargeTransfer::ChargeTransfer(Space &spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
    mol1.cdata.internal = true;
    mol2.cdata.internal = true;
    // cdata1.atoms.resize(numOfAtoms1);
    // cdata2.atoms.resize(numOfAtoms2);
}

ChargeTransfer::ChargeTransfer(Space &spc) : ChargeTransfer(spc, "chargetransfer", "") {}

void QuadrantJump::_to_json(json &j) const {
    j = {{"dir", dir},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molecules[molid].name}};
    _roundjson(j, 3);
}
void QuadrantJump::_from_json(const json &j) {
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
void QuadrantJump::_move(Change &change) {
    assert(molid >= 0);
    assert(!spc.groups.empty());
    assert(spc.geo.getVolume() > 0);

    _sqd = 0.0;

    // pick random group from the system matching molecule type
    // TODO: This can be slow -- implement look-up-table in Space
    auto mollist = spc.findMolecules(molid, Space::ACTIVE); // list of molecules w. 'molid'
    if (not ranges::cpp20::empty(mollist)) {
        auto it = slump.sample(mollist.begin(), mollist.end());
        if (not it->empty()) {
            assert(it->id == molid);
            Point oldcm = it->cm;
            if (index.size() == 2) {
                auto cm_O = Geometry::massCenter(spc.p.begin() + index[0], spc.p.begin() + index[1] + 1,
                                                 spc.geo.getBoundaryFunc());
                it->translate(-2 * spc.geo.vdist(oldcm, cm_O).cwiseProduct(dir.cast<double>()),
                              spc.geo.getBoundaryFunc());
            } else {
                it->translate(-2 * oldcm.cwiseProduct(dir.cast<double>()), spc.geo.getBoundaryFunc());
            }
            _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
            Change::data d;
            d.index = Faunus::distance(spc.groups.begin(), it); // integer *index* of moved group
            d.all = true;                                       // *all* atoms in group were moved
            change.groups.push_back(d);                         // add to list of moved groups

            assert(spc.geo.sqdist(it->cm, Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(),
                                                               -it->cm)) < 1e-9);
        }
    } else
        faunus_logger->warn("{0}: no molecules found", name);
}

QuadrantJump::QuadrantJump(Space &spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
}

QuadrantJump::QuadrantJump(Space &spc) : QuadrantJump(spc, "quadrantjump", "") {}

void AtomicSwapCharge::_to_json(json &j) const {
    j = {{"pH", pH},
         {"pka", pKa},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molname}};
    _roundjson(j, 3);
}
void AtomicSwapCharge::_from_json(const json &j) {
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
typename Space::Tpvec::iterator AtomicSwapCharge::randomAtom() {
    assert(molid >= 0);
    auto mollist = spc.findMolecules(molid); // all `molid` groups
    if (not ranges::cpp20::empty(mollist)) {
        auto git = slump.sample(mollist.begin(), mollist.end()); // random molecule iterator
        if (!git->empty()) {
            auto p = slump.sample(git->begin(), git->end());         // random particle iterator
            cdata.index = Faunus::distance(spc.groups.begin(), git); // integer *index* of moved group
            cdata.atoms[0] = std::distance(git->begin(), p);         // index of particle rel. to group
            return p;
        }
    }
    return spc.p.end();
}
void AtomicSwapCharge::_move(Change &change) {
    _sqd = 0.0;
    auto p = randomAtom();
    if (p != spc.p.end()) {
        // auto &g = spc.groups[cdata.index];
        double oldcharge = p->charge;
        p->charge = fabs(oldcharge - 1);
        _sqd = fabs(oldcharge - 1) - oldcharge;
        change.groups.push_back(cdata);   // add to list of moved groups
        _bias = _sqd * (pH - pKa) * ln10; // one may add bias here...
    }
}
double AtomicSwapCharge::bias(Change &, double, double) { return _bias; }
void AtomicSwapCharge::_accept(Change &) { msqd += _sqd; }
void AtomicSwapCharge::_reject(Change &) { msqd += 0; }

AtomicSwapCharge::AtomicSwapCharge(Space &spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
    cdata.atoms.resize(1);
    cdata.internal = true;
}

AtomicSwapCharge::AtomicSwapCharge(Space &spc) : AtomicSwapCharge(spc, "swapcharge", "") {}

void TranslateRotate::_to_json(json &j) const {
    j = {{"dir", translational_direction},
         {"dp", translational_displacement},
         {"dprot", rotational_displacement},
         {"dirrot", fixed_rotation_axis},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(mean_squared_displacement.avg())},
         {"√⟨θ²⟩/°", std::sqrt(mean_squared_rotation_angle.avg()) / 1.0_deg},
         {"molecule", Faunus::molecules[molid].name}};
    _roundjson(j, 3);
}
void TranslateRotate::_from_json(const json &j) {
    const std::string molname = j.at("molecule");
    const auto molecule = findMoleculeByName(molname);
    if(molecule.atomic) {
        throw ConfigurationError("molecule '{}' cannot be atomic", molname);
    }
    molid = molecule.id();
    translational_direction = j.value("dir", Point(1, 1, 1));
    translational_displacement = j.at("dp").get<double>();
    rotational_displacement = j.at("dprot").get<double>();
    fixed_rotation_axis = j.value("dirrot", Point(0, 0, 0)); // predefined axis of rotation
    if (repeat < 0) {
        repeat = spc.numMolecules<Space::Tgroup::ACTIVE>(molid);
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
std::optional<std::reference_wrapper<Space::Tgroup>> TranslateRotate::findRandomMolecule() const {
    if (auto mollist = spc.findMolecules(molid, Space::ACTIVE); not ranges::cpp20::empty(mollist)) {
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
double TranslateRotate::translateMolecule(Space::Tgroup &group) {
    if (translational_displacement > 0.0) { // translate
        const auto old_mass_center = group.cm;
        const auto displacement_vector =
            Faunus::ranunit(slump, translational_direction) * translational_displacement * slump();

        group.translate(displacement_vector, spc.geo.getBoundaryFunc());
        return spc.geo.sqdist(old_mass_center, group.cm);
    } else {
        return 0.0;
    }
}

/**
 * @param group Group to rotate
 * @return Squared rotation angle around mass-center
 */
double TranslateRotate::rotateMolecule(Space::Tgroup &group) {
    if (rotational_displacement > 0.0) {         // rotate
        Point rotation_axis = Faunus::ranunit(slump); // rotate around random unit vector
        if (fixed_rotation_axis.count() > 0.0) { // OR a fixed, user-defined vector
            rotation_axis = fixed_rotation_axis;
        }
        const auto angle = rotational_displacement * (slump() - 0.5); // random rotation angle
        const Eigen::Quaterniond quaternion(Eigen::AngleAxisd(angle, rotation_axis));
        group.rotate(quaternion, spc.geo.getBoundaryFunc());
        return angle * angle;
    } else {
        return 0.0;
    }
}

void TranslateRotate::_move(Change &change) {
    if (auto group = findRandomMolecule()) { // note that group is of type std::optional
        latest_displacement_squared = translateMolecule(group->get());
        latest_rotation_angle_squared = rotateMolecule(group->get());
        if (latest_displacement_squared > 0.0 || latest_rotation_angle_squared > 0.0) { // report changes
            auto &change_data = change.groups.emplace_back();
            change_data.index = spc.getGroupIndex(group->get());     // integer *index* of moved group
            change_data.all = true;                                  // *all* atoms in group were moved
            change_data.internal = false;                            // internal energy is unchanged
        }
#ifndef NDEBUG
        sanityCheck(group->get());
#endif
    } else {
        latest_displacement_squared = 0.0;   // these are used to track mean squared
        latest_rotation_angle_squared = 0.0; // translational and rotational displacements
    }
}

void TranslateRotate::sanityCheck(const Space::Tgroup &group) const {
    Point cm_recalculated = Geometry::massCenter(group.begin(), group.end(), spc.geo.getBoundaryFunc(), -group.cm);
    double should_be_small = spc.geo.sqdist(group.cm, cm_recalculated);
    if (should_be_small > 1e-6) {
        std::cerr << "cm recalculated: " << cm_recalculated.transpose() << "\n";
        std::cerr << "cm in group:     " << group.cm.transpose() << "\n";
        assert(false);
    }
}

void TranslateRotate::_accept(Change &) {
    mean_squared_displacement += latest_displacement_squared;
    mean_squared_rotation_angle += latest_rotation_angle_squared;
}

void TranslateRotate::_reject(Change &) {
    mean_squared_displacement += 0.0;
    mean_squared_rotation_angle += 0.0;
}

TranslateRotate::TranslateRotate(Space &spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
}

TranslateRotate::TranslateRotate(Space &spc) : TranslateRotate(spc, "moltransrot", "") {}
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

    j = json(mv).at(mv.name);
    CHECK(j.at("molecule") == "A");
    // CHECK(j.at("dir") == Point(0, 1, 0));
    CHECK(j.at("dp") == 1.0);
    CHECK(j.at("repeat") == 2);
    CHECK(j.at("dprot") == 0.5);
}
#endif

namespace Faunus::Move {

void SmartTranslateRotate::_to_json(json &j) const {
    j = {{"Number of counts inside geometry", cntInner},
         {"Number of counts outside geometry", cnt - cntInner},
         {"dir", dir},
         {"dp", dptrans},
         {"dprot", dprot},
         {"p", p},
         {"origo", origo},
         {"rx", r_x},
         {"ry", r_y},
         {"molid", molid},
         {"refid", refid1},
         {"refid", refid2},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molecules[molid].name},
         {"ref1", atoms[refid1].name},
         {"ref2", atoms[refid2].name}};
    _roundjson(j, 3);
}
void SmartTranslateRotate::_from_json(const json &j) {
    assert(!molecules.empty());
    const std::string molname = j.at("molecule");
    const std::string refname1 = j.at("ref1");
    const std::string refname2 = j.at("ref2");
    molid = findMoleculeByName(molname).id();
    refid1 = findAtomByName(refname1).id();
    refid2 = findAtomByName(refname2).id();;
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

void SmartTranslateRotate::_move(Change &change) {
    assert(molid >= 0);
    assert(!spc.groups.empty());
    assert(spc.geo.getVolume() > 0);
    _bias = 0.0;
    _sqd = 0.0;

    // pick random group from the system matching molecule type
    // TODO: This can be slow -- implement look-up-table in Space
    auto mollist = spc.findMolecules(molid, Space::ACTIVE); // list of molecules w. 'molid'
    auto reflist1 = spc.findAtoms(refid1);                  // list of atoms w. 'refid1'
    auto reflist2 = spc.findAtoms(refid2);                  // list of atoms w. 'refid2'
    if (not ranges::cpp20::empty(mollist)) {
        auto it = slump.sample(mollist.begin(), mollist.end()); // chosing random molecule in group of type molname
        auto ref1 = slump.sample(reflist1.begin(), reflist1.end());
        auto ref2 = slump.sample(reflist2.begin(), reflist2.end());
        cylAxis = spc.geo.vdist(ref2->pos, ref1->pos) * 0.5; // half vector between reference atoms
        origo = ref2->pos - cylAxis; // coordinates of middle point between reference atoms: new origo
        if (r_x < cylAxis.norm())    // checking so that a is larger than length of cylAxis
            throw std::runtime_error(
                "specified radius of ellipsoid along the axis connecting reference atoms (rx) must be larger or equal "
                "to half the distance between reference atoms. Specified radius is " +
                std::to_string(r_x) + " Å whereas half the distance between reference atoms is " +
                std::to_string(cylAxis.norm()) + "Å");

        if (not it->empty()) { // checking so that molecule exists
            assert(it->id == molid);

            randNbr = slump();                   // assigning random number in range [0,1]
            molV = spc.geo.vdist(it->cm, origo); // vector between selected molecule and center of geometry
            cosTheta = molV.dot(cylAxis) / molV.norm() / cylAxis.norm(); // cosinus of angle between coordinate vector
                                                                         // of selected molecule and axis connecting
                                                                         // reference atoms
            theta = acos(cosTheta);       // angle between coordinate vector of selected molecule and axis connecting
                                          // reference atoms
            x = cosTheta * molV.norm();   // x coordinate of selected molecule with respect to center of geometry
                                          // (in plane including vectors molV and cylAxis)
            y = sin(theta) * molV.norm(); // y coordinate of selected molecule with respect to center of geometry
                                          // (in plane including vectors molV and cylAxis)
            coord = x * x / (r_x * r_x) + y * y / (r_y * r_y); // calculating normalized coordinate with respect to
                                                               // dimensions of geometry (>1.0 → outside, <1.0 → inside)
            if (not(coord > 1.0 && p < randNbr)) {

                if (coord <= 1.0)
                    cntInner += 1; // counting number of times a molecule is found inside geometry

                cnt += 1; // total number of counts

                if (findBias ==
                    true) { // continuing to adjust bias according to number of molecules inside and outside geometry
                    countNin = 0.0;  // counter keeping track of number of molecules inside geometry
                    countNout = 0.0; // counter keeping track of number of molecules outside geometry
                    Ntot = 0.0;      // total number of particles
                    for (auto &g : mollist) {
                        Ntot += 1.0;
                        molV = spc.geo.vdist(g.cm, origo);
                        cosTheta = molV.dot(cylAxis) / molV.norm() / cylAxis.norm();
                        theta = acos(cosTheta);
                        x = cosTheta * molV.norm();
                        y = sin(theta) * molV.norm();
                        coordTemp = x * x / (r_x * r_x) + y * y / (r_y * r_y);
                        if (coordTemp <= 1.0)
                            countNin += 1.0;
                        else
                            countNout += 1.0;
                    }

                    countNin_avg += countNin;   // appending number of molecules inside geometry
                                                // (since it has type Average)
                    countNout_avg += countNout; // appending number of molecules outside geometry
                                                // (since it has type Average)

                    if (cnt % 100 == 0) {
                        countNin_avgBlocks += countNin_avg.avg();   // appending average number of molecules inside
                                                                    // geometry (type Average)
                        countNout_avgBlocks += countNout_avg.avg(); // appending average number of molecules outside
                                                                    // geometry (type Average)
                    }

                    if (cnt % 100000 == 0) {
                        Nin = countNin_avgBlocks.avg(); // block average number of molecules inside geometry
                        if (countNin_avgBlocks.stdev() / Nin < rsd) {
                            // if block standard deviation is below specified threshold
                            std::cout << "Bias found with rsd = " << countNin_avgBlocks.stdev() / Nin << " < " << rsd
                                      << "\n\n";
                            std::cout << "Average # of water molecules inside sphere: " << Nin << "\n";
                            findBias = false; // stop updating bias, use constant value
                        }
                    }
                }
                if (dptrans > 0) { // translate
                    Point oldcm = it->cm;
                    Point dp = ranunit(slump, dir) * dptrans * slump();

                    it->translate(dp, spc.geo.getBoundaryFunc());
                    _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
                }

                if (dprot > 0) { // rotate
                    Point u = ranunit(slump);
                    double angle = dprot * (slump() - 0.5);
                    Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, u));
                    it->rotate(Q, spc.geo.getBoundaryFunc());
                }

                if (dptrans > 0 || dprot > 0) { // define changes
                    Change::data d;
                    d.index = Faunus::distance(spc.groups.begin(), it); // integer *index* of moved group
                    d.all = true;                                       // *all* atoms in group were moved
                    change.groups.push_back(d);                         // add to list of moved groups
                }
                assert(spc.geo.sqdist(it->cm, Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(),
                                                                   -it->cm)) < 1e-6);
                molV = spc.geo.vdist(it->cm, origo);
                cosTheta = molV.dot(cylAxis) / molV.norm() / cylAxis.norm();
                theta = acos(cosTheta);
                x = cosTheta * molV.norm();
                y = sin(theta) * molV.norm();
                coordNew = x * x / (r_x * r_x) + y * y / (r_y * r_y);

                if (findBias == true) {                 // if using constantly updated bias
                    if (coord <= 1.0 && coordNew > 1.0) // if molecule goes from inside to outside geometry
                        // use corresponding bias, based on this cycle's number of molecules inside
                        _bias = -log(p / (1 - (1 - p) / (p * Ntot + (1 - p) * countNin)));
                    else if (coord > 1.0 && coordNew <= 1.0) // if molecules goes from outside to inside geometry
                        // use corresponding bias, based on this cycle's number of molecules inside
                        _bias = -log(1 / (1 + (1 - p) / (p * Ntot + (1 - p) * countNin)));
                }

                else {                                  // if constant bias has been found
                    if (coord <= 1.0 && coordNew > 1.0) // if molecule goes from inside to outside geometry
                        // use corresponding bias based on average, constant value Nin
                        _bias = -log(p / (1 - (1 - p) / (p * Ntot + (1 - p) * Nin)));
                    else if (coord > 1.0 && coordNew <= 1.0) // if molecule goes from outside to inside geometry
                        // use corresponding bias based on average, constant value Nin
                        _bias = -log(1 / (1 + (1 - p) / (p * Ntot + (1 - p) * Nin)));
                }
            }
        }
    }
}

double SmartTranslateRotate::bias(Change &, double, double) { return _bias; }

SmartTranslateRotate::SmartTranslateRotate(Space &spc, std::string name, std::string cite) : MoveBase(spc, name, cite) {
    repeat = -1; // meaning repeat N times
}

SmartTranslateRotate::SmartTranslateRotate(Space &spc) : SmartTranslateRotate(spc, "smartmoltransrot", "") {}

void ConformationSwap::_to_json(json &j) const {
    j = {{"molid", molid}, {"molecule", Faunus::molecules.at(molid).name}, {"keeppos", inserter.keep_positions}};
    _roundjson(j, 3);
}
void ConformationSwap::_from_json(const json &j) {
    const std::string molecule_name = j.at("molecule");
    const auto molecule = findMoleculeByName(molecule_name);
    if (molecule.conformations.size() < 2) {
        throw ConfigurationError("molecule '{}': minimum two conformations required", molecule_name);
    }
    molid = molecule.id();
    inserter.keep_positions = j.value("keeppos", false);
    setRepeat();
}

void ConformationSwap::setRepeat() {
    assert(molid >= 0);
    if (repeat < 0) { // negative value signals repeat = N number of molecules
        auto groups = spc.findMolecules(molid, Space::ALL);
        repeat = std::distance(groups.begin(), groups.end());
        if (repeat == 0) {
            faunus_logger->warn("{}: no molecules found; repeat set to ZERO", name, repeat);
        }
    }
}
void ConformationSwap::_move(Change &change) {
    if (auto groups = spc.findMolecules(molid, Space::ACTIVE); !ranges::cpp20::empty(groups)) {
        if (auto &group = *slump.sample(groups.begin(), groups.end()); !group.empty()) { // pick random molecule
            inserter.offset = group.cm;                                                  // insert on top of mass center
            const auto particles = inserter(spc.geo, Faunus::molecules[molid], spc.p);   // new conformation
            if (particles.size() == group.size()) {
                checkMassCenterDrift(group.cm, particles);                            // throws if not OK
                std::copy(particles.begin(), particles.end(), group.begin());         // override w. new conformation
                group.confid = Faunus::molecules[molid].conformations.getLastIndex(); // store conformation id
                registerChanges(change, group);                                       // update change object
            } else {
                throw std::runtime_error(name + ": conformation atom count mismatch");
            }
        }
    }
}
void ConformationSwap::registerChanges(Change &change, const Space::Tgroup &group) const {
    auto &group_change = change.groups.emplace_back();
    group_change.index = spc.getGroupIndex(group); // index of moved group
    group_change.all = true;                       // all atoms in group were moved
    group_change.internal = false;                 // skip internal energy calculation
}
/**
 * @throw if there's a mass-center drift
 *
 * Move shouldn't move mass centers, so let's check if this is true
 */
void ConformationSwap::checkMassCenterDrift(const Point &old_mass_center, const ParticleVector &particles) {
    const auto max_allowed_distance = 1.0e-6;
    const auto new_mass_center =
        Geometry::massCenter(particles.begin(), particles.end(), spc.geo.getBoundaryFunc(), -old_mass_center);
    if ((new_mass_center - old_mass_center).norm() > max_allowed_distance) {
        throw std::runtime_error(name + ": unexpected mass center movement");
    }
}

ConformationSwap::ConformationSwap(Space &spc, const std::string &name, const std::string &cite)
    : MoveBase(spc, name, cite) {}

ConformationSwap::ConformationSwap(Space &spc) : ConformationSwap(spc, "conformationswap", "") {
    repeat = -1; // meaning repeat n times
    inserter.dir = {0, 0, 0};
    inserter.rotate = true;
    inserter.allow_overlap = true;
}

} // namespace Faunus::Move
