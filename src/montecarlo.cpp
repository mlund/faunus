#include "montecarlo.h"
#include "speciation.h"
#include "spdlog/spdlog.h"

namespace Faunus {

/**
 * @param du Energy change in units of kT
 * @return True if accepted, false of rejected
 */
bool MetropolisMonteCarlo::metropolis(double du) const {
    if (std::isnan(du)) {
        throw std::runtime_error("Metropolis error: energy cannot be NaN");
    } else if (du < 0) {
        return true;
    } else if (-du > pc::max_exp_argument) {
        mcloop_logger->warn("large metropolis energy");
    }
    return Move::Movebase::slump() <= std::exp(-du);
}

/**
 * This performas the following tasks:
 * - syncs the two states
 * - syncs the two Hamiltonians
 * - resets the sum of energy changes to zero
 * - recalculates the initial energy
 */
void MetropolisMonteCarlo::init() {
    sum_of_energy_changes = 0;
    Change c;
    c.all = true;

    old_state.pot.key = Energy::Energybase::OLD_MONTE_CARLO_STATE; // this is the old energy (current, accepted)
    new_state.pot.key = Energy::Energybase::NEW_MONTE_CARLO_STATE; // this is the new energy (trial)

    old_state.pot.init();
    double u1 = old_state.pot.energy(c);
    initial_energy = u1;

    new_state.sync(old_state, c); // copy all information from state1 into state2
    new_state.pot.init();
    double u2 = new_state.pot.energy(c);

    // check that the energies in the two states are *identical*
    if (std::isfinite(u1) and std::isfinite(u2)) {
        if (std::fabs((u1 - u2) / u1) > 1e-6) {
            faunus_logger->error("u_old = {}, u_new = {}", u1, u2);
            throw std::runtime_error("error aligning energies - this could be a bug...");
        }
    }

    // inject reference to state1 in SpeciationMove (needed to calc. *differences*
    // in ideal excess chem. potentials)
    for (auto speciation_move : moves.moves().find<Move::SpeciationMove>()) {
        speciation_move->setOther(old_state.spc);
    }
}

double MetropolisMonteCarlo::relativeEnergyDrift() {
    Change c;
    c.all = true;
    double final_energy = old_state.pot.energy(c);
    double du = final_energy - initial_energy;
    if (std::isfinite(du)) {
        if (std::fabs(du) <= pc::epsilon_dbl) {
            return 0.0;
        } else if (initial_energy != 0) {
            return (final_energy - (initial_energy + sum_of_energy_changes)) / initial_energy;
        } else if (final_energy != 0)
            return (final_energy - (initial_energy + sum_of_energy_changes)) / final_energy;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

/*
 * We need to construct two identical State objects and to avoid duplicate logs, we
 * temporarily disable the logger for the second object by the arcane _comma operator_
 */
MetropolisMonteCarlo::MetropolisMonteCarlo(const json &j, MPI::MPIController &mpi)
    : original_log_level(faunus_logger->level()), old_state(j),
      new_state((faunus_logger->set_level(spdlog::level::off), j)),
      moves((faunus_logger->set_level(original_log_level), j), new_state.spc, new_state.pot, mpi) {
    init();
}

void MetropolisMonteCarlo::restore(const json &j) {
    try {
        old_state.spc = j; // old/accepted state
        new_state.spc = j; // trial state
        if (j.count("random-move") == 1) {
            Move::Movebase::slump = j["random-move"]; // restore move random number generator
        }
        if (j.count("random-global") == 1) {
            Faunus::random = j["random-global"]; // restore global random number generator
        }
        reactions = j.at("reactionlist").get<decltype(reactions)>(); // should be handled by space
        init();
    } catch (std::exception &e) {
        throw std::runtime_error("error initialising simulation: "s + e.what());
    }
}

void MetropolisMonteCarlo::move() {
    for (int i = 0; i < moves.repeat(); i++) {
        if (auto move = moves.sample(); move != moves.end()) { // pick random move
            Change change;                                     // stores proposed changes due to move
            (**move).move(change);
#ifndef NDEBUG
            // check if atom index indeed belong to the group (index)
            if (not change.sanityCheck(old_state.spc)) {
                throw std::runtime_error("insane change object\n" + json(change).dump(4));
            }
#endif
            if (change) {
                latest_move = *move;
                double unew = new_state.pot.energy(change);      // old potential energy (kT)
                double uold = old_state.pot.energy(change);      // new potential energy (kT)
                double du = unew - uold;                         // potential energy change (kT)
                if (std::isnan(uold) and not std::isnan(unew)) { // if NaN --> finite energy change
                    du = pc::neg_infty;                          // ...always accept
                } else if (std::isnan(unew)) {                   // if moving to NaN, e.g. division by zero,
                    du = pc::infty;                              // ...always reject
                } else if (std::isnan(du)) {                     // if difference is NaN, e.g. infinity - infinity,
                    du = 0.0;                                    // ...always accept
                }
                double bias = (*move)->bias(change, uold, unew); // moves *may* add bias (kT)
                double ideal = TranslationalEntropy(new_state.spc, old_state.spc).energy(change);
                if (std::isnan(du + bias)) {
                    faunus_logger->error("NaN energy du + bias in {} move.", (*move)->name);
                    // throw here?
                }
                if (metropolis(du + bias + ideal)) { // accept move
                    old_state.sync(new_state, change);
                    (*move)->accept(change);
                } else { // reject move
                    new_state.sync(old_state, change);
                    (*move)->reject(change);
                    du = 0.0;
                }
                sum_of_energy_changes += du; // sum of all energy changes
            }
        }
    }
}

void MetropolisMonteCarlo::to_json(json &j) {
    j = old_state.spc.info();
    j["temperature"] = pc::temperature / 1.0_K;
    j["moves"] = moves;
    j["energy"].push_back(old_state.pot);
    j["last move"] = latest_move->name;
}
Energy::Hamiltonian &MetropolisMonteCarlo::pot() { return old_state.pot; }
const Energy::Hamiltonian &MetropolisMonteCarlo::pot() const { return old_state.pot; }

Space &MetropolisMonteCarlo::space() { return old_state.spc; }
const Space &MetropolisMonteCarlo::space() const { return old_state.spc; }

const Space::Tgeometry &MetropolisMonteCarlo::getGeometry() const { return old_state.spc.geo; }

MetropolisMonteCarlo::State::State(const json &j) : spc(j), pot(spc, j.at("energy")) {}

void MetropolisMonteCarlo::State::sync(MetropolisMonteCarlo::State &other, Change &change) {
    spc.sync(other.spc, change);
    pot.sync(&other.pot, change);
}

void to_json(json &j, MetropolisMonteCarlo &mc) { mc.to_json(j); }

TranslationalEntropy::TranslationalEntropy(Space &new_space, Space &old_space)
    : spc_new(new_space), spc_old(old_space) {}

/**
 * @param N_new Number of atoms or molecules after move
 * @param N_old Number of atoms or molecular before move
 * @return Energy contribution (kT) to be added to MC trial energy
 */
double TranslationalEntropy::accumulate(int N_new, int N_old) const {
    double energy = 0.0;
    if (int dN = N_new - N_old; dN > 0) { // atoms or molecules were added
        double V_new = spc_new.geo.getVolume();
        for (int n = 0; n < dN; n++) {
            energy += std::log((N_old + 1 + n) / (V_new * 1.0_molar));
        }
    } else if (dN < 0) { // atoms or molecules were removed
        double V_old = spc_old.geo.getVolume();
        for (int n = 0; n < (-dN); n++) {
            energy -= std::log((N_old - n) / (V_old * 1.0_molar));
        }
    }
    return energy; // kT
}

double TranslationalEntropy::atomSwapEnergy(const Change::data &data) {
    assert(data.dNswap);
    assert(data.atoms.size() == 1);
    double energy = 0.0;
    int id1 = spc_new.groups[data.index][data.atoms.front()].id;
    int id2 = spc_old.groups[data.index][data.atoms.front()].id;
    for (int atomid : {id1, id2}) {
        auto atoms_new = spc_new.findAtoms(atomid);
        auto atoms_old = spc_old.findAtoms(atomid);
        int N_new = range_size(atoms_new); // number of atoms after change
        int N_old = range_size(atoms_old); // number of atoms before change
        energy += accumulate(N_new, N_old);
    }
    return energy; // kT
}

double TranslationalEntropy::atomChangeEnergy(int molid) {
    auto mollist_new = spc_new.findMolecules(molid, Space::ALL); // "ALL" because "ACTIVE"
    auto mollist_old = spc_old.findMolecules(molid, Space::ALL); // ...returns only full groups
    if (range_size(mollist_new) > 1 || range_size(mollist_old) > 1) {
        throw std::runtime_error("multiple atomic groups of the same type is not allowed");
    }
    int N_new = mollist_new.begin()->size(); // number of atoms after move
    int N_old = mollist_old.begin()->size(); // number of atoms before move
    return accumulate(N_new, N_old);
}

double TranslationalEntropy::moleculeChangeEnergy(int molid) {
    auto mollist_new = spc_new.findMolecules(molid, Space::ACTIVE);
    auto mollist_old = spc_old.findMolecules(molid, Space::ACTIVE);
    int N_new = range_size(mollist_new); // number of molecules after move
    int N_old = range_size(mollist_old); // number of molecules before move
    return accumulate(N_new, N_old);
}

/**
 * @param change Change due to latest Monte Carlo move
 * @return Logarithm of the bias for the Metropolis criterion (units of kT)
 */
double TranslationalEntropy::energy(const Change &change) {
    double energy_change = 0.0;
    if (change.dN) {
        std::set<int> already_processed;                 // ignore future encounters of these molid's
        for (const Change::data &data : change.groups) { // loop over each change group
            if (data.dNswap) {                           // number of atoms has changed as a result of a swap move
                energy_change += atomSwapEnergy(data);
            } else { // it is not a swap move
                int molid = spc_new.groups.at(data.index).id;
                assert(molid == spc_old.groups.at(data.index).id);
                if (data.dNatomic and Faunus::molecules[molid].atomic) { // an atomic group has been changed
                    energy_change += atomChangeEnergy(molid);
                } else if (already_processed.count(molid) == 0) { // a molecule has been inserted or deleted
                    energy_change += moleculeChangeEnergy(molid);
                    already_processed.insert(molid); // ignore future encounters of molid
                }
            }
        }
    }
    return energy_change; // kT
}
} // namespace Faunus
