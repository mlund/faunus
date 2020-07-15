#include "montecarlo.h"
#include "speciation.h"
#include "energy.h"
#include "move.h"
#include "spdlog/spdlog.h"

namespace Faunus {

/**
 * @param du Energy change in units of kT
 * @return True if accepted, false of rejected
 */
bool MetropolisMonteCarlo::metropolis(double du) const {
    if (std::isnan(du)) {
        throw std::runtime_error("Metropolis error: energy cannot be NaN");
    }
    if (du < 0) {
        return true;
    } else {
        if (-du > pc::max_exp_argument) {
            mcloop_logger->warn("large negative metropolis energy");
        }
        return Move::Movebase::slump() <= std::exp(-du);
    }
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
    Change change;
    change.all = true;

    state->pot->key = Energy::Energybase::ACCEPTED_MONTE_CARLO_STATE;    // this is the old energy (current, accepted)
    trial_state->pot->key = Energy::Energybase::TRIAL_MONTE_CARLO_STATE; // this is the new energy (trial)

    state->pot->init();
    double energy = state->pot->energy(change);
    initial_energy = energy;

    trial_state->sync(*state, change); // copy all information into trial state
    trial_state->pot->init();
    double trial_energy = trial_state->pot->energy(change);

    // check that the energies in the two states are *identical*
    if (std::isfinite(energy) and std::isfinite(trial_energy)) {
        if (std::fabs((energy - trial_energy) / energy) > 1e-6) {
            faunus_logger->error("u/kT = {}, u_trial/kT = {}", energy, trial_energy);
            throw std::runtime_error("error aligning energies - this could be a bug...");
        }
    }

    // Inject reference to Space into `SpeciationMove`
    // Needed to calc. differences in ideal excess chem. potentials
    for (auto speciation_move : moves->moves().find<Move::SpeciationMove>()) {
        speciation_move->setOther(*state->spc);
    }
}

/**
 * This is a powerful test that compares the final energy
 * with the initial energy plus the accumulated list of
 * changes. Ideally this difference should be very small.
 * If a MC move improperly captures the resulting energy change,
 * a large drift will usually appear. This test has a minimal
 * computational overhead as only the initial and final energies
 * need to be calculated; the accumulated change is sampled during
 * simulation for free.
 *
 * @return Relative energy drift
 */
double MetropolisMonteCarlo::relativeEnergyDrift() {
    Change change;     // change object where ...
    change.all = true; // ... we want to calculate the total energy
    double energy = state->pot->energy(change);
    double du = energy - initial_energy;
    if (std::isfinite(du)) {
        if (std::fabs(du) <= pc::epsilon_dbl) {
            return 0.0;
        } else {
            return (energy - (initial_energy + sum_of_energy_changes)) /
                   (initial_energy != 0 ? initial_energy : energy);
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

MetropolisMonteCarlo::MetropolisMonteCarlo(const json &j, MPI::MPIController &mpi)
    : original_log_level(faunus_logger->level()) {
    state = std::make_shared<State>(j);
    faunus_logger->set_level(spdlog::level::off); // do not duplicate log info
    trial_state = std::make_shared<State>(j);     // ...for the trial state
    faunus_logger->set_level(original_log_level); // restore original log level
    moves = std::make_shared<Move::Propagator>(j, *trial_state->spc, *trial_state->pot, mpi);
    init();
}

/**
 * @todo Too many responsibilities; tidy up!
 */
void MetropolisMonteCarlo::restore(const json &j) {
    try {
        *state->spc = j;       // default, accepted state
        *trial_state->spc = j; // trial state
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

/**
 * This propagates the system using a random MC move.
 * Flow:
 *
 * 1. Pick random MC move
 * 2. Perform move
 * 3. Calculate potential energy change
 * 4. Add bias from move if appropriate
 * 5. Add bias from density fluctuations (grand canonical, speciation)
 * 6. Use Metropolis criterion to either accept or reject move
 * 7. If accepted, copy changes from trial state to state
 * 8. If rejected, restored changes from state to trial state
 * 9. After the move the state and the trial state should be identical!
 */
void MetropolisMonteCarlo::move() {
    assert(moves);
    for (int i = 0; i < moves->repeat(); i++) {
        if (auto move_it = moves->sample(); move_it != moves->end()) { // pick random move
            Change change;                                             // stores proposed changes due to move
            auto move = *move_it;                                      // more readable like this
            move->move(change);
#ifndef NDEBUG
            // check if atom index indeed belong to the group (index)
            if (not change.sanityCheck(*state->spc)) {
                throw std::runtime_error("insane change object\n" + json(change).dump(4));
            }
#endif
            if (change) {
                latest_move = move;
                double trial_energy = trial_state->pot->energy(change);    // trial potential energy (kT)
                double energy = state->pot->energy(change);                // potential energy before move (kT)
                double du = trial_energy - energy;                         // potential energy change (kT)
                if (std::isnan(energy) and not std::isnan(trial_energy)) { // if NaN --> finite energy change
                    du = pc::neg_infty;                                    // ...always accept
                } else if (std::isnan(trial_energy)) {                     // if moving to NaN, e.g. division by zero,
                    du = pc::infty;                                        // ...always reject
                } else if (std::isnan(du)) { // if difference is NaN, e.g. infinity - infinity,
                    du = 0.0;                // ...always accept
                }
                double move_bias = move->bias(change, energy, trial_energy); // moves *may* add bias (kT)
                double density_bias = TranslationalEntropy(*trial_state->spc, *state->spc).energy(change);
                if (std::isnan(du + move_bias)) {
                    faunus_logger->error("NaN energy change in {} move.", move->name);
                    // throw exception here?
                }
                if (metropolis(du + move_bias + density_bias)) { // accept move
                    state->sync(*trial_state, change);
                    move->accept(change);
                } else { // reject move
                    trial_state->sync(*state, change);
                    move->reject(change);
                    du = 0.0;
                }
                sum_of_energy_changes += du;                              // sum of all energy changes
                if (std::isfinite(initial_energy)) {
                    average_energy += initial_energy + sum_of_energy_changes; // update average potential energy
                }
            }
        }
    }
}

Energy::Hamiltonian &MetropolisMonteCarlo::getHamiltonian() { return *state->pot; }

Space &MetropolisMonteCarlo::getSpace() { return *state->spc; }

void from_json(const json &j, MetropolisMonteCarlo::State &state) {
    state.spc = std::make_shared<Space>(j);
    state.pot = std::make_shared<Energy::Hamiltonian>(*state.spc, j.at("energy"));
}

void MetropolisMonteCarlo::State::sync(MetropolisMonteCarlo::State &other, Change &change) {
    spc->sync(*other.spc, change);
    pot->sync(&*other.pot, change);
}

void to_json(json &j, const MetropolisMonteCarlo &mc) {
    j = mc.state->spc->info();
    j["temperature"] = pc::temperature / 1.0_K;
    j["moves"] = *mc.moves;
    j["energy"].push_back(*mc.state->pot);
    if (mc.average_energy.cnt > 0) {
        j["montecarlo"] = {{"average potential energy (kT)", mc.average_energy.avg()},
                           {"last move", mc.latest_move->name}};
    }
}

TranslationalEntropy::TranslationalEntropy(Space &trial_space, Space &space) : trial_spc(trial_space), spc(space) {}

/**
 * @param trial_count Number of atoms or molecules after move
 * @param count Number of atoms or molecular before move
 * @return Energy contribution (kT) to be added to MC trial energy
 */
double TranslationalEntropy::bias(int trial_count, int count) const {
    double energy = 0.0;
    if (int dN = trial_count - count; dN > 0) { // atoms or molecules were added
        double V_trial = trial_spc.geo.getVolume();
        for (int n = 0; n < dN; n++) {
            energy += std::log((count + 1 + n) / (V_trial * 1.0_molar));
        }
    } else if (dN < 0) { // atoms or molecules were removed
        double V = spc.geo.getVolume();
        for (int n = 0; n < (-dN); n++) {
            energy -= std::log((count - n) / (V * 1.0_molar));
        }
    }
    return energy; // kT
}

double TranslationalEntropy::atomSwapEnergy(const Change::data &data) {
    assert(data.dNswap);
    assert(data.atoms.size() == 1);
    double energy = 0.0;
    int id1 = trial_spc.groups[data.index][data.atoms.front()].id;
    int id2 = spc.groups[data.index][data.atoms.front()].id;
    for (auto atomid : {id1, id2}) {
        auto atoms_new = trial_spc.findAtoms(atomid);
        auto atoms_old = spc.findAtoms(atomid);
        int N_new = range_size(atoms_new); // number of atoms after change
        int N_old = range_size(atoms_old); // number of atoms before change
        energy += bias(N_new, N_old);
    }
    return energy; // kT
}

double TranslationalEntropy::atomChangeEnergy(int molid) {
    auto mollist_new = trial_spc.findMolecules(molid, Space::ALL); // "ALL" because "ACTIVE"
    auto mollist_old = spc.findMolecules(molid, Space::ALL);       // ...returns only full groups
    if (range_size(mollist_new) > 1 || range_size(mollist_old) > 1) {
        throw std::runtime_error("multiple atomic groups of the same type is not allowed");
    }
    int N_new = mollist_new.begin()->size(); // number of atoms after move
    int N_old = mollist_old.begin()->size(); // number of atoms before move
    return bias(N_new, N_old);
}

double TranslationalEntropy::moleculeChangeEnergy(int molid) {
    auto mollist_new = trial_spc.findMolecules(molid, Space::ACTIVE);
    auto mollist_old = spc.findMolecules(molid, Space::ACTIVE);
    int N_new = range_size(mollist_new); // number of molecules after move
    int N_old = range_size(mollist_old); // number of molecules before move
    return bias(N_new, N_old);
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
                int molid = trial_spc.groups.at(data.index).id;
                assert(molid == spc.groups.at(data.index).id);
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
