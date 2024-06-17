#include "montecarlo.h"
#include "speciation.h"
#include "energy.h"
#include "move.h"
#include <spdlog/spdlog.h>
#include <range/v3/algorithm/for_each.hpp>

namespace Faunus {

/**
 * @param energy_change Energy change, (new minus old) in units of kT
 * @return True if accepted, false of rejected
 * @note Regardless of outcome the random number generator should be incremented. This is important
 *       when using some MPI schemes where the simulations must be in sync.
 */
bool MetropolisMonteCarlo::metropolisCriterion(const double energy_change) {
    static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");
    if (std::isnan(energy_change)) {
        throw std::runtime_error("Metropolis error: energy cannot be NaN");
    }
    const auto random_number_between_zero_and_one = move::Move::slump();     // engine *must* be propagated!
    if (std::isinf(energy_change) && energy_change < 0.0) {                  // if negative infinity -> quietly accept
        return true;
    }
    if (-energy_change > pc::max_exp_argument) { // if large negative value -> accept with warning
        mcloop_logger->warn("humongous negative energy change");
        return true;
    }
    return random_number_between_zero_and_one <= std::exp(-energy_change);
}

/**
 * This performas the following tasks:
 * - syncs the two states
 * - syncs the two Hamiltonians
 * - resets the sum of energy changes to zero
 * - recalculates the initial energy
 */
void MetropolisMonteCarlo::init() {
    sum_of_energy_changes = 0.0;
    Change change;
    change.everything = true;

    state->pot->state = Energy::EnergyTerm::MonteCarloState::ACCEPTED;    // this is the old energy (current, accepted)
    trial_state->pot->state = Energy::EnergyTerm::MonteCarloState::TRIAL; // this is the new energy (trial)

    state->pot->init();
    auto energy = state->pot->energy(change);
    initial_energy = energy;
    faunus_logger->log(std::isfinite(initial_energy) ? spdlog::level::info : spdlog::level::warn,
                       "initial energy = {:.6E} kT", initial_energy);

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
    change.everything = true; // ... we want to calculate the total energy
    double energy = state->pot->energy(change);
    double du = energy - initial_energy;
    if (std::isfinite(du)) {
        if (std::fabs(du) <= pc::epsilon_dbl) {
            return 0.0;
        }
        return (energy - (initial_energy + sum_of_energy_changes)) /
               (std::fabs(initial_energy) > pc::epsilon_dbl ? initial_energy : energy);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

MetropolisMonteCarlo::MetropolisMonteCarlo(const json &j)
    : original_log_level(faunus_logger->level()) {
    state = std::make_unique<State>(j);
    faunus_logger->set_level(spdlog::level::off); // do not duplicate log info
    trial_state = std::make_unique<State>(j);     // ...for the trial state
    faunus_logger->set_level(original_log_level); // restore original log level
    moves = std::make_unique<move::MoveCollection>(j.at("moves"), *trial_state->spc, *trial_state->pot, *state->spc);
    init();
}

MetropolisMonteCarlo::~MetropolisMonteCarlo() = default;

/**
 * @todo Too many responsibilities; tidy up!
 */
void MetropolisMonteCarlo::restore(const json &j) {
    try {
        from_json(j, *state->spc);       // default, accepted state
        from_json(j, *trial_state->spc); // trial state
        if (j.contains("random-move")) {
            move::Move::slump = j["random-move"]; // restore move random number generator
        }
        if (j.contains("random-global")) {
            Faunus::random = j["random-global"]; // restore global random number generator
        }
        if (j.contains("reactionlist")) {
            faunus_logger->warn("'reactionlist' in state file is deprecated and will be ignored");
        }
        init();
    } catch (std::exception &e) {
        throw std::runtime_error("error initialising simulation: "s + e.what());
    }
}

void MetropolisMonteCarlo::performMove(move::Move& move) {
    Change change;
    move.move(change);
#ifndef NDEBUG
    try {
        change.sanityCheck(state->spc->groups);
    } catch (std::exception &e) {
        throw std::runtime_error(e.what());
    }
#endif
    if (change) {
        latest_move_name = move.getName();
        trial_state->pot->updateState(change);                    // update energy terms to reflect change
        const auto new_energy = trial_state->pot->energy(change); // trial potential energy (kT)
        const auto old_energy = state->pot->energy(change);       // potential energy before move (kT)

        auto energy_change = getEnergyChange(new_energy, old_energy);

        const auto energy_bias = move.bias(change, old_energy, new_energy) +
                                 TranslationalEntropy(*trial_state->spc, *state->spc).energy(change);

        const auto total_trial_energy = energy_change + energy_bias;
        if (std::isnan(total_trial_energy)) {
            faunus_logger->error("NaN energy change in {} move.", move.getName());
        }
        if (metropolisCriterion(total_trial_energy)) { // accept move
            state->sync(*trial_state, change);
            move.accept(change);
        } else { // reject move
            trial_state->sync(*state, change);
            move.reject(change);
            energy_change = 0.0;
        }
        sum_of_energy_changes += energy_change; // sum of all energy changes
        if (std::isfinite(initial_energy)) {
            average_energy += initial_energy + sum_of_energy_changes; // update average potential energy
        }
    } else {
        // The `metropolis()` function propagates the engine and we need to stay in sync
        // Alternatively, we could use `engine.discard()`
        move::Move::slump();
    }
}

/**
 * Policies for infinite/nan energy changes
 * @return modified energy change, new_energy - old_energy
 */
double MetropolisMonteCarlo::getEnergyChange(const double new_energy, const double old_energy) const {
    if (std::isnan(old_energy) and !std::isnan(new_energy)) { // if NaN --> finite energy change
        return pc::neg_infty;                                 // ...always accept
    }
    if (std::isnan(new_energy)) { // if moving to NaN, e.g. division by zero,
        return pc::infty;         // ...always reject
    }
    if (new_energy > 0.0 && std::isinf(new_energy)) { // if positive infinity
        return pc::infty;                             //...always reject
    }
    const auto energy_change = new_energy - old_energy; // potential energy change (kT)
    if (std::isnan(energy_change)) {                    // if difference is NaN, e.g. infinity - infinity,
        return 0.0;
    }
    return energy_change;
}

/**
 * This "sweeps" over all registered MC moves respecting, with the probability
 * of picking a move given by `Movebase::weight`. First stochastic moves
 * are randomly picked and, if needed, repeated (randomly).
 * Next, static moves are performed. Currently static moves
 * are defined by setting `weight=0`.
 *
 * @todo using `weight` to mark as move as static is ugly
 */
void MetropolisMonteCarlo::sweep() {
    assert(moves);
    number_of_sweeps++;
    auto perform_single_move = [&](auto& move) { performMove(*move); };
    ranges::cpp20::for_each(moves->repeatedStochasticMoves(), perform_single_move);
    ranges::cpp20::for_each(moves->constantIntervalMoves(number_of_sweeps), perform_single_move);
}

Energy::Hamiltonian &MetropolisMonteCarlo::getHamiltonian() { return *state->pot; }

Space &MetropolisMonteCarlo::getSpace() { return *state->spc; }

Space& MetropolisMonteCarlo::getTrialSpace() { return *trial_state->spc; }

void from_json(const json &j, MetropolisMonteCarlo::State &state) {
    state.spc = std::make_unique<Space>(j);
    state.pot = std::make_unique<Energy::Hamiltonian>(*state.spc, j.at("energy"));
}

void MetropolisMonteCarlo::State::sync(const State& other, const Change& change) {
    spc->sync(*other.spc, change);
    pot->sync(&*other.pot, change);
}

void to_json(json& j, const MetropolisMonteCarlo& monte_carlo) {
    j = monte_carlo.state->spc->info();
    j["temperature"] = pc::temperature / 1.0_K;
    if (monte_carlo.moves) {
        j["moves"] = *monte_carlo.moves;
    }
    j["number of sweeps"] = monte_carlo.number_of_sweeps;
    j["energy"].push_back(*monte_carlo.state->pot);
    if (!monte_carlo.average_energy.empty()) {
        j["montecarlo"] = {{"average potential energy (kT)", monte_carlo.average_energy.avg()},
                           {"last move", monte_carlo.latest_move_name}};
    }
}

TranslationalEntropy::TranslationalEntropy(const Space& trial_space, const Space& space)
    : trial_spc(trial_space)
    , spc(space) {}

/**
 * @param trial_count Number of atoms or molecules after move
 * @param count Number of atoms or molecular before move
 * @return Energy contribution (kT) to be added to MC trial energy
 */
double TranslationalEntropy::bias(int trial_count, int count) const {
    double energy = 0.0;
    if (int dN = trial_count - count; dN > 0) { // atoms or molecules were added
        double V_trial = trial_spc.geometry.getVolume();
        for (int n = 0; n < dN; n++) {
            energy += std::log((count + 1 + n) / (V_trial * 1.0_molar));
        }
    } else if (dN < 0) { // atoms or molecules were removed
        double V = spc.geometry.getVolume();
        for (int n = 0; n < (-dN); n++) {
            energy -= std::log((count - n) / (V * 1.0_molar));
        }
    }
    return energy; // kT
}

double TranslationalEntropy::atomSwapEnergy(const Change::GroupChange& group_change) const {
    assert(group_change.dNswap);
    assert(group_change.relative_atom_indices.size() == 1);
    double energy = 0.0;
    int id1 = trial_spc.groups.at(group_change.group_index).at(group_change.relative_atom_indices.front()).id;
    int id2 = spc.groups.at(group_change.group_index).at(group_change.relative_atom_indices.front()).id;
    for (auto atomid : {id1, id2}) {
        auto atoms_new = trial_spc.findAtoms(atomid);
        auto atoms_old = spc.findAtoms(atomid);
        int N_new = range_size(atoms_new); // number of atoms after change
        int N_old = range_size(atoms_old); // number of atoms before change
        energy += bias(N_new, N_old);
    }
    return energy; // kT
}

double TranslationalEntropy::atomChangeEnergy(const int molid) const {
    auto mollist_new = trial_spc.findMolecules(molid, Space::Selection::ALL); // "ALL" because "ACTIVE"
    auto mollist_old = spc.findMolecules(molid, Space::Selection::ALL);       // ...returns only full groups
    if (range_size(mollist_new) > 1 || range_size(mollist_old) > 1) {
        throw std::runtime_error("multiple atomic groups of the same type is not allowed");
    }
    int N_new = mollist_new.begin()->size(); // number of atoms after move
    int N_old = mollist_old.begin()->size(); // number of atoms before move
    return bias(N_new, N_old);
}

double TranslationalEntropy::moleculeChangeEnergy(const int molid) const {
    auto mollist_new = trial_spc.findMolecules(molid, Space::Selection::ACTIVE);
    auto mollist_old = spc.findMolecules(molid, Space::Selection::ACTIVE);
    int N_new = range_size(mollist_new); // number of molecules after move
    int N_old = range_size(mollist_old); // number of molecules before move
    return bias(N_new, N_old);
}

/**
 * @param change Change due to latest Monte Carlo move
 * @return Logarithm of the bias for the Metropolis criterion (units of kT)
 */
double TranslationalEntropy::energy(const Change& change) {
    double energy_change = 0.0;
    if (!change.matter_change || change.disable_translational_entropy) {
        return energy_change;
    }
    std::set<MoleculeData::index_type> already_processed;   // ignore future encounters of these molid's
    for (const Change::GroupChange& data : change.groups) { // loop over each change group
        if (data.dNswap) {                                  // number of atoms has changed as a result of a swap move
            energy_change += atomSwapEnergy(data);
        } else { // it is not a swap move
            const auto molid = trial_spc.groups.at(data.group_index).id;
            assert(molid == spc.groups.at(data.group_index).id);
            if (data.dNatomic and Faunus::molecules.at(molid).isAtomic()) { // an atomic group has been changed
                energy_change += atomChangeEnergy(molid);
            } else if (!already_processed.contains(molid)) { // a molecule has been inserted or deleted
                energy_change += moleculeChangeEnergy(molid);
                already_processed.insert(molid); // ignore future encounters of molid
            }
        }
    }
    return energy_change;
}
} // namespace Faunus
