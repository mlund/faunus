#pragma once
#ifndef FAUNUS_MONTECARLO_H
#define FAUNUS_MONTECARLO_H

#include "space.h"
#include <memory>

namespace Faunus {

namespace Energy {
class Hamiltonian;
}

namespace move {
class Move;
class MoveCollection;
} // namespace move

/**
 * @brief Class to handle Monte Carlo moves
 *
 * This class maintains a list of MC moves and takes care of
 * randomly selecting a move; perform the move; decides if the move
 * is accepted; and checks if there is a drift in the potential energy.
 *
 * The moves always operate on a trial state (`new_state`). If accepted,
 * the changes will be copied into the `old_state`. The state here refers
 * to the particle states (positions etc.), the simulation geometry (size, volume),
 * and the state of the Hamiltonian (wave-vectors for Ewald etc.).
 *
 * @todo
 * The class has too many responsibilities, particularly in setting up the
 * system.
 */
class MetropolisMonteCarlo
{
  public:
    /**
     * @brief Class to describe a system state
     *
     * The state stores:
     * - per particle properties (positions, charges, etc.)
     * - per molecule properties (mass center, conformations)
     * - the simulation geometry (dimensions)
     * - state of the Hamiltonian (wave-vectors for Ewald etc.)
     *
     * A MC simulation consists of two states, the "default", old state
     * and the "trial" state. Moves operate directly on the latter.
     * After the move, the two states can be synchronized using the
     * `sync()` function.
     */
    struct State
    {
        std::unique_ptr<Space> spc; //!< Simulation space (positions, geometry, molecules)
        std::unique_ptr<Energy::Hamiltonian> pot; //!< Hamiltonian for calc. potential energy
        void
        sync(const State& other,
             const Change& change); //!< Sync with another state (the other state is not modified)
    };

  private:
    spdlog::level::level_enum original_log_level; //!< Storage for original loglevel
    std::unique_ptr<State> state;                 //!< The accepted MC state
    std::unique_ptr<State> trial_state;           //!< Proposed or trial MC state
    std::unique_ptr<move::MoveCollection> moves;  //!< Storage for all registered MC moves
    std::string latest_move_name;                 //!< Name of latest MC move
    double sum_of_energy_changes = 0.0;           //!< Sum of all potential energy changes
    double initial_energy = 0.0;                  //!< Initial potential energy
    Average<double> average_energy;               //!< Average potential energy of the system
    void init();                                  //!< Reset state
    void performMove(move::Move& move);           //!< Perform move using given move implementation
    double getEnergyChange(double new_energy, double old_energy) const;
    friend void to_json(json&, const MetropolisMonteCarlo&); //!< Write information to JSON object
    unsigned int number_of_sweeps = 0; //!< Number of MC sweeps, e.g. calls to sweep()

  public:
    MetropolisMonteCarlo(const json& j);
    Energy::Hamiltonian& getHamiltonian(); //!< Get Hamiltonian of accepted (default) state
    Space& getSpace();                     //!< Access to space in accepted (default) state
    Space& getTrialSpace();                //!< Access to trial space
    double relativeEnergyDrift();          //!< Relative energy drift from initial configuration
    void sweep();                          //!< Perform all moves (stochastic and static)
    void restore(const json& j);           //!< Restores system from previously store json object
    static bool metropolisCriterion(double energy_change); //!< Metropolis criterion
    ~MetropolisMonteCarlo(); //!< Required due to unique_ptr to incomplete type
};

void from_json(const json&, MetropolisMonteCarlo::State&); //!< Build state from json object
void to_json(json&, const MetropolisMonteCarlo&);

/**
 * @brief Entropy change due to particle fluctuations
 *
 * Entropic contribution to the change associated with a density fluctuation.
 * This far tested with particle fluctuations, only, i.e. when the number of
 * particles change due to a grand canonical or speciation move.
 *
 * @f[
 *     \beta \Delta U = - \sum \ln ( N_o!/N_n! V^{N_n - N_o} )
 * @f]
 *
 * The sum runs over all products and reactants for the old (o) and
 * new (n) configurations. The "energy" is not a real energy but merely the
 * logarithm of the bias to be included in the Metropolis criterion, so that it
 * can be *added* to the potential energy.
 *
 * @note `Space::findMolecules` has O(N) complexity for non-atomic molecules; otherwise constant
 * @todo
 * - [ ] Move to Energy namespace?
 * - [ ] Verify with volume fluctuations which would make `Energy::Isobaric` redundant
 */
class TranslationalEntropy
{
  private:
    const Space& trial_spc;                        //!< Space after proposed MC move ("trial")
    const Space& spc;                              //!< Space before MC move ("default")
    double bias(int trial_count, int count) const; //!< Bias due to change in atom/molecule numbers
    double atomSwapEnergy(
        const Change::GroupChange& group_change) const; //!< Contribution from atomic swap move
    double atomChangeEnergy(int molid) const; //!< Contribution from size-change of atomic group
    double moleculeChangeEnergy(
        int molid) const; //!< Contribution frin change in number of molecular groups

  public:
    TranslationalEntropy(const Space& trial_space, const Space& space);
    double energy(const Change& change); //!< Entropic contribution to MC trial energy
};

} // namespace Faunus
#endif // FAUNUS_MONTECARLO_H
