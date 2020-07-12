#pragma once
#ifndef FAUNUS_MONTECARLO_H
#define FAUNUS_MONTECARLO_H

#include "energy.h"
#include "move.h"

namespace Faunus {

/**
 * @brief Class to handle Monte Carlo moves
 *
 * This class maintains a list of MC moves and takes care of
 * randomly selecting a moves, perform the move, decide if the move
 * if accepted, and also checks if there is a drift in the potential
 * energy.
 *
 * The moves always operate on a trial state (`new_state`). If accepted,
 * the changes will be copied into the `old_state`. The state here refers
 * to the particle states (positions etc.), the simulation geometry (size, volume),
 * and the state of the Hamiltonian (wave-vectors for Ewald etc.).
 *
 * @todo
 * The class has too many responsibilities, particularly in setting up the
 * system. It for example takes care of setting the global temperature(!).
 */
class MetropolisMonteCarlo {
  private:
    struct State {
        Space spc;
        Energy::Hamiltonian pot;
        State(const json &);
        void sync(State &other, Change &change);
    }; //!< Class to describe a "state" incl. Space and Hamiltonian

    spdlog::level::level_enum original_log_level; //!< Storage for original loglevel
    State old_state;                              //!< Old state representing the accepted MC state
    State new_state;                              //!< New state representing the MC trial state
    std::shared_ptr<Move::Movebase> latest_move;  //!< Pointer to latest MC move
    double initial_energy = 0;                    //!< Initial potential energy
    double sum_of_energy_changes = 0;             //!< Sum of all potential energy changes
    Average<double> average_energy;               //!< Average system energy
    void init();                                  //!< Reset state
    Move::Propagator moves;                       //!< Storage for all registered MC moves
    bool metropolis(double du) const;             //!< Metropolis criterion

  public:
    MetropolisMonteCarlo(const json &, MPI::MPIController &);
    Energy::Hamiltonian &getHamiltonian();             //!< Get Hamiltonian of accepted (old) state
    const Energy::Hamiltonian &getHamiltonian() const; //!< Get Hamiltonian of accepted (old) state
    Space &getSpace();                                 //!< Access to space in accepted (old) state
    const Space &getSpace() const;                     //!< Access to space in accepted (old) state
    double relativeEnergyDrift();                      //!< Relative energy drift from initial configuration
    void move();                                       //!< Perform random Monte Carlo move
    void restore(const json &);                        //!< restore system from previously store json object
    void to_json(json &);                              //!< Write information to JSON object
};

void to_json(json &, MetropolisMonteCarlo &);

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
 * new (n) configurations. The "energy" is not really an energy but merely the
 * logarithm of the bias to be included in the Metropolis criterion, i.e. added
 * to the potential energy.
 *
 * @note `Space::findMolecules` has O(N) complexity for non-atomic molecules; otherwise constant
 * @todo
 * - [ ] Move to Energy namespace?
 * - [ ] Verify with volume fluctuations which would make `Energy::Isobaric` redundant
 */
class TranslationalEntropy {
  private:
    Space &spc_new;                              //!< Space after proposed MC move ("new" or "trial")
    Space &spc_old;                              //!< Space before MC move ("old")
    double bias(int N_new, int N_old) const;     //!< Bias due to change in atom/molecule numbers
    double atomSwapEnergy(const Change::data &); //!< Contribution from atomic swap move
    double atomChangeEnergy(int molid);          //!< Contribution from size-change of atomic group
    double moleculeChangeEnergy(int molid);      //!< Contribution frin change in number of molecular groups

  public:
    TranslationalEntropy(Space &new_space, Space &old_space);
    double energy(const Change &); //!< Entropic contribution to MC trial energy
};

} // namespace Faunus
#endif // FAUNUS_MONTECARLO_H
