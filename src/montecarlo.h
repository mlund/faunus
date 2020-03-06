#pragma once
#ifndef FAUNUS_MONTECARLO_H
#define FAUNUS_MONTECARLO_H

#include "energy.h"
#include "move.h"

namespace Faunus {
class MCSimulation {
  private:
    typedef typename Space::Tpvec Tpvec;

    std::string lastMoveName; //!< name of latest move

    spdlog::level::level_enum log_level; //!< Storage for original loglevel

    bool metropolis(double du) const; //!< Metropolis criterion (true=accept)

    struct State {
        Space spc;
        Energy::Hamiltonian pot;
        State(const json &j);

        void sync(State &other, Change &change);
    }; //!< Contains everything to describe a state

    State state1, // old state (accepted)
        state2;   // new state (trial)
    double uinit = 0, dusum = 0;
    Average<double> uavg;

    void init();

  public:
    Move::Propagator moves;

    auto &pot() { return state1.pot; }
    auto &space() { return state1.spc; }
    const auto &pot() const { return state1.pot; }
    const auto &space() const { return state1.spc; }
    const auto &geometry() const { return state1.spc.geo; }
    const auto &particles() const { return state1.spc.p; }

    MCSimulation(const json &, MPI::MPIController &);
    double drift(); //!< Calculates the relative energy drift from initial configuration

    /* currently unused -- see Analysis::SaveState.
                    void store(json &j) const {
                        j = state1.spc;
                        j["random-move"] = Move::Movebase::slump;
                        j["random-global"] = Faunus::random;
                    } // store system to json object
    */
    void restore(const json &j); //!< restore system from previously store json object
    void move();
    void to_json(json &j);
};

void to_json(json &j, MCSimulation &mc);

/**
 * @brief Ideal energy contribution of a speciation move
 *
 * Entropic, contribution to the change associated with a density fluctuation.
 * This far tested with particle fluctuations, only.
 *
 * @f[
 *     \beta \Delta U = - \sum \ln ( N_o!/N_n! V^{N_n - N_o} )
 * @f]
 *
 * where the sum runs over all products and reactants for the old (o) and
 * new (n) configuration.
 *
 * @todo
 * - [ ] Rename to something more meaningful
 * - [ ] Split into a function object with private functions for different change types
 * - [ ] `findMolecules` has O(N) complexity for non-atomic molecules; otherwise constant
 * - [ ] Move to energy.h?
 */

double IdealTerm(Space &, Space &, const Change &);
} // namespace Faunus
#endif // FAUNUS_MONTECARLO_H
