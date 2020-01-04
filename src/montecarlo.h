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
 * This funciton calculates the contribution to the energy change arising from the
 * change in concentration of reactant and products in the current and in the trial state.
 *
 * @f[
 *     \beta \Delta U = - \sum \ln ( N_o!/N_n! V^{N_n - N_o} )
 * @f]
 *
 * where the sum runs over all products and reactants.
 *
 * @todo
 * - use exception message to suggest how to fix the problem
 * - seems out of place; move to another file?
 * - find a better name
 */
double IdealTerm(Space &spc_n, Space &spc_o, const Change &change);
} // end of Faunus namespace
#endif // FAUNUS_MONTECARLO_H
