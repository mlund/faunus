#include "montecarlo.h"
#include "speciation.h"
#include "spdlog/spdlog.h"

namespace Faunus {

bool MCSimulation::metropolis(double du) const {
    if (std::isnan(du))
        throw std::runtime_error("Metropolis error: energy cannot be NaN");
    if (du < 0)
        return true;
    if (-du > pc::max_exp_argument)
        mcloop_logger->warn("warning: large metropolis energy");
    return (Move::Movebase::slump() > std::exp(-du)) ? false : true;
}

void MCSimulation::init() {
    dusum = 0;
    Change c;
    c.all = true;

    state1.pot.key = Energy::Energybase::OLD; // this is the old energy (current, accepted)
    state2.pot.key = Energy::Energybase::NEW; // this is the new energy (trial)

    state1.pot.init();
    double u1 = state1.pot.energy(c);
    uinit = u1;

    state2.sync(state1, c); // copy all information from state1 into state2
    state2.pot.init();
    double u2 = state2.pot.energy(c);

    // check that the energies in state1 and state2 are *identical*
    if (std::isfinite(u1) and std::isfinite(u2)) {
        if (std::fabs((u1 - u2) / u1) > 1e-3) {
            std::cerr << "u1 = " << u1 << "  u2 = " << u2 << std::endl;
            throw std::runtime_error("error aligning energies - this could be a bug...");
        }
    }

    // inject reference to state1 in SpeciationMove (needed to calc. *differences*
    // in ideal excess chem. potentials)
    for (auto speciation_move : moves.moves().find<Move::SpeciationMove>()) {
        speciation_move->setOther(state1.spc);
    }
}

double MCSimulation::drift() {
    Change c;
    c.all = true;
    double ufinal = state1.pot.energy(c);
    double du = ufinal - uinit;
    if (std::isfinite(du)) {
        if (std::fabs(du) < 1e-10)
            return 0;
        if (uinit != 0)
            return (ufinal - (uinit + dusum)) / uinit;
        else if (ufinal != 0)
            return (ufinal - (uinit + dusum)) / ufinal;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

/*
 * We need to construct two identical State objects and to avoid duplicate logs, we
 * temporarily disable the logger for the second object by the arcane _comma operator_
 */
MCSimulation::MCSimulation(const json &j, MPI::MPIController &mpi)
    : log_level(faunus_logger->level()), state1(j), state2((faunus_logger->set_level(spdlog::level::off), j)),
      moves((faunus_logger->set_level(log_level), j), state2.spc, mpi) {
    init();
}

void MCSimulation::restore(const json &j) {
    try {
        state1.spc = j; // old/accepted state
        state2.spc = j; // trial state
        if (j.count("random-move") == 1)
            Move::Movebase::slump = j["random-move"]; // restore move random number generator
        if (j.count("random-global") == 1)
            Faunus::random = j["random-global"];                     // restore global random number generator
        reactions = j.at("reactionlist").get<decltype(reactions)>(); // should be handled by space
        init();
    } catch (std::exception &e) {
        throw std::runtime_error("error initialising simulation: "s + e.what());
    }
}

void MCSimulation::move() {
    Change change;
    for (int i = 0; i < moves.repeat(); i++) {
        auto mv = moves.sample(); // pick random move
        if (mv != moves.end()) {
            change.clear();
            (**mv).move(change);
#ifndef NDEBUG
            // check if atom index indeed belong to the group (index)
            if (not change.sanityCheck(state1.spc))
                throw std::runtime_error("insane change object\n" + json(change).dump(4));
#endif
            if (change) {
                lastMoveName = (**mv).name; // store name of move for output
                double unew, uold, du;
                //#pragma omp parallel sections
                {
                    //#pragma omp section
                    { unew = state2.pot.energy(change); }
                    //#pragma omp section
                    { uold = state1.pot.energy(change); }
                }

                du = unew - uold;

                // if any energy returns NaN (from i.e. division by zero), the
                // configuration will always be rejected, or if moving from NaN
                // to a finite energy, always accepted.

                if (std::isnan(uold) and not std::isnan(unew))
                    du = -pc::infty; // accept
                else if (std::isnan(unew))
                    du = pc::infty; // reject

                // if the difference in energy is NaN (from i.e. infinity minus infinity), the
                // configuration will always be accepted. This should be
                // noted during equilibration.

                else if (std::isnan(du))
                    du = 0; // accept

                double bias = (**mv).bias(change, uold, unew);
                double ideal = IdealTerm(state2.spc, state1.spc, change);
                if (std::isnan(du + bias))
                    faunus_logger->error("Infinite du + bias in " + lastMoveName + " move.");

                if (metropolis(du + bias + ideal)) { // accept move
                    state1.sync(state2, change);
                    (**mv).accept(change);
                } else { // reject move
                    state2.sync(state1, change);
                    (**mv).reject(change);
                    du = 0;
                }
                dusum += du; // sum of all energy changes
            }
        }
    }
}

void MCSimulation::to_json(json &j) {
    j = state1.spc.info();
    j["temperature"] = pc::temperature / 1.0_K;
    j["moves"] = moves;
    j["energy"].push_back(state1.pot);
    j["last move"] = lastMoveName;
}

MCSimulation::State::State(const json &j) : spc(j), pot(spc, j.at("energy")) {}

void MCSimulation::State::sync(MCSimulation::State &other, Change &change) {
    spc.sync(other.spc, change);
    pot.sync(&other.pot, change);
}

void to_json(json &j, MCSimulation &mc) { mc.to_json(j); }

double IdealTerm(Space &spc_new, Space &spc_old, const Change &change) {
    if (not change.dN) {
        return 0.0;
    }
    double NoverO = 0;
    for (const Change::data &m : change.groups) { // loop over each change group
        assert(not change.empty());
        int N_old = 0;                            // number of molecules/atoms before change
        int N_new = 0;                            // number of molecules/atoms after change
        if (m.dNswap) {                           // the number of atoms has changed as a result of a swap move
            assert(m.atoms.size() == 1);
            auto &new_group = spc_new.groups.at(m.index);
            auto &old_group = spc_old.groups.at(m.index);
            int id1 = (new_group.begin() + m.atoms.front())->id;
            int id2 = (old_group.begin() + m.atoms.front())->id;
            for (int atom_id : {id1, id2}) {
                auto mollist_new = spc_new.findAtoms(atom_id);
                auto mollist_old = spc_old.findAtoms(atom_id);
                N_new = range_size(mollist_new);
                N_old = range_size(mollist_old);
                int dN = N_new - N_old;
                double V_new = spc_new.geo.getVolume();
                double V_old = spc_old.geo.getVolume();
                if (dN > 0) {
                    for (int n = 0; n < dN; n++) {
                        NoverO += std::log((N_old + 1 + n) / (V_new * 1.0_molar));
                    }
                } else {
                    for (int n = 0; n < (-dN); n++)
                        NoverO -= std::log((N_old - n) / (V_old * 1.0_molar));
                }
            }
        } else {              // it is not a swap move
            int molid = spc_new.groups.at(m.index).id;
            assert(molid == spc_old.groups.at(m.index).id);
            if (m.dNatomic) { // the number of atomic molecules has changed
                auto mollist_n = spc_new.findMolecules(molid, Space::ALL); // why not "ACTIVE"?
                auto mollist_o = spc_old.findMolecules(molid, Space::ALL);

                if (range_size(mollist_n) > 1 || range_size(mollist_o) > 1)
                    throw std::runtime_error("Bad definition: One group per atomic molecule!");
                if (not molecules[molid].atomic)
                    throw std::runtime_error("Only atomic molecules!");

                // Below is safe due to the catches above
                // add consistency criteria with m.atoms.size() == N
                N_new = mollist_n.begin()->size();
                N_old = mollist_o.begin()->size();
            } else {
                assert(spc_new.groups[m.index].empty() != spc_old.groups[m.index].empty());
                auto mollist_n = spc_new.findMolecules(molid, Space::ACTIVE);
                auto mollist_o = spc_old.findMolecules(molid, Space::ACTIVE);
                N_new = range_size(mollist_n);
                N_old = range_size(mollist_o);
            }
            int dN = N_new - N_old;
            if (dN != 0) {
                double V = spc_new.geo.getVolume() * 1.0_molar;
                if (dN > 0) {
                    for (int n = 0; n < dN; n++) {
                        NoverO += std::log((N_old + 1 + n) / V);
                    }
                } else {
                    for (int n = 0; n < (-dN); n++) {
                        NoverO -= std::log((N_old - n) / V);
                    }
                }
            }
        }
    }
    return NoverO;
}

} // namespace Faunus
