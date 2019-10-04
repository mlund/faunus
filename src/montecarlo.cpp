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

MCSimulation::MCSimulation(const json &j, MPI::MPIController &mpi) : state1(j), state2(j), moves(j, state2.spc, mpi) {
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

                if (std::isnan(du + bias + ideal)) {
                    std::cout << "NaN energy for " << lastMoveName << " move." << std::endl;
                    std::cout << "du: " << du << std::endl;
                    std::cout << "Bias: " << bias << std::endl;
                    std::cout << "Ideal: " << ideal << std::endl;
                }

                if (metropolis(du + bias + ideal)) { // accept move
                    state1.sync(state2, change);
                    (**mv).accept(change);
                } else { // reject move
                    state2.sync(state1, change);
                    (**mv).reject(change);
                    du = 0;
                }
                dusum += du; // sum of all energy changes
            } else state2.sync(state1, change);
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

double IdealTerm(Space &spc_n, Space &spc_o, const Change &change) {
    using Tpvec = typename Space::Tpvec;
    double NoverO = 0;
    if (change.dN) { // Has the number of any molecules changed?
        for (auto &m : change.groups) {
            int N_o = 0;
            int N_n = 0;
            if (m.dNswap) {
                assert(m.atoms.size() == 1);
                auto &g_n = spc_n.groups.at(m.index);
                int id1 = (g_n.begin() + m.atoms.front())->id;
                auto &g_o = spc_o.groups.at(m.index);
                int id2 = (g_o.begin() + m.atoms.front())->id;
                for (int id : {id1, id2}) {
                    auto atomlist_n = spc_n.findAtoms(id);
                    auto atomlist_o = spc_o.findAtoms(id);
                    N_n = size(atomlist_n);
                    N_o = size(atomlist_o);
                    int dN = N_n - N_o;
                    double V_n = spc_n.geo.getVolume();
                    double V_o = spc_o.geo.getVolume();
                    if (dN > 0)
                        for (int n = 0; n < dN; n++)
                            NoverO += std::log((N_o + 1 + n) / (V_n * 1.0_molar));
                    else
                        for (int n = 0; n < (-dN); n++)
                            NoverO -= std::log((N_o - n) / (V_o * 1.0_molar));
                }
            } else {
                if (m.dNatomic) {
                    auto mollist_n = spc_n.findMolecules(spc_n.groups[m.index].id, Space::ALL);
                    auto mollist_o = spc_o.findMolecules(spc_o.groups[m.index].id, Space::ALL);
                    if (size(mollist_n) > 1 || size(mollist_o) > 1)
                        throw std::runtime_error("Bad definition: One group per atomic molecule!");
                    if (not molecules[spc_n.groups[m.index].id].atomic)
                        throw std::runtime_error("Only atomic molecules!");
                    // Below is safe due to the catches above
                    // add consistency criteria with m.atoms.size() == N
                    N_n = mollist_n.begin()->size();
                    N_o = mollist_o.begin()->size();
                } else {
                    auto mollist_n = spc_n.findMolecules(spc_n.groups[m.index].id, Space::ACTIVE);
                    auto mollist_o = spc_o.findMolecules(spc_o.groups[m.index].id, Space::ACTIVE);
                    N_n = size(mollist_n);
                    N_o = size(mollist_o);
                }
                int dN = N_n - N_o;
                if (dN != 0) {
                    double V = spc_n.geo.getVolume() * 1.0_molar;
                    if (dN > 0)
                        for (int n = 0; n < dN; n++)
                            NoverO += std::log( (N_o + 1 + n) / V );
                    else
                        for (int n = 0; n < (-dN); n++)
                            NoverO -= std::log( (N_o - n) / V );
                }
            }
        }
    }
    return NoverO;
}

} // namespace Faunus
