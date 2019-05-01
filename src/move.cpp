#include "core.h"
#include "move.h"
#include "speciation.h"

namespace Faunus {
namespace Move {

Random Movebase::slump; // static instance of Random (shared for all moves)

void Movebase::from_json(const json &j) {
    auto it = j.find("repeat");
    if (it != j.end()) {
        if (it->is_number())
            repeat = it->get<double>();
        else if (it->is_string())
            if (it->get<std::string>() == "N")
                repeat = -1;
    }
    _from_json(j);
    if (repeat < 0)
        repeat = 0;
}

void Movebase::to_json(json &j) const {
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

void Movebase::move(Change &change) {
    timer.start();
    timer_move.start();
    cnt++;
    change.clear();
    _move(change);
    if (change.empty())
        timer.stop();
    timer_move.stop();
}

void Movebase::accept(Change &c) {
    accepted++;
    _accept(c);
    timer.stop();
}

void Movebase::reject(Change &c) {
    rejected++;
    _reject(c);
    timer.stop();
}

double Movebase::bias(Change &, double, double) {
    return 0; // du
}

void Movebase::_accept(Change &) {}

void Movebase::_reject(Change &) {}

void from_json(const json &j, Movebase &m) { m.from_json(j); }

void to_json(json &j, const Movebase &m) {
    assert(!m.name.empty());
    m.to_json(j[m.name]);
}

void ChainRotationMovebase::_from_json(const json &j) {
    molname = j.at("molecule");
    dprot = j.at("dprot");
    allow_small_box = j.value("skiplarge", true); // todo rename the json attribute and make false default
}

void ChainRotationMovebase::_to_json(json &j) const {
    using namespace u8;
    j = {{"molecule", molname},
         {"dprot", dprot},
         {u8::rootof + u8::bracket("r_cm" + u8::squared), std::sqrt(msqdispl.avg())}};
    if (small_box_encountered > 0) {
        j["skipped"] = double(small_box_encountered) / cnt; // todo rename the json attribute
    }
    _roundjson(j, 3);
}

void ChainRotationMovebase::_move(Change &change) {
    permit_move = true;
    sqdispl = 0;
    if (std::fabs(dprot) > 1e-9) {
        if (select_segment() > 0) {
            double angle = dprot * (slump() - 0.5);
            rotate_segment(angle);
            store_change(change);
        }
    }
}

void ChainRotationMovebase::_accept(Change &) { msqdispl += sqdispl; }
void ChainRotationMovebase::_reject(Change &) { msqdispl += 0; }
double ChainRotationMovebase::bias(Change &, double, double) { return permit_move ? 0 : pc::infty; }

void AtomicTranslateRotate::_to_json(json &j) const {
    j = {{"dir", dir},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molname}};
    _roundjson(j, 3);
}
void AtomicTranslateRotate::_from_json(const json &j) {
    assert(!molecules.empty());
    try {
        assertKeys(j, {"molecule", "dir", "repeat"});
        molname = j.at("molecule");
        auto it = findName(molecules, molname);
        if (it == molecules.end())
            throw std::runtime_error("unknown molecule '" + molname + "'");
        molid = it->id();
        dir = j.value("dir", Point(1, 1, 1));
        if (repeat < 0) {
            auto v = spc.findMolecules(molid, Tspace::ALL);
            repeat = std::distance(v.begin(), v.end()); // repeat for each molecule...
            if (repeat > 0)
                repeat = repeat * v.front().size(); // ...and for each atom
        }
    } catch (std::exception &e) {
        std::cerr << name << ": " << e.what();
        throw;
    }
}
void AtomicTranslateRotate::translateParticle(
    std::vector<Faunus::ParticleTemplate<Faunus::Charge>,
                std::__1::allocator<Faunus::ParticleTemplate<Faunus::Charge>>>::iterator p,
    double dp) {
    auto &g = spc.groups[cdata.index];
    Point oldpos = p->pos;
    p->pos += ranunit(slump, dir) * dp * slump();

    spc.geo.boundary(p->pos);
    _sqd = spc.geo.sqdist(oldpos, p->pos); // squared displacement
    if (not g.atomic) {                    // recalc mass-center for non-molecular groups
        g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.cm);
#ifndef NDEBUG
        Point cmbak = g.cm;                             // backup mass center
        g.translate(-cmbak, spc.geo.getBoundaryFunc()); // translate to {0,0,0}
        double should_be_zero = spc.geo.sqdist({0, 0, 0}, Geometry::massCenter(g.begin(), g.end()));
        if (should_be_zero > 1e-6)
            throw std::runtime_error("atomic move too large");
        else
            g.translate(cmbak, spc.geo.getBoundaryFunc());
#endif
    }
}
void AtomicTranslateRotate::_move(Change &change) {
    auto p = randomAtom();
    if (p not_eq spc.p.end()) {
        double dp = atoms.at(p->id).dp;
        double dprot = atoms.at(p->id).dprot;

        if (dp > 0) // translate
            translateParticle(p, dp);

        if (dprot > 0) { // rotate
            Point u = ranunit(slump);
            double angle = dprot * (slump() - 0.5);
            Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, u));
            p->rotate(Q, Q.toRotationMatrix());
        }

        if (dp > 0 or dprot > 0)
            change.groups.push_back(cdata); // add to list of moved groups
    }
}
void AtomicTranslateRotate::_accept(Change &) { msqd += _sqd; }
void AtomicTranslateRotate::_reject(Change &) { msqd += 0; }
AtomicTranslateRotate::AtomicTranslateRotate(Tspace &spc) : spc(spc) {
    name = "transrot";
    repeat = -1; // meaning repeat N times
    cdata.atoms.resize(1);
    cdata.internal = true;
}
std::vector<Particle>::iterator AtomicTranslateRotate::randomAtom() {
    assert(molid >= 0);
    auto mollist = spc.findMolecules(molid, Tspace::ALL); // all `molid` groups
    if (size(mollist) > 0) {
        auto git = slump.sample(mollist.begin(), mollist.end()); // random molecule iterator
        if (not git->empty()) {
            auto p = slump.sample(git->begin(), git->end());         // random particle iterator
            cdata.index = Faunus::distance(spc.groups.begin(), git); // integer *index* of moved group
            cdata.atoms[0] = std::distance(git->begin(), p);         // index of particle rel. to group
            return p;
        }
    }
    return spc.p.end();
}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
Propagator::Propagator(const json &j, Tspace &spc, MPI::MPIController &mpi) {
#pragma GCC diagnostic pop

    if (j.count("random") == 1) {
        Movebase::slump = j["random"]; // slump is static --> shared for all moves
        Faunus::random = j["random"];
    }

    for (auto &m : j.at("moves")) { // loop over move list
        size_t oldsize = vec.size();
        for (auto it : m.items()) {
            try {
                if (it.key() == "moltransrot")
                    this->template push_back<Move::TranslateRotate<Tspace>>(spc);
                else if (it.key() == "conformationswap")
                    this->template push_back<Move::ConformationSwap<Tspace>>(spc);
                else if (it.key() == "transrot")
                    this->template push_back<Move::AtomicTranslateRotate>(spc);
                else if (it.key() == "pivot")
                    this->template push_back<Move::PivotMove<Tspace>>(spc);
                else if (it.key() == "crankshaft")
                    this->template push_back<Move::CrankshaftMove<Tspace>>(spc);
                else if (it.key() == "volume")
                    this->template push_back<Move::VolumeMove<Tspace>>(spc);
                else if (it.key() == "charge")
                    this->template push_back<Move::ChargeMove<Tspace>>(spc);
                else if (it.key() == "rcmc")
                    this->template push_back<Move::SpeciationMove>(spc);
                else if (it.key() == "quadrantjump")
                    this->template push_back<Move::QuadrantJump<Tspace>>(spc);
                else if (it.key() == "cluster")
                    this->template push_back<Move::Cluster<Tspace>>(spc);
                    // new moves go here...
#ifdef ENABLE_MPI
                else if (it.key() == "temper")
                    this->template push_back<Move::ParallelTempering>(spc, mpi);
                    // new moves requiring MPI go here...
#endif
                if (vec.size() == oldsize + 1) {
                    vec.back()->from_json(it.value());
                    addWeight(vec.back()->repeat);
                } else
                    throw std::runtime_error("unknown move");
            } catch (std::exception &e) {
                throw std::runtime_error("Error adding move '" + it.key() + "': " + e.what() + usageTip[it.key()]);
            }
        }
    }
}

void Propagator::addWeight(double weight) {
    w.push_back(weight);
    dist = std::discrete_distribution<>(w.begin(), w.end());
    _repeat = int(std::accumulate(w.begin(), w.end(), 0.0));
}

#ifdef ENABLE_MPI

void ParallelTempering::findPartner() {
    int dr = 0;
    partner = mpi.rank();
    (mpi.random() > 0.5) ? dr++ : dr--;
    (mpi.rank() % 2 == 0) ? partner += dr : partner -= dr;
}
bool ParallelTempering::goodPartner() {
    assert(partner != mpi.rank() && "Selfpartner! This is not supposed to happen.");
    if (partner >= 0)
        if (partner < mpi.nproc())
            if (partner != mpi.rank())
                return true;
    return false;
}
void ParallelTempering::_to_json(json &j) const {
    j = {{"replicas", mpi.nproc()}, {"datasize", pt.getFormat()}};
    json &_j = j["exchange"];
    _j = json::object();
    for (auto &m : accmap)
        _j[m.first] = {{"attempts", m.second.cnt}, {"acceptance", m.second.avg()}};
}
void ParallelTempering::_move(Change &change) {
    double Vold = spc.geo.getVolume();
    findPartner();
    Tpvec p; // temperary storage
    p.resize(spc.p.size());
    if (goodPartner()) {
        change.all = true;
        pt.sendExtra[VOLUME] = Vold;  // copy current volume for sending
        pt.recv(mpi, partner, p);     // receive particles
        pt.send(mpi, spc.p, partner); // send everything
        pt.waitrecv();
        pt.waitsend();

        double Vnew = pt.recvExtra[VOLUME];
        if (Vnew < 1e-9 || spc.p.size() != p.size())
            MPI_Abort(mpi.comm, 1);

        if (std::fabs(Vnew - Vold) > 1e-9)
            change.dV = true;

        spc.p = p;
        spc.geo.setVolume(Vnew);

        // update mass centers
        for (auto &g : spc.groups)
            if (g.atomic == false)
                g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.begin()->pos);
    }
}
double ParallelTempering::exchangeEnergy(double mydu) {
    std::vector<MPI::FloatTransmitter::floatp> duSelf(1), duPartner;
    duSelf[0] = mydu;
    duPartner = ft.swapf(mpi, duSelf, partner);
    return duPartner.at(0); // return partner energy change
}
double ParallelTempering::bias(Change &, double uold, double unew) {
    return exchangeEnergy(unew - uold); // Exchange dU with partner (MPI)
}
std::string ParallelTempering::id() {
    std::ostringstream o;
    if (mpi.rank() < partner)
        o << mpi.rank() << " <-> " << partner;
    else
        o << partner << " <-> " << mpi.rank();
    return o.str();
}
void ParallelTempering::_accept(Change &) {
    if (goodPartner())
        accmap[id()] += 1;
}
void ParallelTempering::_reject(Change &) {
    if (goodPartner())
        accmap[id()] += 0;
}
void ParallelTempering::_from_json(const json &j) { pt.setFormat(j.value("format", std::string("XYZQI"))); }
ParallelTempering::ParallelTempering(Tspace &spc, MPI::MPIController &mpi) : spc(spc), mpi(mpi) {
    name = "temper";
    partner = -1;
    pt.recvExtra.resize(1);
    pt.sendExtra.resize(1);
}
#endif

} // end of namespace Move

bool MCSimulation::metropolis(double du) const {
    if (std::isnan(du))
        throw std::runtime_error("Metropolis error: energy cannot be NaN");
    if (du < 0)
        return true;
    if (-du > pc::max_exp_argument)
        std::cerr << "warning: large metropolis energy" << std::endl;
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
    if (std::isfinite(u1) and std::isfinite(u2))
        if (std::fabs((u1 - u2) / u1) > 1e-3) {
            std::cerr << "u1 = " << u1 << "  u2 = " << u2 << endl;
            throw std::runtime_error("error aligning energies - this could be a bug...");
        }

    // inject reference to state1 in SpeciationMove (needed to calc. *differences*
    // in ideal excess chem. potentials)
    for (auto base : moves.vec) {
        auto derived = std::dynamic_pointer_cast<Move::SpeciationMove>(base);
        if (derived)
            derived->setOther(state1.spc);
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

                double bias = (**mv).bias(change, uold, unew) + IdealTerm(state2.spc, state1.spc, change);

                if (metropolis(du + bias)) { // accept move
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

MCSimulation::State::State(const json &j) : spc(j), pot(spc, j) {}

void MCSimulation::State::sync(MCSimulation::State &other, Change &change) {
    spc.sync(other.spc, change);
    pot.sync(&other.pot, change);
}

void to_json(json &j, MCSimulation &mc) { mc.to_json(j); }

double IdealTerm(Tspace &spc_n, Tspace &spc_o, const Change &change) {
    using Tpvec = typename Tspace::Tpvec;
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
                    auto mollist_n = spc_n.findMolecules(spc_n.groups[m.index].id, Tspace::ALL);
                    auto mollist_o = spc_o.findMolecules(spc_o.groups[m.index].id, Tspace::ALL);
                    if (size(mollist_n) > 1 || size(mollist_o) > 1)
                        throw std::runtime_error("Bad definition: One group per atomic molecule!");
                    if (not molecules[spc_n.groups[m.index].id].atomic)
                        throw std::runtime_error("Only atomic molecules!");
                    // Below is safe due to the catches above
                    // add consistency criteria with m.atoms.size() == N
                    N_n = mollist_n.begin()->size();
                    N_o = mollist_o.begin()->size();
                } else {
                    auto mollist_n = spc_n.findMolecules(spc_n.groups[m.index].id, Tspace::ACTIVE);
                    auto mollist_o = spc_o.findMolecules(spc_o.groups[m.index].id, Tspace::ACTIVE);
                    N_n = size(mollist_n);
                    N_o = size(mollist_o);
                }
                int dN = N_n - N_o;
                if (dN != 0) {
                    double V_n = spc_n.geo.getVolume();
                    double V_o = spc_o.geo.getVolume();
                    if (dN > 0)
                        for (int n = 0; n < dN; n++)
                            NoverO += std::log((N_o + 1 + n) / (V_n * 1.0_molar));
                    else
                        for (int n = 0; n < (-dN); n++)
                            NoverO -= std::log((N_o - n) / (V_o * 1.0_molar));
                }
            }
        }
    }
    return NoverO;
}
} // namespace Faunus
