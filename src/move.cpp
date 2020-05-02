#include "core.h"
#include "move.h"
#include "speciation.h"
#include "clustermove.h"
#include "chainmove.h"
#include "aux/iteratorsupport.h"
#include "aux/eigensupport.h"
#include "spdlog/spdlog.h"

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
        if (Faunus::molecules[molid].rigid) {
            faunus_logger->warn("structure of rigid molecule {} may be disturbed by {}", molname, name);
        }
        dir = j.value("dir", Point(1, 1, 1));
        if (repeat < 0) {
            auto v = spc.findMolecules(molid, Space::ALL);
            repeat = std::distance(v.begin(), v.end()); // repeat for each molecule...
            if (repeat > 0)
                repeat = repeat * v.front().size(); // ...and for each atom
        }
    } catch (std::exception &e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}
void AtomicTranslateRotate::translateParticle(Space::Tpvec::iterator p, double dp) {
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
AtomicTranslateRotate::AtomicTranslateRotate(Space &spc) : spc(spc) {
    name = "transrot";
    repeat = -1; // meaning repeat N times
    cdata.atoms.resize(1);
    cdata.internal = true;
}
std::vector<Particle>::iterator AtomicTranslateRotate::randomAtom() {
    assert(molid >= 0);
    auto mollist = spc.findMolecules(molid, Space::ALL); // all `molid` groups
    if (not ranges::cpp20::empty(mollist)) {
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
Propagator::Propagator(const json &j, Space &spc, MPI::MPIController &mpi) {
#pragma GCC diagnostic pop

    if (j.count("random") == 1) {
        Movebase::slump = j["random"]; // slump is static --> shared for all moves
        Faunus::random = j["random"];
    }

    for (auto &m : j.at("moves")) { // loop over move list
        size_t oldsize = _moves.size();
        for (auto it : m.items()) {
            try {
                if (it.key() == "moltransrot")
                    _moves.emplace_back<Move::TranslateRotate>(spc);
                else if (it.key() == "smartmoltransrot")
                    _moves.emplace_back<Move::SmartTranslateRotate>(spc);
                else if (it.key() == "conformationswap")
                    _moves.emplace_back<Move::ConformationSwap>(spc);
                else if (it.key() == "transrot")
                    _moves.emplace_back<Move::AtomicTranslateRotate>(spc);
                else if (it.key() == "pivot")
                    _moves.emplace_back<Move::PivotMove>(spc);
                else if (it.key() == "crankshaft")
                    _moves.emplace_back<Move::CrankshaftMove>(spc);
                else if (it.key() == "volume")
                    _moves.emplace_back<Move::VolumeMove>(spc);
                else if (it.key() == "charge")
                    _moves.emplace_back<Move::ChargeMove>(spc);
                else if (it.key() == "chargetransfer")
                    _moves.emplace_back<Move::ChargeTransfer>(spc);
                else if (it.key() == "rcmc")
                    _moves.emplace_back<Move::SpeciationMove>(spc);
                else if (it.key() == "quadrantjump")
                    _moves.emplace_back<Move::QuadrantJump>(spc);
                else if (it.key() == "cluster")
                    _moves.emplace_back<Move::Cluster>(spc);
                    // new moves go here...
#ifdef ENABLE_MPI
                else if (it.key() == "temper")
                    _moves.emplace_back<Move::ParallelTempering>(spc, mpi);
                    // new moves requiring MPI go here...
#endif
                if (_moves.size() == oldsize + 1) {
                    _moves.back()->from_json(it.value());
                    addWeight(_moves.back()->repeat);
                } else
                    throw std::runtime_error("unknown move");
            } catch (std::exception &e) {
                throw std::runtime_error("Error adding move '" + it.key() + "': " + e.what() + usageTip[it.key()]);
            }
        }
    }
}

void Propagator::addWeight(double weight) {
    _weights.push_back(weight);
    distribution = std::discrete_distribution<>(_weights.begin(), _weights.end());
    _repeat = int(std::accumulate(_weights.begin(), _weights.end(), 0.0));
}

void to_json(json &j, const Propagator &propagator) {
    j = propagator._moves;
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
        // store group sizes
        for (auto &g : spc.groups) {
    	    pt.sendExtra.push_back((float)g.size());
        }
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

        size_t i = 0;
        for (auto &g : spc.groups) {
            // assign correct sizes to the groups
            g.resize((int)pt.recvExtra[i+1]);
            if (g.atomic == false) {
                // update mass center of molecular groups
                g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.begin()->pos);
            }
            ++i;
        }
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
ParallelTempering::ParallelTempering(Space &spc, MPI::MPIController &mpi) : spc(spc), mpi(mpi) {
    name = "temper";
    partner = -1;
    pt.recvExtra.resize(1);
    pt.sendExtra.resize(1);
}
#endif

void VolumeMove::_to_json(json &j) const {
    using namespace u8;
    if (cnt > 0) {
        j = {{"dV", dV},
             {"method", method->first},
             {bracket("V"), Vavg.avg()},
             {rootof + bracket(Delta + "V" + squared), std::sqrt(msqd.avg())},
             {cuberoot + rootof + bracket(Delta + "V" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
        _roundjson(j, 3);
    }
}
void VolumeMove::_from_json(const json &j) {
    try {
        method = methods.find(j.value("method", "isotropic"));
        if (method == methods.end())
            std::runtime_error("unknown volume change method");
        dV = j.at("dV");
    } catch (std::exception &e) {
        throw std::runtime_error(e.what());
    }
}
void VolumeMove::_move(Change &change) {
    if (dV > 0) {
        change.dV = true;
        change.all = true;
        Vold = spc.geo.getVolume();
        Vnew = std::exp(std::log(Vold) + (slump() - 0.5) * dV);
        deltaV = Vnew - Vold;
        spc.scaleVolume(Vnew, method->second);
    } else
        deltaV = 0;
}
void VolumeMove::_accept(Change &) {
    msqd += deltaV * deltaV;
    Vavg += spc.geo.getVolume();
}
VolumeMove::VolumeMove(Space &spc) : spc(spc) {
    name = "volume";
    repeat = 1;
}
void VolumeMove::_reject(Change &) {
    msqd += 0;
    Vavg += spc.geo.getVolume();
}

void ChargeMove::_to_json(json &j) const {
    using namespace u8;
    j = {{"index", atomIndex},
         {"dq", dq},
         {rootof + bracket(Delta + "q" + squared), std::sqrt(msqd.avg())},
         {cuberoot + rootof + bracket(Delta + "q" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
    _roundjson(j, 3);
}
void ChargeMove::_from_json(const json &j) {
    dq = j.at("dq").get<double>();
    atomIndex = j.at("index").get<int>();
    auto git = spc.findGroupContaining(spc.p[atomIndex]);                    // group containing atomIndex
    cdata.index = std::distance(spc.groups.begin(), git);                    // integer *index* of moved group
    cdata.atoms[0] = std::distance(git->begin(), spc.p.begin() + atomIndex); // index of particle rel. to group
}
void ChargeMove::_move(Change &change) {
    if (dq > 0) {
        auto &p = spc.p[atomIndex]; // refence to particle
        double qold = p.charge;
        p.charge += dq * (slump() - 0.5);
        deltaq = p.charge - qold;
        change.groups.push_back(cdata); // add to list of moved groups
    } else
        deltaq = 0;
}
void ChargeMove::_accept(Change &) { msqd += deltaq * deltaq; }
void ChargeMove::_reject(Change &) { msqd += 0; }
ChargeMove::ChargeMove(Space &spc) : spc(spc) {
    name = "charge";
    repeat = 1;
    cdata.internal = true; // the group is internally changed
    cdata.atoms.resize(1); // we change exactly one atom
}

void ChargeTransfer::_to_json(json &j) const {
    using namespace u8;
    j = {{"dq", dq},
         {rootof + bracket(Delta + "q" + squared), std::sqrt(msqd.avg())},
         {cuberoot + rootof + bracket(Delta + "q" + squared), std::cbrt(std::sqrt(msqd.avg()))}};
    _roundjson(j, 3);
}
void ChargeTransfer::_from_json(const json &j) {
    try {
        dq = j.at("dq").get<double>();
        mol1.molname = j.at("mol1"); // string containing name of molecule 1
        mol2.molname = j.at("mol2"); // string containing name of molecule 2
        mol1.molrange =
            j.at("molrange1")
                .get<std::vector<double>>(); // vector containing lower and upper limit of total charge of molecule 1
        mol2.molrange =
            j.at("molrange2")
                .get<std::vector<double>>(); // vector containing lower and upper limit of total charge of molecule 2
        mol1.min =
            j.at("min1").get<std::vector<double>>(); // vector containing lower limits of atomic charges in molecule 1
        mol1.max =
            j.at("max1").get<std::vector<double>>(); // vector containing upper limits of atomic charges in molecule 1
        mol2.min =
            j.at("min2").get<std::vector<double>>(); // vector containing lower limits of atomic charges in molecule 2
        mol2.max =
            j.at("max2").get<std::vector<double>>();   // vector containing upper limits of atomic charges in molecule 2
        auto git1 = findName(molecules, mol1.molname); // group containing mol1.molname
        auto git2 = findName(molecules, mol2.molname); // group containing mol2.molname

        if (git1 == molecules.end()) // checking so that molecule1 exists
            throw std::runtime_error("unknown molecule '" + mol1.molname + "'");

        mol1.id = git1->id();
        mol2.id = git2->id();

        if (repeat < 0) {
            auto v = spc.findMolecules(mol1.id);
            repeat = std::distance(v.begin(), v.end());
        }

        if (git2 == molecules.end()) // checking so that molecule2 exists
            throw std::runtime_error("unknown molecule '" + mol2.molname + "'");

        if (repeat < 0) {
            auto v = spc.findMolecules(mol2.id);
            repeat = std::distance(v.begin(), v.end());
        }

        if (mol1.min.size() !=
            mol1.max.size()) // checking so that mol1.min and mol1.max contains equal number of entries
            throw std::runtime_error("mol1.min and mol1.max need to have the same number of entries. mol1.min has " +
                                     std::to_string(mol1.min.size()) + " and mol1.max has " +
                                     std::to_string(mol1.max.size()) + " entries");

        if (mol1.min.size() == 0 || mol1.max.size() == 0) // checking so that mol1.min and mol1.max are not empty
            throw std::runtime_error(
                "mol1.min and mol1.max both need to have nonzero number of entries. mol1.min has " +
                std::to_string(mol1.min.size()) + " and mol1.max has " + std::to_string(mol1.max.size()) + " entries");

        if (mol2.min.size() !=
            mol2.max.size()) // checking so that mol2.min and mol2.max contains equal number of entries
            throw std::runtime_error("mol2.min and mol2.max need to have the same number of entries. mol2.min has " +
                                     std::to_string(mol2.min.size()) + " and mol2.max has " +
                                     std::to_string(mol2.max.size()) + " entries");

        if (mol2.min.size() == 0 || mol2.max.size() == 0) // checking so that mol2.min and mol2.max are not empty
            throw std::runtime_error(
                "mol2.min and mol2.max both need to have nonzero number of entries. mol2.min has " +
                std::to_string(mol2.min.size()) + " and mol2.max has " + std::to_string(mol2.max.size()) + " entries");

    } catch (std::exception &e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}

void ChargeTransfer::_move(Change &change) {

    auto mollist1 = spc.findMolecules(mol1.id, Space::ACTIVE);
    auto mollist2 = spc.findMolecules(mol2.id, Space::ACTIVE);
    if ((not ranges::cpp20::empty(mollist1)) and (not ranges::cpp20::empty(mollist2))) {
        auto git1 = slump.sample(mollist1.begin(), mollist1.end()); // selecting a random molecule of type molecule1
        auto git2 = slump.sample(mollist2.begin(), mollist2.end()); // selecting a random molecule of type molecule2

        if (!git1->empty() && !git2->empty()) { // check that both molecule1 and molecule 2 exist

            if (dq > 0) {
                // change.chargeMove =
                //    true; // setting to true makes the self-energy being computed and added to the total energy
                mol1.numOfAtoms = Faunus::distance(git1->begin(), git1->end());
                mol2.numOfAtoms = Faunus::distance(git2->begin(), git2->end());

                mol1.ratio.clear(); // clearing vector containing ratio of atomic charge ranges and the charge range of
                                    // the whole molecule1
                mol2.ratio.clear(); // clearing vector containing ratio of atomic charge ranges and the charge range of
                                    // the whole molecule2

                for (i = 0; i < mol1.numOfAtoms; i++) {
                    mol1.ratio.push_back(
                        (mol1.max[i] - mol1.min[i]) /
                        (mol1.molrange[1] - mol1.molrange[0])); // calculating ratio of atom i in molecule 1
                }

                for (i = 0; i < mol2.numOfAtoms; i++) {
                    mol2.ratio.push_back(
                        (mol2.max[i] - mol2.min[i]) /
                        (mol2.molrange[1] - mol2.molrange[0])); // calculating ratio of atom i in molecule 2
                }

                mol1.charges = 0; // setting sum of all atomic charges in molecule1 to zero
                mol2.charges = 0; // setting sum of all atomic charges in molecule2 to zero
                deltaq = dq * (slump() - 0.5);
                mol1.changeQ.clear(); // clearing vector containing attempted charge moves on all atoms in molecule1
                mol2.changeQ.clear(); // clearing vector containing attempted charge moves on all atoms in molecule2
                mol1.cdata.index = Faunus::distance(spc.groups.begin(), git1);
                mol2.cdata.index = Faunus::distance(spc.groups.begin(), git2);

                for (i = 0; i < mol1.numOfAtoms; i++) {
                    auto p = git1->begin() + i; // object containing atom i in molecule1
                    mol1.changeQ.push_back(
                        deltaq * mol1.ratio[i]); // assigning attempted charge move of atom i in molecule1 to vector
                    // sumChanges1 += changeQ1[i];
                    mol1.charges +=
                        p->charge + mol1.changeQ[i]; // adding new attempted charge of atom i in molecule1 to sum
                }
                for (i = 0; i < mol2.numOfAtoms; i++) { // Doing the same as above loop but for molecule2
                    auto p = git2->begin() + i;
                    mol2.changeQ.push_back(-deltaq * mol2.ratio[i]);
                    // sumMoves2 += changeQ2[i];
                    mol2.charges += p->charge + mol2.changeQ[i];
                }

                // Torodial boundary conditions
                if (mol1.charges < mol1.molrange[0]) { // Checking if sum of new attempted atomic charges in molecule1
                                                       // will fall below lower limit in molrange1
                    sumTemp = 0;                       // resetting temporary sum of atomic charges
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        sumTemp +=
                            p->charge -
                            (2 * mol1.min[i] -
                             (p->charge +
                              mol1.changeQ[i])); // temporary sum of charge moves attempted on all atoms in molecule1
                        p->charge = 2 * mol1.min[i] -
                                    (p->charge + mol1.changeQ[i]); // new attempted charge of atom i in molecule1,
                                                                   // obeying torodial boundary conditions
                    }
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        p->charge += sumTemp * mol2.ratio[i]; // new attempted charge of atom i in molecule2, obeying
                                                              // torodial boundary conditions
                    }
                }

                else if (mol1.charges >
                         mol1.molrange[1]) { // same procedure as above if statement, but if sum of new atempted charges
                                             // in molecule1 falls above upper limit in molrange1
                    sumTemp = 0;
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        sumTemp += p->charge - (2 * mol1.max[i] - (p->charge + mol1.changeQ[i]));
                        p->charge = 2 * mol1.max[i] - (p->charge + mol1.changeQ[i]);
                    }
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        p->charge += sumTemp * mol2.ratio[i];
                    }
                }

                else if (mol2.charges < mol2.molrange[0]) { // same as first if statement, but with respect to molecule2
                                                            // and its molrange
                    sumTemp = 0;
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        sumTemp += p->charge - (2 * mol2.min[i] - (p->charge + mol2.changeQ[i]));
                        p->charge = 2 * mol2.min[i] - (p->charge + mol2.changeQ[i]);
                    }
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        p->charge += sumTemp * mol1.ratio[i];
                    }
                }

                else if (mol2.charges >
                         mol2.molrange[1]) { // same as previous if statement, but if sum of new attempted charges in
                                             // molecule2 falls above upper limit in molrange2
                    sumTemp = 0;
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        sumTemp += p->charge - (2 * mol2.max[i] - (p->charge + mol2.changeQ[i]));
                        p->charge = 2 * mol2.max[i] - (p->charge + mol2.changeQ[i]);
                    }
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        p->charge += sumTemp * mol1.ratio[i];
                    }
                }

                else { // in case no boundaries were crossed, i.e. all new charges lies within their respective ranges
                    for (i = 0; i < mol1.numOfAtoms; i++) {
                        auto p = git1->begin() + i;
                        p->charge += mol1.changeQ[i];
                    }
                    for (i = 0; i < mol2.numOfAtoms; i++) {
                        auto p = git2->begin() + i;
                        p->charge += mol2.changeQ[i];
                    }
                }
                mol1.cdata.all = true;               // change all atoms in molecule1
                mol2.cdata.all = true;               // change all atoms in molecule2
                change.groups.push_back(mol1.cdata); // add to list of moved groups
                change.groups.push_back(mol2.cdata); // add to list of moved groups

            } else
                deltaq = 0;
        }
    }
}

void ChargeTransfer::_accept(Change &) { msqd += deltaq * deltaq; }
void ChargeTransfer::_reject(Change &) { msqd += 0; }
ChargeTransfer::ChargeTransfer(Space &spc) : spc(spc) {

    name = "chargetransfer";
    repeat = -1; // meaning repeat N times
    mol1.cdata.internal = true;
    mol2.cdata.internal = true;
    // cdata1.atoms.resize(numOfAtoms1);
    // cdata2.atoms.resize(numOfAtoms2);
}

void QuadrantJump::_to_json(json &j) const {
    j = {{"dir", dir},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molecules[molid].name}};
    _roundjson(j, 3);
}
void QuadrantJump::_from_json(const json &j) {
    assert(!molecules.empty());
    try {
        std::string molname = j.at("molecule");
        auto it = findName(molecules, molname);
        if (it == molecules.end())
            throw std::runtime_error("unknown molecule '" + molname + "'");
        molid = it->id();
        dir = j.value("dir", Point(1, 1, 1));
        index = j.value("index", decltype(index)());
        if (repeat < 0) {
            auto v = spc.findMolecules(molid);
            repeat = std::distance(v.begin(), v.end());
        }
    } catch (std::exception &e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}
void QuadrantJump::_move(Change &change) {
    assert(molid >= 0);
    assert(!spc.groups.empty());
    assert(spc.geo.getVolume() > 0);

    _sqd = 0.0;

    // pick random group from the system matching molecule type
    // TODO: This can be slow -- implement look-up-table in Space
    auto mollist = spc.findMolecules(molid, Space::ACTIVE); // list of molecules w. 'molid'
    if (not ranges::cpp20::empty(mollist)) {
        auto it = slump.sample(mollist.begin(), mollist.end());
        if (not it->empty()) {
            assert(it->id == molid);
            Point oldcm = it->cm;
            if (index.size() == 2) {
                auto cm_O = Geometry::massCenter(spc.p.begin() + index[0], spc.p.begin() + index[1] + 1,
                                                 spc.geo.getBoundaryFunc());
                it->translate(-2 * spc.geo.vdist(oldcm, cm_O).cwiseProduct(dir.cast<double>()),
                              spc.geo.getBoundaryFunc());
            } else {
                it->translate(-2 * oldcm.cwiseProduct(dir.cast<double>()), spc.geo.getBoundaryFunc());
            }
            _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
            Change::data d;
            d.index = Faunus::distance(spc.groups.begin(), it); // integer *index* of moved group
            d.all = true;                                       // *all* atoms in group were moved
            change.groups.push_back(d);                         // add to list of moved groups

            assert(spc.geo.sqdist(it->cm, Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(),
                                                               -it->cm)) < 1e-9);
        }
    } else
        faunus_logger->warn("{0}: no molecules found", name);
}
QuadrantJump::QuadrantJump(Space &spc) : spc(spc) {
    name = "quadrantjump";
    repeat = -1; // meaning repeat N times
}
void AtomicSwapCharge::_to_json(json &j) const {
    j = {{"pH", pH},
         {"pka", pKa},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molname}};
    _roundjson(j, 3);
}
void AtomicSwapCharge::_from_json(const json &j) {
    assert(!molecules.empty());
    try {
        molname = j.at("molecule");
        auto it = findName(molecules, molname);
        if (it == molecules.end())
            throw std::runtime_error("unknown molecule '" + molname + "'");
        molid = it->id();
        pH = j.at("pH").get<double>();
        pKa = j.at("pKa").get<double>();
        if (repeat < 0) {
            auto v = spc.findMolecules(molid);
            repeat = std::distance(v.begin(), v.end()); // repeat for each molecule...
            if (repeat > 0)
                repeat = repeat * v.front().size(); // ...and for each atom
        }
    } catch (std::exception &e) {
        throw std::runtime_error(name + ": "s + e.what());
    }
}
typename Space::Tpvec::iterator AtomicSwapCharge::randomAtom() {
    assert(molid >= 0);
    auto mollist = spc.findMolecules(molid); // all `molid` groups
    if (not ranges::cpp20::empty(mollist)) {
        auto git = slump.sample(mollist.begin(), mollist.end()); // random molecule iterator
        if (!git->empty()) {
            auto p = slump.sample(git->begin(), git->end());         // random particle iterator
            cdata.index = Faunus::distance(spc.groups.begin(), git); // integer *index* of moved group
            cdata.atoms[0] = std::distance(git->begin(), p);         // index of particle rel. to group
            return p;
        }
    }
    return spc.p.end();
}
void AtomicSwapCharge::_move(Change &change) {
    _sqd = 0.0;
    auto p = randomAtom();
    if (p != spc.p.end()) {
        // auto &g = spc.groups[cdata.index];
        double oldcharge = p->charge;
        p->charge = fabs(oldcharge - 1);
        _sqd = fabs(oldcharge - 1) - oldcharge;
        change.groups.push_back(cdata);   // add to list of moved groups
        _bias = _sqd * (pH - pKa) * ln10; // one may add bias here...
    }
}
double AtomicSwapCharge::bias(Change &, double, double) { return _bias; }
void AtomicSwapCharge::_accept(Change &) { msqd += _sqd; }
void AtomicSwapCharge::_reject(Change &) { msqd += 0; }
AtomicSwapCharge::AtomicSwapCharge(Space &spc) : spc(spc) {
    name = "swapcharge";
    repeat = -1; // meaning repeat N times
    cdata.atoms.resize(1);
    cdata.internal = true;
}
void TranslateRotate::_to_json(json &j) const {
    j = {{"dir", dir},
         {"dp", dptrans},
         {"dprot", dprot},
         {"dirrot", dirrot},
         {"molid", molid},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molecules[molid].name}};
    _roundjson(j, 3);
}
void TranslateRotate::_from_json(const json &j) {
    assert(!molecules.empty());
    try {
        std::string molname = j.at("molecule");
        auto it = findName(molecules, molname);
        if (it == molecules.end())
            throw std::runtime_error("unknown molecule '" + molname + "'");
        molid = it->id();
        dir = j.value("dir", Point(1, 1, 1));
        dprot = j.at("dprot");
        dirrot = j.value("dirrot", Point(0, 0, 0)); // predefined axis of rotation
        dptrans = j.at("dp");
        if (repeat < 0) {
            auto v = spc.findMolecules(molid);
            repeat = std::distance(v.begin(), v.end());
        }
    } catch (std::exception &e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}
void TranslateRotate::_move(Change &change) {
    assert(molid >= 0);
    assert(!spc.groups.empty());
    assert(spc.geo.getVolume() > 0);

    _sqd = 0;

    // pick random group from the system matching molecule type
    // TODO: This can be slow -- implement look-up-table in Space
    auto mollist = spc.findMolecules(molid, Space::ACTIVE); // list of molecules w. 'molid'
    if (not ranges::cpp20::empty(mollist)) {
        auto it = slump.sample(mollist.begin(), mollist.end());
        if (not it->empty()) {
            assert(it->id == molid);

            if (dptrans > 0) { // translate
                Point oldcm = it->cm;
                Point dp = ranunit(slump, dir) * dptrans * slump();

                it->translate(dp, spc.geo.getBoundaryFunc());
                _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
            }

            if (dprot > 0) { // rotate
                Point u = ranunit(slump);
                if (dirrot.count() > 0)
                    u = dirrot;
                double angle = dprot * (slump() - 0.5);
                Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, u));
                it->rotate(Q, spc.geo.getBoundaryFunc());
            }

            if (dptrans > 0 || dprot > 0) { // define changes
                Change::data d;
                d.index = Faunus::distance(spc.groups.begin(), it); // integer *index* of moved group
                d.all = true;                                       // *all* atoms in group were moved
                change.groups.push_back(d);                         // add to list of moved groups
            }
#ifndef NDEBUG
            // check if mass center is correctly moved and can be re-calculated
            Point cm_recalculated = Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(), -it->cm);
            double should_be_small = spc.geo.sqdist(it->cm, cm_recalculated);
            if (should_be_small > 1e-6) {
                std::cerr << "cm recalculated: " << cm_recalculated.transpose() << "\n";
                std::cerr << "cm in group:     " << it->cm.transpose() << "\n";
                assert(false);
            }
//            assert(spc.geo.sqdist(it->cm, Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(),
//                                                               -it->cm)) < 1e-6);////////////////
#endif
        }
    }
}
TranslateRotate::TranslateRotate(Space &spc) : spc(spc) {
    name = "moltransrot";
    repeat = -1; // meaning repeat N times
}

void SmartTranslateRotate::_to_json(json &j) const {
    j = {{"Number of counts inside geometry", cntInner},
         {"Number of counts outside geometry", cnt - cntInner},
         {"dir", dir},
         {"dp", dptrans},
         {"dprot", dprot},
         {"p", p},
         {"origo", origo},
         {"rx", r_x},
         {"ry", r_y},
         {"molid", molid},
         {"refid", refid1},
         {"refid", refid2},
         {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
         {"molecule", molecules[molid].name},
         {"ref1", atoms[refid1].name},
         {"ref2", atoms[refid2].name}};
    _roundjson(j, 3);
}
void SmartTranslateRotate::_from_json(const json &j) {
    assert(!molecules.empty());
    try {
        std::string molname = j.at("molecule");
        std::string refname1 = j.at("ref1");
        std::string refname2 = j.at("ref2");
        auto it = findName(molecules, molname);
        auto ref1 = findName(atoms, refname1);
        auto ref2 = findName(atoms, refname2);
        if (it == molecules.end())
            throw std::runtime_error("unknown molecule '" + molname + "'");
        if (ref1 == atoms.end())
            throw std::runtime_error("unknown reference atom '" + refname1 + "'");
        if (ref2 == atoms.end())
            throw std::runtime_error("unknown reference atom '" + refname2 + "'");
        molid = it->id();
        refid1 = ref1->id();
        refid2 = ref2->id();
        dir = j.value("dir", Point(1, 1, 1));
        dprot = j.at("dprot");
        dptrans = j.at("dp");
        p = j.at("p");
        r_x = j.at("rx");  // length of ellipsoidal radius along axis connecting reference atoms (in Å)
        r_y = j.at("ry");  // length of ellipsoidal radius perpendicular to axis connecting reference atoms (in Å)
        rsd = j.at("rsd"); // threshold for relative standard deviation of molecules inside geometry. When it goes below
                           // this value, a constant bias is used

        if (repeat < 0) {
            auto v = spc.findMolecules(molid);
            repeat = std::distance(v.begin(), v.end());
        }
    } catch (std::exception &e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}

void SmartTranslateRotate::_move(Change &change) {
    assert(molid >= 0);
    assert(!spc.groups.empty());
    assert(spc.geo.getVolume() > 0);
    _bias = 0.0;
    _sqd = 0.0;

    // pick random group from the system matching molecule type
    // TODO: This can be slow -- implement look-up-table in Space
    auto mollist = spc.findMolecules(molid, Space::ACTIVE); // list of molecules w. 'molid'
    auto reflist1 = spc.findAtoms(refid1);                  // list of atoms w. 'refid1'
    auto reflist2 = spc.findAtoms(refid2);                  // list of atoms w. 'refid2'
    if (not ranges::cpp20::empty(mollist)) {
        auto it = slump.sample(mollist.begin(), mollist.end()); // chosing random molecule in group of type molname
        auto ref1 = slump.sample(reflist1.begin(), reflist1.end());
        auto ref2 = slump.sample(reflist2.begin(), reflist2.end());
        cylAxis = spc.geo.vdist(ref2->pos, ref1->pos) * 0.5; // half vector between reference atoms
        origo = ref2->pos - cylAxis; // coordinates of middle point between reference atoms: new origo
        if (r_x < cylAxis.norm())    // checking so that a is larger than length of cylAxis
            throw std::runtime_error(
                "specified radius of ellipsoid along the axis connecting reference atoms (rx) must be larger or equal "
                "to half the distance between reference atoms. Specified radius is " +
                std::to_string(r_x) + " Å whereas half the distance between reference atoms is " +
                std::to_string(cylAxis.norm()) + "Å");

        if (not it->empty()) { // checking so that molecule exists
            assert(it->id == molid);

            randNbr = slump();                   // assigning random number in range [0,1]
            molV = spc.geo.vdist(it->cm, origo); // vector between selected molecule and center of geometry
            cosTheta = molV.dot(cylAxis) / molV.norm() /
                       cylAxis.norm(); // cosinus of angle between coordinate vector of selected molecule and axis
                                       // connecting reference atoms
            theta = acos(
                cosTheta); // angle between coordinate vector of selected molecule and axis connecting reference atoms
            x = cosTheta * molV.norm();   // x coordinate of selected molecule with respect to center of geometry (in
                                          // plane including vectors molV and cylAxis)
            y = sin(theta) * molV.norm(); // y coordinate of selected molecule with respect to center of geometry (in
                                          // plane including vectors molV and cylAxis)
            coord =
                x * x / (r_x * r_x) + y * y / (r_y * r_y); // calculating normalized coordinate with respect to
                                                           // dimensions of geometry (> 1.0 -> outside, < 1.0 -> inside)

            if (not(coord > 1.0 && p < randNbr)) {

                if (coord <= 1.0)
                    cntInner += 1; // counting number of times a molecule is found inside geometry

                cnt += 1; // total number of counts

                if (findBias ==
                    true) { // continuing to adjust bias according to number of molecules inside and outside geometry
                    countNin = 0.0;  // counter keeping track of number of molecules inside geometry
                    countNout = 0.0; // counter keeping track of number of molecules outside geometry
                    Ntot = 0.0;      // total number of particles
                    for (auto &g : mollist) {
                        Ntot += 1.0;
                        molV = spc.geo.vdist(g.cm, origo);
                        cosTheta = molV.dot(cylAxis) / molV.norm() / cylAxis.norm();
                        theta = acos(cosTheta);
                        x = cosTheta * molV.norm();
                        y = sin(theta) * molV.norm();
                        coordTemp = x * x / (r_x * r_x) + y * y / (r_y * r_y);
                        if (coordTemp <= 1.0)
                            countNin += 1.0;
                        else
                            countNout += 1.0;
                    }

                    countNin_avg +=
                        countNin; // appending number of molecules inside geometry (since it has type Average)
                    countNout_avg +=
                        countNout; // appending number of molecules outside geometry (since it has type Average)

                    if (cnt % 100 == 0) {
                        countNin_avgBlocks +=
                            countNin_avg.avg(); // appending average number of molecules inside geometry (type Average)
                        countNout_avgBlocks +=
                            countNout_avg
                                .avg(); // appending average number of molecules outside geometry (type Average)
                    }

                    if (cnt % 100000 == 0) {
                        Nin = countNin_avgBlocks.avg(); // block average number of molecules inside geometry
                        if (countNin_avgBlocks.stdev() / Nin <
                            rsd) { // if block standard deviation is below specified threshold
                            std::cout << "Bias found with rsd = " << countNin_avgBlocks.stdev() / Nin << " < " << rsd
                                 << "\n\n";
                            std::cout << "Average # of water molecules inside sphere: " << Nin << "\n";
                            findBias = false; // stop updating bias, use constant value
                        }
                    }
                }
                if (dptrans > 0) { // translate
                    Point oldcm = it->cm;
                    Point dp = ranunit(slump, dir) * dptrans * slump();

                    it->translate(dp, spc.geo.getBoundaryFunc());
                    _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
                }

                if (dprot > 0) { // rotate
                    Point u = ranunit(slump);
                    double angle = dprot * (slump() - 0.5);
                    Eigen::Quaterniond Q(Eigen::AngleAxisd(angle, u));
                    it->rotate(Q, spc.geo.getBoundaryFunc());
                }

                if (dptrans > 0 || dprot > 0) { // define changes
                    Change::data d;
                    d.index = Faunus::distance(spc.groups.begin(), it); // integer *index* of moved group
                    d.all = true;                                       // *all* atoms in group were moved
                    change.groups.push_back(d);                         // add to list of moved groups
                }
                assert(spc.geo.sqdist(it->cm, Geometry::massCenter(it->begin(), it->end(), spc.geo.getBoundaryFunc(),
                                                                   -it->cm)) < 1e-6);
                molV = spc.geo.vdist(it->cm, origo);
                cosTheta = molV.dot(cylAxis) / molV.norm() / cylAxis.norm();
                theta = acos(cosTheta);
                x = cosTheta * molV.norm();
                y = sin(theta) * molV.norm();
                coordNew = x * x / (r_x * r_x) + y * y / (r_y * r_y);

                if (findBias == true) {                 // if using constantly updated bias
                    if (coord <= 1.0 && coordNew > 1.0) // if molecule goes from inside to outside geometry
                        _bias = -log(p / (1 - (1 - p) / (p * Ntot +
                                                         (1 - p) * countNin))); // use corresponding bias, based on this
                                                                                // cycle's number of molecules inside
                    else if (coord > 1.0 && coordNew <= 1.0) // if molecules goes from outside to inside geometry
                        _bias = -log(1 / (1 + (1 - p) / (p * Ntot +
                                                         (1 - p) * countNin))); // use corresponding bias, based on this
                                                                                // cycle's number of molecules inside
                }

                else {                                  // if constant bias has been found
                    if (coord <= 1.0 && coordNew > 1.0) // if molecule goes from inside to outside geometry
                        _bias = -log(p / (1 - (1 - p) / (p * Ntot + (1 - p) * Nin))); // use corresponding bias based on
                                                                                      // average, constant value Nin
                    else if (coord > 1.0 && coordNew <= 1.0) // if molecule goes from outside to inside geometry
                        _bias = -log(1 / (1 + (1 - p) / (p * Ntot + (1 - p) * Nin))); // use corresponding bias based on
                                                                                      // average, constant value Nin
                }
            }
        }
    }
}

double SmartTranslateRotate::bias(Change &, double, double) { return _bias; }

SmartTranslateRotate::SmartTranslateRotate(Space &spc) : spc(spc) {
    name = "smartmoltransrot";
    repeat = -1; // meaning repeat N times
}

void ConformationSwap::_to_json(json &j) const {
    j = {{"molid", molid}, {"molecule", molecules[molid].name}, {"keeppos", inserter.keep_positions}};
    _roundjson(j, 3);
}
void ConformationSwap::_from_json(const json &j) {
    assert(!molecules.empty());
    try {
        std::string molname = j.at("molecule");
        inserter.keep_positions = j.value("keeppos", false);
        auto it = findName(molecules, molname);
        if (it == molecules.end())
            throw std::runtime_error("unknown molecule '" + molname + "'");
        molid = it->id();
        if (molecules[molid].conformations.size() < 2)
            throw std::runtime_error("minimum two conformations required");
        if (repeat < 0) {
            auto v = spc.findMolecules(molid);
            repeat = std::distance(v.begin(), v.end());
        }
    } catch (std::exception &e) {
        throw std::runtime_error(name + ": " + e.what());
    }
}
void ConformationSwap::_move(Change &change) {
    assert(molid >= 0);
    assert(change.empty());

    auto mollist = spc.findMolecules(molid, Space::ACTIVE); // list of molecules w. 'molid'
    if (not ranges::cpp20::empty(mollist)) {
        auto g = slump.sample(mollist.begin(), mollist.end());
        if (not g->empty()) {
            inserter.offset = g->cm;

            // Get a new conformation that should be properly wrapped around the boundaries
            // (if applicable) and have the same mass-center as "g->cm".
            Tpvec p = inserter(spc.geo, spc.p, molecules[molid]);
            if (p.size() not_eq g->size())
                throw std::runtime_error(name + ": conformation atom count mismatch");

            newconfid = molecules[molid].conformations.getLastIndex();

            std::copy(p.begin(), p.end(), g->begin()); // override w. new conformation
#ifndef NDEBUG
            // this move shouldn't move mass centers, so let's check if this is true:
            Point newcm = Geometry::massCenter(p.begin(), p.end(), spc.geo.getBoundaryFunc(), -g->cm);
            if ((newcm - g->cm).norm() > 1e-6)
                throw std::runtime_error(name + ": unexpected mass center movement");
#endif
            Change::data d;
            d.index = Faunus::distance(spc.groups.begin(), g); // integer *index* of moved group
            d.all = true;                                      // *all* atoms in group were moved
            d.internal = false;                                // we *don't* want to calculate the internal energy
            change.groups.push_back(d);                        // add to list of moved groups
        }
    }
}
void ConformationSwap::_accept(Change &change) {
    assert(change.groups.size() == 1);
    spc.groups[change.groups.front().index].confid = newconfid;
}
ConformationSwap::ConformationSwap(Space &spc) : spc(spc) {
    name = "conformationswap";
    repeat = -1; // meaning repeat n times
    inserter.dir = {0, 0, 0};
    inserter.rotate = true;
    inserter.allow_overlap = true;
}

ForceMove::ForceMove() {
    // resize forces and velocities to mathc spc.p
}

} // namespace Move

} // namespace Faunus
