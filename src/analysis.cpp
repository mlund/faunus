#include "analysis.h"

namespace Faunus {

namespace Analysis {

void to_json(json &j, const Analysisbase &base) { base.to_json(j); }

Analysisbase::~Analysisbase() {}

/*
 * This is always called by the destructor, but may be called
 * any number of times earlier.
 */
void Analysisbase::_to_disk() {}

void Analysisbase::to_disk() { _to_disk(); }

void Analysisbase::sample() {
    totstepcnt++;
    stepcnt++;
    if (stepcnt == steps) {
        stepcnt = 0;
        if (totstepcnt > nskip) {
            cnt++;
            timer.start();
            _sample();
            timer.stop();
        }
    }
}

void Analysisbase::from_json(const json &j) {
    steps = j.value("nstep", 0);
    nskip = j.value("nskip", 0);
    _from_json(j);
}

void Analysisbase::to_json(json &j) const {
    assert(not name.empty());
    auto &_j = j[name];
    _to_json(_j);
    if (cnt > 0) {
        if (timer.result() > 0.01) // only print if more than 1% of the time
            _j["relative time"] = _round(timer.result());
        _j["nstep"] = steps;
        _j["samples"] = cnt;
        if (nskip > 0)
            _j["nskip"] = nskip;
    }
    if (not cite.empty())
        _j["reference"] = cite;
}

void Analysisbase::_to_json(json &) const {}

void Analysisbase::_from_json(const json &) {}

void SystemEnergy::normalize() {
    // assert(V.cnt>0);
    double sum = ehist.sumy();
    for (auto &i : ehist.getMap()) {
        i.second = i.second / sum;
    }
}

void SystemEnergy::_sample() {
    auto ulist = energyFunc();
    double tot = std::accumulate(ulist.begin(), ulist.end(), 0.0);
    if (not std::isinf(tot)) {
        uavg += tot;
        u2avg += tot * tot;
    }
    f << cnt * steps << sep << tot;
    for (auto u : ulist)
        f << sep << u;
    f << "\n";
    // ehist(tot)++;
}

void SystemEnergy::_to_json(json &j) const {
    j = {{"file", file}, {"init", uinit}, {"final", energyFunc()}};
    if (cnt > 0) {
        j["mean"] = uavg.avg();
        j["Cv/kB"] = u2avg.avg() - std::pow(uavg.avg(), 2);
    }
    _roundjson(j, 5);
    // normalize();
    // ehist.save( "distofstates.dat" );
}

void SystemEnergy::_from_json(const json &j) {
    file = MPI::prefix + j.at("file").get<std::string>();
    if (f)
        f.close();
    f.open(file);
    if (!f)
        throw std::runtime_error(name + ": cannot open output file " + file);
    assert(!names.empty());
    std::string suffix = file.substr(file.find_last_of(".") + 1);
    if (suffix == "csv")
        sep = ",";
    else {
        sep = " ";
        f << "#";
    }
    f << "total";
    for (auto &n : names)
        f << sep << n;
    f << "\n";
}

void SaveState::_to_json(json &j) const { j["file"] = file; }

void SaveState::_sample() { writeFunc(file); }

SaveState::~SaveState() {
    if (steps == -1)
        _sample();
}

PairFunctionBase::PairFunctionBase(const json &j) { from_json(j); }

PairFunctionBase::~PairFunctionBase() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        double Vr = 1, sum = hist.sumy();
        hist.stream_decorator = [&](std::ostream &o, double r, double N) {
            if (dim == 3)
                Vr = 4 * pc::pi * std::pow(r, 2) * dr;
            else if (dim == 2) {
                Vr = 2 * pc::pi * r * dr;
                if (Rhypersphere > 0)
                    Vr = 2.0 * pc::pi * Rhypersphere * std::sin(r / Rhypersphere) * dr;
            } else if (dim == 1)
                Vr = dr;
            if (Vr > 0)
                o << r << " " << N * V / (Vr * sum) << "\n";
        };
        f << hist;
    }
}

void PairFunctionBase::_to_json(json &j) const {
    j = {{"dr", dr / 1.0_angstrom}, {"name1", name1},        {"name2", name2}, {"file", file}, {"dim", dim},
         {"slicedir", slicedir},    {"thickness", thickness}};
    if (Rhypersphere > 0)
        j["Rhyper"] = Rhypersphere;
}

void PairFunctionBase::_from_json(const json &j) {
    assertKeys(j, {"file", "name1", "name2", "dim", "dr", "Rhyper", "nstep", "nskip", "slicedir", "thickness"});
    file = j.at("file");
    name1 = j.at("name1");
    name2 = j.at("name2");
    dim = j.value("dim", 3);
    dr = j.value("dr", 0.1) * 1.0_angstrom;
    slicedir = j.value("slicedir", slicedir);
    thickness = j.value("thickness", 0);
    hist.setResolution(dr, 0);
    Rhypersphere = j.value("Rhyper", -1.0);
}

void VirtualVolume::_sample() {
    if (fabs(dV) > 1e-10) {
        double Vold = getVolume(), Uold = pot.energy(c);
        scaleVolume(Vold + dV);
        double Unew = pot.energy(c);
        scaleVolume(Vold);

        // check if energy change is too big for exp()
        double x = pc::infty, du = Unew - Uold;
        if (-du < pc::max_exp_argument)
            x = std::exp(-du);
        if (std::isinf(x)) {
            std::cerr << name + ": skipping sample event due to excessive energy likely due to overlaps.";
            cnt--; // cnt is incremented by sample() so we need to decrease
        } else {
            assert(not std::isnan(x));
            duexp += x;
#ifndef NDEBUG
            // check volume and particle positions are properly restored
            double err = std::fabs((Uold - pot.energy(c)) / Uold); // must be ~zero!
            if (std::isfinite(err))                                // catch if Uold==0
                assert(err < 1e-4);
#endif
        }
    }
}

void VirtualVolume::_from_json(const json &j) { dV = j.at("dV"); }

void VirtualVolume::_to_json(json &j) const {
    double pex = log(duexp.avg()) / dV; // excess pressure
    j = {{"dV", dV}, {"Pex/mM", pex / 1.0_mM}, {"Pex/Pa", pex / 1.0_Pa}, {"Pex/kT/" + u8::angstrom + u8::cubed, pex}};
    _roundjson(j, 5);
}

void QRtraj::_sample() { write_to_file(); }

void QRtraj::_to_json(json &j) const { j = {{"file", file}}; }

void CombinedAnalysis::sample() {
    for (auto i : this->vec)
        i->sample();
}

CombinedAnalysis::~CombinedAnalysis() {
    for (auto i : this->vec)
        i->to_disk();
}
void FileReactionCoordinate::_to_json(json &j) const {
    j = *rc;
    j["type"] = type;
    j["file"] = filename;
    j.erase("range");      // these are for penalty function
    j.erase("resolution"); // use only, so no need to show
    if (cnt > 0)
        j["average"] = avg.avg();
}
void FileReactionCoordinate::_sample() {
    if (file) {
        double val = (*rc)();
        avg += val;
        file << cnt * steps << " " << val << " " << avg.avg() << "\n";
    }
}
FileReactionCoordinate::FileReactionCoordinate(const json &j, Tspace &spc) {
    using namespace ReactionCoordinate;
    from_json(j);
    name = "reactioncoordinate";
    filename = MPI::prefix + j.at("file").get<std::string>();
    file.open(filename); // output file
    type = j.at("type").get<std::string>();
    try {
        if (type == "atom")
            rc = std::make_shared<AtomProperty>(j, spc);
        else if (type == "molecule")
            rc = std::make_shared<MoleculeProperty>(j, spc);
        else if (type == "system")
            rc = std::make_shared<SystemProperty>(j, spc);
        else if (type == "cmcm")
            rc = std::make_shared<MassCenterSeparation>(j, spc);
        if (rc == nullptr)
            throw std::runtime_error("unknown coordinate type");

    } catch (std::exception &e) {
        throw std::runtime_error("error for reaction coordinate '" + type + "': " + e.what() +
                                 usageTip["coords=[" + type + "]"]);
    }
}

void WidomInsertion::_sample() {
    if (!change.empty()) {
        Tpvec pin;
        auto &g = spc.groups.at(change.groups.at(0).index);
        assert(g.empty());
        g.resize(g.capacity()); // active group
        for (int i = 0; i < ninsert; ++i) {
            pin = rins(spc.geo, spc.p, molecules.at(molid));
            if (!pin.empty()) {
                if (absolute_z)
                    for (auto &p : pin)
                        p.pos.z() = std::fabs(p.pos.z());

                assert(pin.size() == g.size());
                spc.geo.randompos(pin[0].pos, random);
                spc.geo.randompos(pin[1].pos, random);

                std::copy(pin.begin(), pin.end(), g.begin()); // copy into ghost group
                if (!g.atomic)                                // update molecular mass-center
                    g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.begin()->pos);

                expu += exp(-pot->energy(change)); // widom average
            }
        }
        g.resize(0); // deactive molecule
    }
}
void WidomInsertion::_to_json(json &j) const {
    double excess = -std::log(expu.avg());
    j = {{"dir", rins.dir},
         {"molecule", molname},
         {"insertions", expu.cnt},
         {"absz", absolute_z},
         {u8::mu + "/kT", {{"excess", excess}}}};
}
void WidomInsertion::_from_json(const json &j) {
    ninsert = j.at("ninsert");
    molname = j.at("molecule");
    absolute_z = j.value("absz", false);
    rins.dir = j.value("dir", Point({1, 1, 1}));

    auto it = findName(molecules, molname); // loop for molecule in topology
    if (it != molecules.end()) {
        molid = it->id();
        auto m = spc.findMolecules(molid, Tspace::INACTIVE); // look for molecules in space
        if (size(m) > 0) {                                   // did we find any?
            if (m.begin()->size() == 0) {                    // pick the first and check if it's really inactive
                change.clear();
                Change::data d;                                    // construct change object
                d.index = distance(spc.groups.begin(), m.begin()); // group index
                d.all = true;
                d.internal = m.begin()->atomic;
                change.groups.push_back(d); // add to change object
                return;
            }
        }
    }
    throw std::runtime_error(name + ": no inactive '" + molname + "' groups found");
}
void Density::_sample() {
    // count atom and groups of individual id's
    Nmol.clear();
    Natom.clear();

    // make sure all atom counts are initially zero
    for (auto &g : spc.groups) {
        if (g.atomic)
            for (auto p = g.begin(); p < g.trueend(); ++p)
                Natom[p->id] = 0;
        else
            Nmol[g.id] = 0;
    }

    double V = spc.geo.getVolume();
    Vavg += V;
    Lavg += std::cbrt(V);
    invVavg += 1 / V;

    for (auto &g : spc.groups)
        if (g.atomic) {
            for (auto &p : g)
                Natom[p.id]++;
            atmdhist[g.id](g.size())++;
        } else if (not g.empty())
            Nmol[g.id]++;

    for (auto &i : Nmol) {
        rho_mol[i.first] += i.second / V;
        moldhist[i.first](i.second)++;
    }

    for (auto &i : Natom)
        rho_atom[i.first] += i.second / V;

    if (Faunus::reactions.size() > 0) { // in case of reactions involving atoms (swap moves)
        for (auto &rit : reactions) {
            for (auto pid : rit._prodid_a) {
                auto atomlist = spc.findAtoms(pid.first);
                swpdhist[pid.first](size(atomlist))++;
            }
            for (auto rid : rit._reagid_a) {
                auto atomlist = spc.findAtoms(rid.first);
                swpdhist[rid.first](size(atomlist))++;
            }
        }
    }
}
void Density::_to_json(json &j) const {
    using namespace u8;
    j[bracket("V")] = Vavg.avg();
    j[bracket("1/V")] = invVavg.avg();
    j[bracket(cuberoot + "V")] = Lavg.avg();
    j[cuberoot + bracket("V")] = std::cbrt(Vavg.avg());

    auto &_j = j["atomic"];
    for (auto &i : rho_atom)
        if (i.second.cnt > 0)
            _j[atoms.at(i.first).name] = json({{"c/M", i.second.avg() / 1.0_molar}});

    auto &_jj = j["molecular"];
    for (auto &i : rho_mol)
        if (i.second.cnt > 0)
            _jj[molecules.at(i.first).name] = json({{"c/M", i.second.avg() / 1.0_molar}});
    _roundjson(j, 4);
}
Density::Density(const json &j, Tspace &spc) : spc(spc) {
    from_json(j);
    name = "density";
    for (auto &m : molecules) {
        if (m.atomic)
            atmdhist[m.id()].setResolution(1, 0);
        else
            moldhist[m.id()].setResolution(1, 0);
    }
    if (Faunus::reactions.size() > 0) { // in case of reactions involving atoms (swap moves)
        for (auto &rit : reactions) {
            for (auto pid : rit._prodid_a) {
                swpdhist[pid.first].setResolution(1, 0);
            }
            for (auto rid : rit._reagid_a) {
                swpdhist[rid.first].setResolution(1, 0);
            }
        }
    }
}
Density::~Density() {
    for (auto &m : atmdhist) { // atomic molecules
        std::string file = "rho-"s + molecules.at(m.first).name + ".dat";
        std::ofstream f(MPI::prefix + file);
        if (f) {
            m.second.stream_decorator = [&](std::ostream &o, int N, double samplings) {
                double sum = m.second.sumy();
                if (samplings > 0)
                    o << N << " " << samplings << " " << samplings / sum << "\n";
            };
            f << "# N samplings P\n" << m.second;
        }
    }
    for (auto &m : moldhist) { // polyatomic molecules
        std::string file = "rho-"s + molecules.at(m.first).name + ".dat";
        std::ofstream f(MPI::prefix + file);
        if (f) {
            m.second.stream_decorator = [&](std::ostream &o, int N, double samplings) {
                double sum = m.second.sumy();
                if (samplings > 0)
                    o << N << " " << samplings << " " << samplings / sum << "\n";
            };
            f << "# N samplings P\n" << m.second;
        }
    }
    if (Faunus::reactions.size() > 0) { // in case of reactions involving atoms (swap moves)
        for (auto &rit : reactions) {
            for (auto pid : rit._prodid_a) {
                std::string file = "rho-"s + atoms.at(pid.first).name + ".dat";
                std::ofstream f(MPI::prefix + file);
                if (f) {
                    swpdhist.at(pid.first).stream_decorator = [&](std::ostream &o, int N, double samplings) {
                        double sum = swpdhist.at(pid.first).sumy();
                        if (samplings > 0)
                            o << N << " " << samplings << " " << samplings / sum << "\n";
                    };
                    f << "# N samplings P\n" << swpdhist.at(pid.first);
                }
            }
            for (auto rid : rit._reagid_a) {
                std::string file = "rho-"s + atoms.at(rid.first).name + ".dat";
                std::ofstream f(MPI::prefix + file);
                if (f) {
                    swpdhist.at(rid.first).stream_decorator = [&](std::ostream &o, int N, double samplings) {
                        double sum = swpdhist.at(rid.first).sumy();
                        if (samplings > 0)
                            o << N << " " << samplings << " " << samplings / sum << "\n";
                    };
                    f << "# N samplings P\n" << swpdhist.at(rid.first);
                }
            }
        }
    }
}
} // namespace Analysis
} // namespace Faunus
