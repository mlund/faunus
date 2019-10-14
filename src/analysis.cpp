#include "analysis.h"
#include "move.h"
#include "energy.h"
#include "reactioncoordinate.h"
#include "multipole.h"
#include "aux/iteratorsupport.h"
#include "aux/eigensupport.h"
#include "spdlog/spdlog.h"

#include <iomanip>
#include <iostream>

namespace Faunus {

namespace Analysis {

void to_json(json &j, const Analysisbase &base) { base.to_json(j); }

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
SystemEnergy::SystemEnergy(const json &j, Energy::Hamiltonian &pot) {
    for (auto i : pot.vec)
        names.push_back(i->name);
    name = "systemenergy";
    from_json(j);
    energyFunc = [&pot]() {
        Change change;
        change.all = true;
        std::vector<double> u;
        u.reserve(pot.vec.size());
        for (auto i : pot.vec)
            u.push_back(i->energy(change));
        return u;
    };
    ehist.setResolution(0.25);
    auto u = energyFunc();
    uinit = std::accumulate(u.begin(), u.end(), 0.0); // initial energy
}

void SaveState::_to_json(json &j) const { j["file"] = file; }

void SaveState::_sample() { writeFunc(file); }

SaveState::~SaveState() {
    if (steps == -1)
        _sample();
}
SaveState::SaveState(const json &j, Space &spc) {
    using std::ref;
    using std::placeholders::_1;
    from_json(j);
    name = "savestate";
    steps = j.value("nstep", -1);
    saverandom = j.value("saverandom", false);
    file = j.at("file");
    std::string suffix = file.substr(file.find_last_of(".") + 1);
    file = MPI::prefix + file;

    if (suffix == "aam")
        writeFunc = std::bind([](std::string file, Space &s) { FormatAAM::save(file, s.p); }, _1, std::ref(spc));

    else if (suffix == "gro")
        writeFunc = std::bind([](std::string file, Space &s) { FormatGRO::save(file, s); }, _1, std::ref(spc));

    else if (suffix == "pqr")
        writeFunc = std::bind([](std::string file, Space &s) { FormatPQR::save(file, s.p, s.geo.getLength()); }, _1,
                              std::ref(spc));

    else if (suffix == "xyz")
        writeFunc = std::bind([](std::string file, Space &s) { FormatXYZ::save(file, s.p, s.geo.getLength()); }, _1,
                              std::ref(spc));

    else if (suffix == "json") // JSON state file
        writeFunc = [&spc, this](const std::string &file) {
            std::ofstream f(file);
            if (f) {
                json j;
                Faunus::to_json(j, spc);
                if (this->saverandom) {
                    j["random-move"] = Move::Movebase::slump;
                    j["random-global"] = Faunus::random;
                }
                f << std::setw(2) << j;
            }
        };

    else if (suffix == "ubj") // Universal Binary JSON state file
        writeFunc = [&spc, this](const std::string &file) {
            std::ofstream f(file, std::ios::binary);
            if (f) {
                json j;
                Faunus::to_json(j, spc);
                if (this->saverandom) {
                    j["random-move"] = Move::Movebase::slump;
                    j["random-global"] = Faunus::random;
                }
                auto v = json::to_ubjson(j); // json --> binary
                f.write((const char *)v.data(), v.size() * sizeof(decltype(v)::value_type));
            }
        };

    if (writeFunc == nullptr)
        throw std::runtime_error("unknown file extension for '" + file + "'");
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

PairAngleFunctionBase::PairAngleFunctionBase(const json &j) : PairFunctionBase(j) { from_json(j); }

PairAngleFunctionBase::~PairAngleFunctionBase() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        hist2.stream_decorator = [&](std::ostream &o, double r, double N) {
            o << r << " " << N << "\n";
        };
        f << hist2;
    }
    file = file+".hist.dat"; // make sure that the file is not overwritten by base-destructor
}

void PairAngleFunctionBase::_from_json(const json &) { hist2.setResolution(dr, 0); }

void VirtualVolume::_sample() {
    if (fabs(dV) > 1e-10) {
        // store old volume and energy
        double Vold = getVolume(), Uold = pot.energy(c);
        // scale entire system to new volume
        scaleVolume(Vold + dV);
        double Unew = pot.energy(c);
        // restore saved system
        scaleVolume(Vold);

        // check if energy change is too big for exp()
        double x = pc::infty;
        double du = Unew - Uold; // system energy change
        if (-du < pc::max_exp_argument)
            x = std::exp(-du);
        if (std::isinf(x)) {
            faunus_logger->warn("{0}: skipping sample event due to excessive energy likely due to overlaps.", name);
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
VirtualVolume::VirtualVolume(const json &j, Space &spc, Energy::Energybase &pot) : pot(pot) {
    from_json(j);
    c.dV = true;
    c.all = true;
    name = "virtualvolume";
    cite = "doi:10.1063/1.472721";
    getVolume = [&spc]() { return spc.geo.getVolume(); };
    scaleVolume = [&spc](double Vnew) { spc.scaleVolume(Vnew); };
}

void QRtraj::_sample() { write_to_file(); }

void QRtraj::_to_json(json &j) const { j = {{"file", file}}; }
QRtraj::QRtraj(const json &j, Space &spc) {
    from_json(j);
    name = "qrfile";
    file = j.value("file", "qrtraj.dat"s);
    f.open(MPI::prefix + file);
    if (not f)
        throw std::runtime_error("error opening "s + file);
    f.precision(6);
    write_to_file = [&groups = spc.groups, &f = f]() {
        for (auto &g : groups) {
            for (auto it = g.begin(); it != g.trueend(); ++it) { // loop over *all* particles
                if (it < g.end())
                    f << it->charge << " " << atoms[it->id].sigma * 0.5 << " ";
                else
                    f << "0 0 "; // zero charge and radii for inactive particles
            }
        }
        f << "\n";               // newline for every frame
    };
}

void CombinedAnalysis::sample() {
    for (auto &ptr : this->vec)
        ptr->sample();
}

CombinedAnalysis::~CombinedAnalysis() {
    // this is really a hack; the constructor should not be in charge of this
    for (auto &ptr : this->vec)
        ptr->to_disk();
}

CombinedAnalysis::CombinedAnalysis(const json &j, Space &spc, Energy::Hamiltonian &pot) {
    if (j.is_array()) {
        for (auto &m : j) {
            for (auto it = m.begin(); it != m.end(); ++it) {
                if (it->is_object()) {
                    try {
                        size_t oldsize = this->vec.size();
                        if (it.key() == "atomprofile")
                            emplace_back<AtomProfile>(it.value(), spc);
                        else if (it.key() == "atomrdf")
                            emplace_back<AtomRDF>(it.value(), spc);
                        else if (it.key() == "atomdipdipcorr")
                            emplace_back<AtomDipDipCorr>(it.value(), spc);
                        else if (it.key() == "density")
                            emplace_back<Density>(it.value(), spc);
                        else if (it.key() == "chargefluctuations")
                            emplace_back<ChargeFluctuations>(it.value(), spc);
                        else if (it.key() == "molrdf")
                            emplace_back<MoleculeRDF>(it.value(), spc);
                        else if (it.key() == "multipole")
                            emplace_back<Multipole>(it.value(), spc);
                        else if (it.key() == "multipoledist")
                            emplace_back<MultipoleDistribution>(it.value(), spc);
                        else if (it.key() == "polymershape")
                            emplace_back<PolymerShape>(it.value(), spc);
                        else if (it.key() == "qrfile")
                            emplace_back<QRtraj>(it.value(), spc);
                        else if (it.key() == "reactioncoordinate")
                            emplace_back<FileReactionCoordinate>(it.value(), spc);
                        else if (it.key() == "sanity")
                            emplace_back<SanityCheck>(it.value(), spc);
                        else if (it.key() == "savestate")
                            emplace_back<SaveState>(it.value(), spc);
                        else if (it.key() == "scatter")
                            emplace_back<ScatteringFunction>(it.value(), spc);
                        else if (it.key() == "sliceddensity")
                            emplace_back<SlicedDensity>(it.value(), spc);
                        else if (it.key() == "systemenergy")
                            emplace_back<SystemEnergy>(it.value(), pot);
                        else if (it.key() == "virtualvolume")
                            emplace_back<VirtualVolume>(it.value(), spc, pot);
                        else if (it.key() == "widom")
                            emplace_back<WidomInsertion>(it.value(), spc, pot);
                        else if (it.key() == "xtcfile")
                            emplace_back<XTCtraj>(it.value(), spc);
                        // additional analysis go here...

                        if (this->vec.size() == oldsize)
                            throw std::runtime_error("unknown analysis: "s + it.key());

                    } catch (std::exception &e) {
                        throw std::runtime_error(e.what() + usageTip[it.key()]);
                    }
                }
            }
        }
    }
}

void FileReactionCoordinate::_to_json(json &j) const {
    json rcjson = *rc; // envoke to_json(...)
    if (rcjson.count(type) == 0)
        throw std::runtime_error("error writing json for reaction coordinate");
    j = rcjson[type];
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

FileReactionCoordinate::FileReactionCoordinate(const json &j, Space &spc) {
    from_json(j);
    name = "reactioncoordinate";
    filename = MPI::prefix + j.at("file").get<std::string>();
    file.open(filename); // output file
    type = j.at("type").get<std::string>();
    rc = ReactionCoordinate::createReactionCoordinate({{type, j}}, spc);
}

void WidomInsertion::_sample() {
    if (!change.empty()) {
        ParticleVector pin;
        auto &g = spc.groups.at(change.groups.at(0).index);
        assert(g.empty() && g.capacity() > 0);
        g.resize(g.capacity()); // active group
        for (int i = 0; i < ninsert; ++i) {
            pin = rins(spc.geo, spc.p, molecules.at(molid));
            if (not pin.empty()) {
                if (absolute_z) {
                    for (auto &p : pin)
                        p.pos.z() = std::fabs(p.pos.z());
                }

                assert(pin.size() == g.size());

                std::copy(pin.begin(), pin.end(), g.begin()); // copy into ghost group
                if (not g.atomic)                             // update molecular mass-center
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
        auto m = spc.findMolecules(molid, Space::INACTIVE);  // look for inactive molecules in space
        if (size(m) > 0) {                                   // did we find any?
            if (m.begin()->size() == 0) {                    // pick the first and check if it's really inactive
                if (m.begin()->capacity() > 0) {             // and it must have a non-zero capacity
                    change.clear();
                    Change::data d;                                    // construct change object
                    d.index = Faunus::distance(spc.groups.begin(), m.begin()); // group index
                    d.all = true;
                    d.internal = m.begin()->atomic; // calculate internal energy of non-molecular groups only
                    change.groups.push_back(d);     // add to change object
                    return;
                }
            }
        }
    }
    throw std::runtime_error(name + ": no inactive '" + molname + "' groups found");
}

WidomInsertion::WidomInsertion(const json &j, Space &spc, Energy::Hamiltonian &pot) : spc(spc), pot(&pot) {
    from_json(j);
    name = "widom";
    cite = "doi:10/dkv4s6";
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

    for (auto &g : spc.groups) {
        if (g.atomic) {
            for (auto &p : g)
                Natom[p.id]++;
            atmdhist[g.id](g.size())++;
        } else if (not g.empty())
            Nmol[g.id]++;
    }

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
Density::Density(const json &j, Space &spc) : spc(spc) {
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
void SanityCheck::_sample() {
    // loop over all groups
    for (auto &g : spc.groups) {
        // check if particles are inside container
        for (auto &i : g) // loop over active particles
            if (spc.geo.collision(i.pos))
                throw std::runtime_error("step "s + std::to_string(cnt) + ": index " +
                                         std::to_string(&i - &(*g.begin())) + " of group " +
                                         std::to_string(std::distance(spc.groups.begin(), spc.findGroupContaining(i))) +
                                         " outside container");

        // The groups must exactly contain all particles in `p`
        size_t i = 0;
        for (auto &g : spc.groups)
            for (auto it = g.begin(); it != g.trueend(); ++it)
                if (&*it != &spc.p.at(i++))
                    throw std::runtime_error("group vector out of sync");
        assert(i == spc.p.size());

        // check if molecular mass centers are correct
        if (not g.atomic)
            if (not g.empty()) {
                Point cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.cm);
                double sqd = spc.geo.sqdist(g.cm, cm);
                if (sqd > 1e-6) {
                    std::cerr << "step:      " << cnt << std::endl
                              << "molecule:  " << &g - &*spc.groups.begin() << std::endl
                              << "dist:      " << sqrt(sqd) << std::endl
                              << "g.cm:      " << g.cm.transpose() << std::endl
                              << "actual cm: " << cm.transpose() << std::endl;
                    FormatPQR::save(MPI::prefix + "sanity-" + std::to_string(cnt) + ".pqr", spc.p, spc.geo.getLength());
                    throw std::runtime_error("mass center-out-of-sync");
                }
            }
    }
}
SanityCheck::SanityCheck(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "sanity";
    steps = j.value("nstep", -1);
}
void AtomRDF::_sample() {
    V += spc.geo.getVolume(dim);
    auto active = spc.activeParticles();
    for (auto i = active.begin(); i != active.end(); ++i)
        for (auto j = i; ++j != active.end();)
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                Point rvec = spc.geo.vdist(i->pos, j->pos);
                if (slicedir.sum() > 0) {
                    if (rvec.cwiseProduct(slicedir.cast<double>()).norm() < thickness) {
                        // rvec = rvec.cwiseProduct( Point(1.,1.,1.) - slice.cast<double>() );
                        hist(rvec.norm())++;
                    }
                } else {
                    hist(rvec.norm())++;
                }
            }
}
AtomRDF::AtomRDF(const json &j, Space &spc) : PairFunctionBase(j), spc(spc) {
    name = "atomrdf";
    auto it = findName(atoms, name1);
    if (it == atoms.end())
        throw std::runtime_error("unknown atom '" + name1 + "'");
    id1 = it->id();

    it = findName(atoms, name2);
    if (it == atoms.end())
        throw std::runtime_error("unknown atom '" + name2 + "'");
    id2 = it->id();
}
void MoleculeRDF::_sample() {
    V += spc.geo.getVolume(dim);
    auto mollist1 = spc.findMolecules(id1, Space::ACTIVE);
    auto mollist2 = spc.findMolecules(id2, Space::ACTIVE);
    auto mollist = ranges::view::concat(mollist1, mollist2);
    for (auto i = mollist.begin(); i != mollist.end(); ++i) {
        for (auto j = i; ++j != mollist.end();) {
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                double r = std::sqrt(spc.geo.sqdist(i->cm, j->cm));
                hist(r)++;
            }
        }
    }
}
MoleculeRDF::MoleculeRDF(const json &j, Space &spc) : PairFunctionBase(j), spc(spc) {
    name = "molrdf";

    auto it = findName(molecules, name1);
    if (it == molecules.end())
        throw std::runtime_error(name + ": unknown molecule '" + name1 + "'\n");
    id1 = it->id();

    it = findName(molecules, name2);
    if (it == molecules.end())
        throw std::runtime_error(name + ": unknown molecule '" + name2 + "'\n");
    id2 = it->id();
    assert(id1 >= 0 && id2 >= 0);
}

void AtomDipDipCorr::_sample() {
    V += spc.geo.getVolume(dim);
    auto active = spc.activeParticles();
    for (auto i = active.begin(); i != active.end(); ++i) {
        for (auto j = i; ++j != active.end();) {
            if ((i->id == id1 && j->id == id2) || (i->id == id2 && j->id == id1)) {
                Point rvec = spc.geo.vdist(i->pos, j->pos);
                if (slicedir.sum() > 0) {
                    if (rvec.cwiseProduct(slicedir.cast<double>()).norm() < thickness) {
                        double dipdip = i->getExt().mu.dot(j->getExt().mu);
                        double r1 = rvec.norm();
                        hist2(r1) += dipdip;
                        hist(r1)++; // get g(r) for free
                    }
                } else {
                    double dipdip = i->getExt().mu.dot(j->getExt().mu);
                    double r1 = rvec.norm();
                    hist2(r1) += dipdip;
                    hist(r1)++; // get g(r) for free
                }
            }
        }
    }
}
AtomDipDipCorr::AtomDipDipCorr(const json &j, Space &spc) : PairAngleFunctionBase(j), spc(spc) {
    name = "atomdipdipcorr";
    auto it = findName(atoms, name1);
    if (it == atoms.end())
        throw std::runtime_error("unknown atom '" + name1 + "'");
    id1 = it->id();

    it = findName(atoms, name2);
    if (it == atoms.end())
        throw std::runtime_error("unknown atom '" + name2 + "'");
    id2 = it->id();
}

// =============== XTCtraj ===============

XTCtraj::XTCtraj(const json &j, Space &s) : filter([](Particle &) { return true; }), xtc(1e6), spc(s) {
    from_json(j);
    name = "xtcfile";
    assert(filter); // filter must be callable
}

void XTCtraj::_to_json(json &j) const {
    j["file"] = file;
    if (not names.empty())
        j["molecules"] = names;
}

void XTCtraj::_from_json(const json &j) {
    file = MPI::prefix + j.at("file").get<std::string>();

    // By default, *all* active and inactive groups are saved,
    // but here allow for a user defined list of molecule ids
    names = j.value("molecules", std::vector<std::string>());
    if (not names.empty()) {
        molids = Faunus::names2ids(Faunus::molecules, names); // molecule types to save
        if (not molids.empty())
            filter = [&](Particle &i) {
                for (auto &g : spc.groups)     // loop over all active and inactive groups
                    if (g.contains(i, true)) { // does group contain particle?
                        if (std::find(molids.begin(), molids.end(), g.id) != molids.end())
                            return true;
                        else
                            return false;
                    }
                return false;
            };
    };
}

void XTCtraj::_sample() {
    xtc.setbox(spc.geo.getLength()); // set box dimensions for frame

    // On some gcc/clang and certain ubuntu/macos combinations,
    // the ranges::view::filter(rng,unaryp) clears the `filter` function.
    // Using the ranges piping seem to solve the issue.
    assert(filter);
    auto particles = spc.p | ranges::view::filter(filter);
    assert(filter);
    bool rc = xtc.save(file, particles.begin(), particles.end());
    if (rc == false)
        faunus_logger->warn("error saving xtc");
}

// =============== MultipoleDistribution ===============

double MultipoleDistribution::g2g(const MultipoleDistribution::Tgroup &g1, const MultipoleDistribution::Tgroup &g2) {
    double u = 0;
    for (auto &i : g1)
        for (auto &j : g2)
            u += i.charge * j.charge / spc.geo.vdist(i.pos, j.pos).norm();
    return u;
}
void MultipoleDistribution::save() const {
    using std::setw;
    if (cnt > 0) {
        std::ofstream f(MPI::prefix + filename.c_str());
        if (f) {
            char w = 12;
            f.precision(4);
            f << "# Multipolar energies (kT/lB)\n"
              << std::left << setw(w) << "# R/AA" << std::right << setw(w) << "exact" << setw(w) << "total" << setw(w)
              << "ionion" << setw(w) << "iondip" << setw(w) << "dipdip" << setw(w) << "ionquad" << setw(w)
              << "mucorr\n";
            for (auto &i : m)
                f << std::left << setw(w) << i.first * dr << std::right << setw(w) << i.second.tot << setw(w)
                  << i.second.ii.avg() + i.second.id.avg() + i.second.dd.avg() + i.second.iq.avg() << setw(w)
                  << i.second.ii << setw(w) << i.second.id << setw(w) << i.second.dd << setw(w) << i.second.iq
                  << setw(w) << i.second.mucorr << "\n";
        }
    }
}
void MultipoleDistribution::_sample() {
    for (auto &gi : spc.findMolecules(ids[0]))     // find active molecules
        for (auto &gj : spc.findMolecules(ids[1])) // find active molecules
            if (gi != gj) {
                auto a = Faunus::toMultipole(gi, spc.geo.getBoundaryFunc());
                auto b = Faunus::toMultipole(gj, spc.geo.getBoundaryFunc());
                Point R = spc.geo.vdist(gi.cm, gj.cm);
                auto &d = m[to_bin(R.norm(), dr)];
                d.tot += g2g(gi, gj);
                d.ii += a.charge * b.charge / R.norm();
                d.id += q2mu(a.charge * b.getExt().mulen, b.getExt().mu, b.charge * a.getExt().mulen, a.getExt().mu, R);
                d.dd += mu2mu(a.getExt().mu, b.getExt().mu, a.getExt().mulen * b.getExt().mulen, R);
                d.iq += q2quad(a.charge, b.getExt().Q, b.charge, a.getExt().Q, R);
                d.mucorr += a.getExt().mu.dot(b.getExt().mu);
            }
}

void MultipoleDistribution::_to_json(json &j) const { j = {{"molecules", names}, {"file", filename}, {"dr", dr}}; }

MultipoleDistribution::MultipoleDistribution(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "Multipole Distribution";
    dr = j.value("dr", 0.2);
    filename = j.at("file").get<std::string>();
    names = j.at("molecules").get<decltype(names)>(); // molecule names
    ids = names2ids(molecules, names);                // names --> molids
    if (ids.size() != 2)
        throw std::runtime_error("specify exactly two molecules");
}

MultipoleDistribution::~MultipoleDistribution() { save(); }

// =============== PolymerShape ===============

void PolymerShape::_to_json(json &j) const {
    using namespace u8;
    json &k = j["molecules"];
    for (int i : ids)
        k[molecules[i].name] = {
            {bracket("Rg"), Rg.at(i).avg()},
            {bracket("Re"), Re.at(i).avg()},
            {bracket("Rg" + squared), Rg2.at(i).avg()},
            {bracket("Rg" + squared) + "-" + bracket("Rg") + squared, Rg2.at(i).avg() - std::pow(Rg.at(i).avg(), 2.0)},
            {bracket("Re" + squared) + "/" + bracket("Rg" + squared), Re2.at(i).avg() / Rg2.at(i).avg()},
            {rootof + bracket("Rg" + squared), sqrt(Rg2.at(i).avg())},
            {rootof + bracket("Re" + squared), sqrt(Re2.at(i).avg())},
            {rootof + bracket("Rgxyz" + squared),
             {sqrt(Rg2x.at(i).avg()), sqrt(Rg2y.at(i).avg()), sqrt(Rg2z.at(i).avg())}}};
}
Point PolymerShape::vectorgyrationRadiusSquared(typename Space::Tgroup &g) const {
    double sum = 0;
    Point t, r2(0, 0, 0);
    for (auto &i : g) {
        double mw = atoms.at(i.id).mw;
        t = i.pos - g.cm;
        spc.geo.boundary(t);
        r2.x() += mw * t.x() * t.x();
        r2.y() += mw * t.y() * t.y();
        r2.z() += mw * t.z() * t.z();
        sum += mw;
    }
    assert(sum > 0 && "Zero molecular weight not allowed.");
    return r2 * (1. / sum);
}
void PolymerShape::_sample() {
    for (int i : ids)
        for (auto &g : spc.findMolecules(i))
            if (g.size() > 1) {
                Point r2 = vectorgyrationRadiusSquared(g);
                double rg2 = r2.sum();
                double re2 = spc.geo.sqdist(g.begin()->pos, (g.end() - 1)->pos);
                Rg[i] += sqrt(rg2);
                Re[i] += sqrt(re2);
                Rg2[i] += rg2;
                Re2[i] += re2; // end-2-end squared
                Rg2x[i] += r2.x();
                Rg2y[i] += r2.y();
                Rg2z[i] += r2.z();
                double rs = Re2[i].avg() / Rg2[i].avg(); // fluctuations in shape factor
                Rs[i] += rs;
                Rs2[i] += rs * rs;
            }
}
PolymerShape::PolymerShape(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "Polymer Shape";
    auto names = j.at("molecules").get<std::vector<std::string>>(); // molecule names
    ids = names2ids(molecules, names);                              // names --> molids
}
void AtomProfile::_from_json(const json &j) {
    ref = j.value("origo", Point(0, 0, 0));
    dir = j.value("dir", dir);
    file = j.at("file").get<std::string>();
    names = j.at("atoms").get<decltype(names)>();              // atom names
    auto vec_of_ids = names2ids(Faunus::atoms, names);         // names --> molids
    ids = std::set<int>(vec_of_ids.begin(), vec_of_ids.end()); // copy vector to set
    dr = j.value("dr", 0.1);
    tbl.setResolution(dr, 0);
    count_charge = j.value("charge", false);
}
void AtomProfile::_to_json(json &j) const {
    j = {{"origo", ref}, {"dir", dir}, {"atoms", names}, {"file", file}, {"dr", dr}, {"charge", count_charge}};
}
void AtomProfile::_sample() {
    for (auto &g : spc.groups)
        for (auto &p : g)
            if (ids.count(p.id) > 0) {
                Point rvec = spc.geo.vdist(p.pos, ref);
                double r = rvec.cwiseProduct(dir.cast<double>()).norm();
                if (count_charge)
                    tbl(r) += p.charge; // count charges
                else
                    tbl(r) += 1; // count atoms
            }
}
AtomProfile::AtomProfile(const json &j, Space &spc) : spc(spc) {
    name = "atomprofile";
    from_json(j);
}
void AtomProfile::_to_disk() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        tbl.stream_decorator = [&](std::ostream &o, double r, double N) {
            double Vr = 1;
            int dim = dir.sum();
            switch (dim) {
            case 3:
                Vr = 4 * pc::pi * std::pow(r, 2) * dr;
                break;
            case 2:
                Vr = 2 * pc::pi * r * dr;
                break;
            case 1:
                Vr = dr;
                break;
            default:
                throw std::runtime_error("bad dimension");
            }
            if (Vr < dr) // take care of the case where Vr=0
                Vr = dr; // then the volume element is simply dr

            N = N / double(cnt);                                          // average number of particles/charges
            o << r << " " << N << " " << N / Vr * 1e27 / pc::Nav << "\n"; // ... and molar concentration
        };
        f << "# r N rho/M\n" << tbl;
    }
}
void SlicedDensity::_from_json(const json &j) {
    file = j.at("file").get<std::string>();
    names = j.at("atoms").get<decltype(names)>(); // molecule names
    ids = names2ids(atoms, names);                // names --> molids
    dz = j.value("dz", 0.1);
    if (j.count("atomcom") == 1) {
        atomCOM = j.at("atomcom");
        idCOM = findName(atoms, atomCOM)->id();
    }
    N.setResolution(dz);
}
void SlicedDensity::_to_json(json &j) const {
    j = {{"atoms", names}, {"file", file}, {"dz", dz}, {"atomcom", atomCOM}};
}

void SlicedDensity::_sample() {
    Group<Particle> all(spc.p.begin(), spc.p.end());
    double zcm = 0;
    if (idCOM >= 0) { // calc. mass center of selected atoms
        auto slice = all.find_id(idCOM);
        zcm = Geometry::massCenter(slice.begin(), slice.end(), spc.geo.getBoundaryFunc()).z();
    }
    // count atoms in slices
    for (auto &g : spc.groups) // loop over all groups
        for (auto &i : g)      // loop over active particles
            if (std::find(ids.begin(), ids.end(), i.id) not_eq ids.end())
                N(i.pos.z() - zcm)++;
}
SlicedDensity::SlicedDensity(const json &j, Space &spc) : spc(spc) {
    name = "sliceddensity";
    from_json(j);
}
SlicedDensity::~SlicedDensity() {
    std::ofstream f(MPI::prefix + file);
    if (f and cnt > 0) {
        f << "# z rho/M\n";
        Point L = spc.geo.getLength();
        double halfz = 0.5 * L.z();
        double volume = L.x() * L.y() * dz;
        for (double z = -halfz; z <= halfz; z += dz)
            f << z << " " << N(z) / volume / cnt * 1e27 / pc::Nav << "\n";
    }
}
void ChargeFluctuations::_sample() {
    for (auto &g : spc.findMolecules(mol_iter->id(), Space::ACTIVE)) {
        size_t cnt = 0;
        for (auto &p : g) {
            idcnt[cnt][p.id]++;
            charge[cnt] += p.charge;
            cnt++;
        }
    }
}
void ChargeFluctuations::_to_json(json &j) const {
    std::vector<std::string> mainname; // main name of atom with fluctuating charge
    std::vector<double> qavg;          // average charge
    std::vector<double> qstdev;        // standard deviation of the charge
    for (size_t i = 0; i < idcnt.size(); ++i) {
        qavg.push_back(charge.at(i).avg());
        qstdev.push_back(charge.at(i).stdev());
        // we look for the id that was sampled most often
        auto id_max = std::max_element(std::begin(idcnt.at(i)), std::end(idcnt.at(i)),
                                       [](auto &p1, auto &p2) { return p1.second < p2.second; });
        mainname.push_back(atoms.at(id_max->first).name);
    }
    if (verbose)
        j = {{"dominant atoms", mainname}, {"<q>", qavg}, {"std", qstdev}};
    j["molecule"] = mol_iter->name;
    if (not file.empty())
        j["pqrfile"] = file;
}
void ChargeFluctuations::_to_disk() {
    if (not file.empty()) {
        auto molecules = spc.findMolecules(mol_iter->id(), Space::ALL);
        if (Faunus::size(molecules) > 0) {
            auto &g = *molecules.begin();
            ParticleVector pvec;        // temporary particle vector
            pvec.reserve(g.capacity()); // allocate required memory already now
            size_t cnt = 0;
            for (auto p = g.begin(); p < g.trueend(); ++p) {
                // we look for the id that was sampled most often
                auto id_max = std::max_element(std::begin(idcnt.at(cnt)), std::end(idcnt.at(cnt)),
                                               [](auto &p1, auto &p2) { return p1.second < p2.second; });
                pvec.push_back(Faunus::atoms.at(id_max->first));
                pvec.back().charge = charge.at(cnt).avg();
                pvec.back().pos = p->pos - g.cm;
                spc.geo.boundary(pvec.back().pos);
                cnt++;
            }
            FormatPQR::save(MPI::prefix + file, pvec, spc.geo.getLength());
        }
    }
}
ChargeFluctuations::ChargeFluctuations(const json &j, Space &spc) : spc(spc) {
    from_json(j);
    name = "chargefluctuations";
    file = j.value("pqrfile", ""s);
    verbose = j.value("verbose", true);
    auto molname = j.at("molecule").get<std::string>(); // molecule name
    mol_iter = findName(Faunus::molecules, molname);
    if (mol_iter == Faunus::molecules.end())
        throw std::runtime_error("unknown species '" + molname + "'");
    if (mol_iter->atomic)
        throw std::runtime_error("only molecular groups allowed");
    idcnt.resize(mol_iter->atoms.size());
    charge.resize(mol_iter->atoms.size());
}
void Multipole::_sample() {
    for (auto &g : spc.groups)
        if (!g.atomic) {
            auto &d = _map[g.id];
            Particle p = Faunus::toMultipole(g, spc.geo.getBoundaryFunc());
            d.Z += p.charge;
            d.mu += p.getExt().mulen;
            d.Z2 += p.charge * p.charge;
            d.mu2 += p.getExt().mulen * p.getExt().mulen;
        }
}
void Multipole::_to_json(json &j) const {
    json &k = j["molecules"];
    for (auto &d : _map)
        k[molecules[d.first].name] = {{"Z", d.second.Z.avg()},
                                      {"Z2", d.second.Z2.avg()},
                                      {"C", d.second.Z2.avg() - std::pow(d.second.Z.avg(), 2)},
                                      {u8::mu, d.second.mu.avg()},
                                      {u8::mu + u8::squared, d.second.mu2.avg()}};
}
Multipole::Multipole(const json &j, const Space &spc) : spc(spc) {
    from_json(j);
    name = "multipole";
}
void ScatteringFunction::_sample() {
    p.clear();
    for (int id : ids) { // loop over molecule names
        auto groups = spc.findMolecules(id);
        for (auto &g : groups) // loop over groups
            if (usecom && !g.atomic)
                p.push_back(g.cm);
            else
                for (auto &i : g) // loop over particle index in group
                    p.push_back(i.pos);
    }
    debye.sample(p, spc.geo.getVolume());
}
void ScatteringFunction::_to_json(json &j) const { j = {{"molecules", names}, {"com", usecom}}; }

ScatteringFunction::ScatteringFunction(const json &j, Space &spc) try : spc(spc), debye(j) {
    from_json(j);
    name = "scatter";
    usecom = j.value("com", true);
    filename = j.at("file").get<std::string>();
    names = j.at("molecules").get<decltype(names)>(); // molecule names
    ids = names2ids(molecules, names);                // names --> molids
} catch (std::exception &e) {
    throw std::runtime_error("debye formula: "s + e.what());
}

ScatteringFunction::~ScatteringFunction() { debye.save(filename); }

} // namespace Analysis
} // namespace Faunus
