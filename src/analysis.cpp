#include "analysis.h"

void Faunus::Analysis::to_json(Faunus::json &j, const Faunus::Analysis::Analysisbase &base) {
    base.to_json( j );
}

Faunus::Analysis::Analysisbase::~Analysisbase() {}

/*
 * This is always called by the destructor, but may be called
 * any number of times earlier.
 */
void Faunus::Analysis::Analysisbase::_to_disk() {}

void Faunus::Analysis::Analysisbase::to_disk() {
    _to_disk();
}

void Faunus::Analysis::Analysisbase::sample() {
    totstepcnt++;
    stepcnt++;
    if ( stepcnt == steps ) {
        stepcnt = 0;
        if ( totstepcnt > nskip ) {
            cnt++;
            timer.start();
            _sample();
            timer.stop();
        }
    }
}

void Faunus::Analysis::Analysisbase::from_json(const Faunus::json &j) {
    steps = j.value("nstep", 0);
    nskip = j.value("nskip", 0);
    _from_json(j);
}

void Faunus::Analysis::Analysisbase::to_json(Faunus::json &j) const {
    assert( not name.empty() );
    auto &_j = j[name];
    _to_json(_j);
    if (cnt>0) {
        if (timer.result() > 0.01) // only print if more than 1% of the time
            _j["relative time"] = _round( timer.result() );
        _j["nstep"] = steps;
        _j["samples"] = cnt;
        if (nskip>0)
            _j["nskip"] = nskip;
    }
    if (not cite.empty())
        _j["reference"] = cite;
}

void Faunus::Analysis::Analysisbase::_to_json(Faunus::json&) const {}

void Faunus::Analysis::Analysisbase::_from_json(const Faunus::json&) {}

void Faunus::Analysis::SystemEnergy::normalize() {
    //assert(V.cnt>0);
    double sum = ehist.sumy();
    for (auto &i : ehist.getMap()) {
        i.second = i.second/sum ;
    }
}

void Faunus::Analysis::SystemEnergy::_sample() {
    auto ulist = energyFunc();
    double tot = std::accumulate(ulist.begin(), ulist.end(), 0.0);
    if (not std::isinf(tot)) {
        uavg+=tot;
        u2avg+=tot*tot;
    }
    f << cnt*steps << sep << tot;
    for (auto u : ulist)
        f << sep << u;
    f << "\n";
    //ehist(tot)++;
}

void Faunus::Analysis::SystemEnergy::_to_json(Faunus::json &j) const {
    j = { {"file", file}, {"init",uinit}, {"final", energyFunc()} };
    if (cnt>0) {
        j["mean"] = uavg.avg();
        j["Cv/kT"+u8::squared] = (u2avg.avg() - std::pow(uavg.avg(),2)) / pc::temperature;
    }
    _roundjson(j,5);
    //normalize();
    //ehist.save( "distofstates.dat" );

}

void Faunus::Analysis::SystemEnergy::_from_json(const Faunus::json &j) {
    file = MPI::prefix + j.at("file").get<std::string>();
    if (f)
        f.close();
    f.open(file);
    if (!f)
        throw std::runtime_error(name + ": cannot open output file " + file);
    assert(!names.empty());
    std::string suffix = file.substr(file.find_last_of(".") + 1);
    if (suffix=="csv")
        sep=",";
    else {
        sep=" ";
        f << "#";
    }
    f << "total";
    for (auto &n : names)
        f << sep << n;
    f << "\n";
}

void Faunus::Analysis::SaveState::_to_json(Faunus::json &j) const {
    j["file"] = file;
}

void Faunus::Analysis::SaveState::_sample() {
    writeFunc(file);
}

Faunus::Analysis::SaveState::~SaveState() {
    if (steps==-1)
        _sample();
}

Faunus::Analysis::PairFunctionBase::PairFunctionBase(const Faunus::json &j) { from_json(j); }

Faunus::Analysis::PairFunctionBase::~PairFunctionBase() {
    std::ofstream f(MPI::prefix + file);
    if (f) {
        double Vr=1, sum = hist.sumy();
        hist.stream_decorator = [&](std::ostream &o, double r, double N) {
            if (dim==3)
                Vr = 4 * pc::pi * std::pow(r,2) * dr;
            else if (dim==2) {
                Vr = 2 * pc::pi * r * dr;
                if ( Rhypersphere > 0)
                    Vr = 2.0*pc::pi*Rhypersphere*std::sin(r/Rhypersphere) * dr;
            }
            else if (dim==1)
                Vr = dr;
            if (Vr>0)
                o << r << " " << N*V/(Vr*sum) << "\n";
        };
        f << hist;
    }
}

void Faunus::Analysis::PairFunctionBase::_to_json(Faunus::json &j) const {
    j = {
        {"dr", dr/1.0_angstrom},
        {"name1", name1},
        {"name2", name2},
        {"file", file},
        {"dim", dim}
    };
    if (Rhypersphere>0)
        j["Rhyper"] = Rhypersphere;
}

void Faunus::Analysis::PairFunctionBase::_from_json(const Faunus::json &j) {
    assertKeys(j, {
            "file", "name1", "name2", "dim", "dr", "Rhyper", "nstep", "nskip"
            });
    file = j.at("file");
    name1 = j.at("name1");
    name2 = j.at("name2");
    dim = j.value("dim", 3);
    dr = j.value("dr", 0.1) * 1.0_angstrom;
    hist.setResolution(dr, 0);
    Rhypersphere = j.value("Rhyper", -1.0);
}

void Faunus::Analysis::VirtualVolume::_sample() {
    if (fabs(dV)>1e-10) {
        double Vold = getVolume(), Uold = pot.energy(c);
        scaleVolume(Vold + dV);
        double Unew = pot.energy(c);
        scaleVolume(Vold);

        // check if energy change is too big for exp()
        double x=pc::infty, du = Unew-Uold;
        if (-du < pc::max_exp_argument)
            x = std::exp(-du);
        if (std::isinf(x)) {
            std::cerr << name+": skipping sample event due to excessive energy likely due to overlaps.";
            cnt--; // cnt is incremented by sample() so we need to decrease
        } else {
            assert(not std::isnan(x));
            duexp += x;
#ifndef NDEBUG
            // check volume and particle positions are properly restored
            double err = std::fabs((Uold-pot.energy(c))/Uold); // must be ~zero!
            if (std::isfinite(err)) // catch if Uold==0
                assert(err < 1e-4);
#endif
        }
    }
}

void Faunus::Analysis::VirtualVolume::_from_json(const Faunus::json &j) { dV = j.at("dV"); }

void Faunus::Analysis::VirtualVolume::_to_json(Faunus::json &j) const {
    double pex = log(duexp.avg()) / dV; // excess pressure
    j = {
        {"dV", dV}, {"Pex/mM", pex/1.0_mM},
        {"Pex/Pa", pex/1.0_Pa}, {"Pex/kT/"+u8::angstrom+u8::cubed, pex}
    };
    _roundjson(j,5);
}

void Faunus::Analysis::QRtraj::_sample() { write_to_file(); }

void Faunus::Analysis::QRtraj::_to_json(Faunus::json &j) const { j = {{"file", file}}; }

void Faunus::Analysis::CombinedAnalysis::sample() {
    for (auto i : this->vec) i->sample();
}

Faunus::Analysis::CombinedAnalysis::~CombinedAnalysis() {
    for (auto i : this->vec) i->to_disk();
}
