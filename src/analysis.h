#pragma once

#include <numeric>
#include "space.h"
#include "io.h"
#include "energy.h"
#include "mpi.h"

namespace Faunus {

    namespace Analysis {

        class Analysisbase {
            private:
                inline virtual void _to_json(json &j) const {};
                inline virtual void _from_json(const json &j) {};
                virtual void _sample()=0;
                int stepcnt=0;
                TimeRelativeOfTotal<std::chrono::microseconds> timer;

            protected:
                int steps=0; //!< Sample interval (do not modify)
                int cnt=0;   //!< number of samples

            public:
                std::string name; //!< descriptive name
                std::string cite; //!< reference, url, doi etc. describing the analysis

                inline void to_json(json &j) const {
                    assert( !name.empty() );
                    auto &_j = j[name];
                    _to_json(_j);
                    if (cnt>0) {
                        _j["relative time"] = _round( timer.result() );
                        _j["nstep"] = steps;
                        _j["samples"] = cnt;
                    }
                    if (!cite.empty())
                        _j["reference"] = cite;
                } //!< JSON report w. statistics, output etc.

                inline void from_json(const json &j) {
                    steps = j.value("nstep", 0);
                    _from_json(j);
                } //!< configure from json object

                inline virtual void sample()
                {
                    stepcnt++;
                    if ( stepcnt == steps )
                    {
                        cnt++;
                        stepcnt = 0;
                        timer.start();
                        _sample();
                        timer.stop();
                    }
                }
        };

        inline void to_json(json &j, const Analysisbase &base) {
            base.to_json( j );
        }

        /**
         * @brief Excess chemical potential of molecules
         *
         * `widom`       | Description
         * :------------ | :-----------------------------------------
         * `molecule`    | Name of molecule to insert
         * `ninsert`     | Number of insertions per sample event
         * `dir=[1,1,1]` | Inserting directions
         * `absz=false`  | Apply `std::fabs` on all z-coordinates of inserted molecule
         *
         * This will insert a non-perturbing ghost molecule into
         * the system and calculate a Widom average to measure
         * the free energy of the insertion process, _i.e._ the
         * excess chemical potential:
         *
         * $$ \mu^{excess} = -k_BT\ln \langle e^{-\delta u/k_BT} \rangle_0 $$
         *
         * where $$delta u$$ is the energy change of the perturbation and the
         * average runs over the _unperturbed_ ensemble.
         *
         * The position and molecular rotation is random.
         * For use with rod-like particles on surfaces, the `absz`
         * keyword may be used to ensure orientations on only one
         * half-sphere.
         */
        template<typename Tspace>
            class WidomInsertion : public Analysisbase
        {
            private:
                typedef typename Tspace::Tparticle Tparticle;
                typedef typename Tspace::Tpvec Tpvec;
                typedef MoleculeData<Tpvec> TMoleculeData;

                Tspace& spc;
                Energy::Energybase* pot;
                RandomInserter<TMoleculeData> rins;
                std::string molname; // molecule name
                int ninsert;
                int molid;        // molecule id
                bool absolute_z=false;
                Average<double> expu;
                Change change;

                void _sample() override
                {
                    if (!change.empty()) {
                        Tpvec pin;
                        auto &g = spc.groups.at( change.groups.at(0).index );
                        assert(g.empty());
                        g.resize(g.capacity());
                        for ( int i = 0; i < ninsert; ++i )
                        {
                            pin = rins(spc.geo, spc.p, molecules<Tpvec>.at(molid));
                            if (!pin.empty()) {
                                if (absolute_z)
                                    for (auto &p : pin)
                                        p.pos.z() = std::fabs(p.pos.z());
                                double du = pot->energy(change);
                                expu += exp(-du); // widom average
                            }
                        }
                        g.resize(0);
                    }
                }

                void _to_json(json &j) const override
                {
                    using namespace u8;
                    double excess = -std::log(expu.avg());
                    j = {
                        { "dir", rins.dir }, { "molecule", molname },
                        { "insertions", expu.cnt }, { "absz", absolute_z },
                        { mu+"/kT",
                            {
                                { "excess", excess }
                            }
                        }
                    };
                }

                void _from_json(const json &j) override {
                    ninsert = j.at("ninsert");
                    molname = j.at("molecule");
                    absolute_z = j.value("absz", false);
                    rins.dir = j.value("dir", Point({1,1,1}) );

                    auto it = findName( molecules<Tpvec>, molname); // loop for molecule in topology
                    if (it!=molecules<Tpvec>.end()) {
                        molid = it->id();
                        auto m = spc.findMolecules(molid, Tspace::INACTIVE); // look for molecules in space
                        if (size(m)>0) { // did we find any?
                            if (m.begin()->size()==0) { // pick the first and check if it's really inactive
                                change.clear();
                                Change::data d; // construct change object
                                d.index = distance(spc.groups.begin(), m.begin()); // group index
                                d.all = true;
                                change.groups.push_back(d); // add to change object
                                return;
                            }
                        }
                    }
                    throw std::runtime_error(name+": no inactive '" + molname + "' groups found");
                }

            public:

                template<class Tenergy>
                    WidomInsertion( const json &j, Tspace &spc, Tenergy &pot ) : spc(spc), pot(&pot) {
                        from_json(j);
                        name = "widom";
                        cite = "doi:10/dkv4s6";
                        rins.checkOverlap = false;
                    }
        };

        template<class Tspace>
            class AtomProfile : public Analysisbase {
                private:
                    Tspace &spc;
                    typedef typename Tspace::Tparticle Tparticle;
                    Table2D<double,unsigned int> tbl;
                    std::vector<std::string> names;
                    std::vector<int> ids;
                    std::string file;
                    Point ref={0,0,0};
                    double dr=0.1;

                    void _from_json(const json &j) override {
                        ref = j.value("origo", Point(0,0,0));
                        file = j.at("file").get<std::string>();
                        names = j.at("atoms").get<decltype(names)>(); // molecule names
                        ids = names2ids(atoms<Tparticle>, names);     // names --> molids
                        dr = j.value("dr", 0.1);
                        tbl.setResolution(dr);
                    }

                    void _to_json(json &j) const override {
                        j = {{"origo", ref}, {"atoms", names}, {"file", file}, {"dr", dr}};
                    }

                    void _sample() override {
                        for (auto &g : spc.groups)
                            for (auto &p : g)
                                if (std::find(ids.begin(), ids.end(), p.id)!=ids.end()) {
                                    double r = spc.geo.vdist(p.pos, ref).norm();
                                    tbl(r)++;
                                }
                    }

                public:

                    AtomProfile(const json &j, Tspace &spc) : spc(spc) {
                        name = "atomprofile";
                        from_json(j);
                    }

                    ~AtomProfile() {
                        std::ofstream f(file);
                        if (f)
                            f << "# r N rho/M\n";
                            for (auto &m : tbl.getMap()) {
                                double r = m.first;
                                double N = m.second/double(cnt);
                                f << r << " " << N << " " << N/(4*pc::pi*r*r*dr)*1e27/pc::Nav << "\n";
                            }
                    }
            };

        /**
         * @brief Analysis of particle densities
         */
        template<class Tspace>
            class Density : public Analysisbase {
                private:
                    Tspace& spc;
                    typedef typename Tspace::Tparticle Tparticle;
                    typedef typename Tspace::Tpvec Tpvec;

                    std::map<int, Average<double>> rho_mol, rho_atom;
                    std::map<int,int> Nmol, Natom;
                    Average<double> Lavg, Vavg, invVavg;

                    inline void _sample() override {
                        // count atom and groups of individual id's
                        Nmol.clear();
                        Natom.clear();
                        for (auto &g : spc.groups)
                            if (g.atomic)
                                for (auto &p : g)
                                    Natom[p.id]++;
                            else
                                if (!g.empty())
                                    Nmol[g.id]++;

                        double V = spc.geo.getVolume();
                        Vavg += V;
                        Lavg += std::cbrt(V);
                        invVavg += 1/V;

                        for (auto &i : Nmol)
                            rho_mol[i.first] += i.second/V;

                        for (auto &i : Natom)
                            rho_atom[i.first] += i.second/V;
                    }

                    inline void _to_json(json &j) const override {
                        using namespace u8;
                        j[ bracket( "V" ) ] = Vavg.avg();
                        j[ bracket( "1/V" ) ] = invVavg.avg();
                        j[ bracket( cuberoot + "V" ) ] = Lavg.avg();
                        j[ cuberoot + bracket("V") ] = std::cbrt(Vavg.avg());

                        auto &_j = j["atomic"];
                        for (auto &i : rho_atom)
                            _j[ atoms<Tparticle>.at(i.first).name ] = json({{ "c/M", _round(i.second.avg() / 1.0_molar) }});

                        auto &_jj = j["molecular"];
                        for (auto &i : rho_mol)
                            _jj[ molecules<Tpvec>.at(i.first).name ] = json({{ "c/M", _round(i.second.avg() / 1.0_molar) }});
                        _roundjson(j,4);
                    }

                public:
                    Density( const json &j, Tspace &spc ) : spc(spc) {
                        from_json(j);
                        name = "density";
                    }
            };

        class SystemEnergy : public Analysisbase {
            private:
                std::string file, sep=" ";
                std::ofstream f;
                std::function<std::vector<double>()> energyFunc;
                Average<double> uavg; //!< mean energy
                std::vector<std::string> names;
                double uinit;

                inline void _sample() override {
                    auto ulist = energyFunc();
                    double tot = std::accumulate(ulist.begin(), ulist.end(), 0.0);
                    uavg+=tot;
                    f << cnt*steps << sep << tot;
                    for (auto u : ulist)
                        f << sep << u;
                    f << "\n";
                }

                inline void _to_json(json &j) const override {
                    j["file"] = file;
                    j["init"] = uinit;
                    j["final"] = energyFunc();
                    if (cnt>0)
                        j["mean"] = uavg.avg();
                    _roundjson(j,5);
                }

                inline void _from_json(const json &j) override {
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

            public:
                template<class Tspace>
                    SystemEnergy( const json &j, Energy::Hamiltonian<Tspace> &pot ) {
                        for (auto i : pot.vec)
                            names.push_back(i->name);
                        from_json(j);
                        name = "systemenergy";
                        energyFunc = [&pot]() {
                            Change change;
                            change.all = true;
                            std::vector<double> u;
                            u.reserve(pot.vec.size());
                            for (auto i : pot.vec)
                                u.push_back( i->energy(change) );
                            return u;
                        };
                        auto u = energyFunc();
                        uinit = std::accumulate(u.begin(), u.end(), 0.0); // initial energy
                    }
        }; //!< Save system energy to disk. Keywords: `nstep`, `file`.

        class SaveState : public Analysisbase {
            private:
                std::function<void(std::string)> writeFunc = nullptr;
                std::string file;
                inline void _to_json(json &j) const override {
                    j["file"] = file;
                }
                void _sample() override {
                    writeFunc(file);
                }
            public:
                template<class Tspace>
                    SaveState(const json &j, Tspace &spc) {
                        using std::ref;
                        using std::placeholders::_1;

                        from_json(j);
                        name = "savestate";
                        steps = j.value("nstep", -1);
                        file = j.at("file");
                        std::string suffix = file.substr(file.find_last_of(".") + 1);
                        file = MPI::prefix + file;
                        if ( suffix == "aam" )
                            writeFunc = std::bind(
                                    []( std::string file, Tspace &s ) { FormatAAM::save(file, s.p); },
                                    _1, std::ref(spc));
                        if ( suffix == "pqr" )
                            writeFunc = std::bind(
                                    []( std::string file, Tspace &s ) { FormatPQR::save(file, s.p, s.geo.getLength()); },
                                    _1, std::ref(spc));
                        if ( suffix == "state" )
                            writeFunc = [&spc](const std::string &file) {
                                std::ofstream f(file);
                                if (f) {
                                    f << std::setw(4) << json(spc);
                                }
                            };
                    }

                ~SaveState() {
                    if (steps==-1)
                        _sample();
                }
        };

        /**
         * @brief Base class for distribution functions etc.
         */
        class PairFunctionBase : public Analysisbase {
            protected:
                int dim;
                int id1=-1, id2=-1; // particle id (mol or atom)
                double dr;
                Table2D<double,double> hist;
                Table2D<double,Average<double>> hist2;
                std::string name1, name2, file, file2;
                double Rhypersphere; // Radius of 2D hypersphere
                Average<double> V;   // average volume (angstrom^3)
                virtual void normalize();
            private:

                inline void _from_json(const json &j) override {
                    file = j.at("file");
                    file2 = j.value("file2", file+".avg");
                    name1 = j.at("name1");
                    name2 = j.at("name2");
                    dim = j.value("dim", 3);
                    dr = j.value("dr", 0.1);
                    hist.setResolution(dr);
                    hist2.setResolution(dr);
                    Rhypersphere = j.value("Rhyper", -1.0);
                }

                inline void _to_json(json &j) const override {
                    j = {
                        { "dr", dr },
                        { "name1", name1 },
                        { "name2", name2 },
                        { "file", file },
                        { "file2", file2 },
                        { "dim", dim },
                        { "Rhyper", Rhypersphere }
                    };
                }

            public:
                inline PairFunctionBase(const json &j) {
                    from_json(j);
                }

                inline virtual ~PairFunctionBase()
                {
                    normalize();
                    hist.save( file );
                }
        };

        inline void PairFunctionBase::normalize()
        {
            //assert(V.cnt>0);
            double Vr=1, sum = hist.sumy();
            for (auto &i : hist.getMap()) {
                if (dim==3)
                    Vr = 4 * pc::pi * std::pow(i.first,2) * dr;
                if (dim==2) {
                    Vr = 2 * pc::pi * i.first * dr;
                    if ( Rhypersphere > 0)
                        Vr = 2.0*pc::pi*Rhypersphere*std::sin(i.first/Rhypersphere) * dr;
                }
                if (dim==1)
                    Vr = dr;
                i.second = i.second/sum * V/Vr;
            }
        }

        /**
         * @brief Atomic radial distribution function, g(r)
         *
         * We sample the pair correlation function between atom id's _i_ and _j_,
         * \f[
         * g_{ij}(r) = \frac{ N_{ij}(r) }{ \sum_{r=0}^{\infty} N_{ij}(r) } \cdot \frac{ \langle V \rangle }{ V(r) }
         * \f]
         *
         * where \f$ N_{ij}(r) \f$ is the number of observed pairs, accumulated over the
         * entire ensemble,  in the separation
         * interval \f$[r, r+dr] \f$ and \f$ V(r) \f$ is the corresponding volume element
         * which depends on dimensionality:
         *
         * \f$ V(r) \f$               | Dimensions (`dim`)
         * :------------------------- | :----------------------------------------
         * \f$ 4\pi r^2 dr \f$        | 3 (for particles in free space, default)
         * \f$ 2\pi r dr \f$          | 2 (for particles confined on a plane)
         * \f$ 2\pi R sin(r/R) dr \f$ | 2 (for particles confined on a 2D hypersphere surface, also needs input `Rhypersphere`)
         * \f$ dr \f$                 | 1 (for particles confined on a line)
         *
         * Example JSON input:
         *
         *     { "nstep":20, "pairs":
         *        [
         *          { "name1":"Na", "name2":"Cl", "dim":3, "dr":0.1, "file":"rdf-nacl.dat"},
         *          { "name1":"Na", "name2":"Na", "dim":3, "dr":0.1, "file":"rdf-nana.dat"}
         *        ]
         *     }
         */
        template<class Tspace>
            class AtomRDF : public PairFunctionBase {
                private:
                    Tspace &spc;

                    void _sample() override
                    {
                        V += spc.geo.getVolume( dim );
                        for ( auto i = spc.p.begin(); i != spc.p.end(); ++i )
                            for ( auto j=i; ++j != spc.p.end(); )
                                if (
                                        ( i->id==id1 && j->id==id2 ) ||
                                        ( i->id==id2 && j->id==id1 )
                                   )
                                {
                                    double r = std::sqrt( spc.geo.sqdist( i->pos, j->pos ) );
                                    hist(r)++;
                                }
                    }

                public:
                    AtomRDF( const json &j, Tspace &spc ) : PairFunctionBase(j), spc(spc) {
                        typedef typename Tspace::Tparticle Tparticle;
                        name = "atomrdf";

                        auto it = findName( atoms<Tparticle>, name1 );
                        if ( it == atoms<Tparticle>.end() )
                            throw std::runtime_error("unknown atom '" + name1 + "'");
                        id1 = it->id();

                        it = findName( atoms<Tparticle>, name2 );
                        if ( it == atoms<Tparticle>.end() )
                            throw std::runtime_error("unknown atom '" + name2 + "'");
                        id2 = it->id();
                    }
            };

        /** @brief Same as `AtomRDF` but for molecules. Identical input. */
        template<class Tspace>
            class MoleculeRDF : public PairFunctionBase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace &spc;

                    void _sample() override
                    {
                        V += spc.geo.getVolume( dim );
                        for ( auto i = spc.groups.begin(); i != spc.groups.end(); ++i )
                            for ( auto j=i; ++j != spc.groups.end(); )
                                if (
                                        ( i->id==id1 && j->id==id2 ) ||
                                        ( i->id==id2 && j->id==id1 )
                                   )
                                {
                                    double r = std::sqrt( spc.geo.sqdist( i->cm, j->cm ) );
                                    hist(r)++;
                                }
                    }

                public:
                    MoleculeRDF( const json &j, Tspace &spc ) : PairFunctionBase(j), spc(spc) {
                        name = "molrdf";

                        auto it = findName( molecules<Tpvec>, name1 );
                        if (it == molecules<Tpvec>.end())
                            throw std::runtime_error(name + ": unknown molecule '" + name1 + "'\n");
                        id1 = it->id();

                        it = findName( molecules<Tpvec>, name2 );
                        if (it == molecules<Tpvec>.end())
                            throw std::runtime_error(name + ": unknown molecule '" + name2 + "'\n");
                        id2 = it->id();
                        assert(id1>=0 && id2>=0);
                    }
            };

        /**
         * @brief Write XTC trajectory file
         *
         * JSON keywords: `nstep`, `file`
         */
        template<class Tspace>
            class XTCtraj : public Analysisbase
        {
            private:
                typedef typename Tspace::Tparticle Tparticle;
                void _to_json(json &j) const override {
                    j["file"] = file;
                }
                void _from_json(const json &j) override {
                    file = MPI::prefix + j.at("file").get<std::string>();
                }

                FormatXTC xtc;
                Tspace &spc;
                std::string file;

                void _sample() override
                {
                    xtc.setbox( spc.geo.getLength() );
                    bool rc = xtc.save(file, spc.p, Group<Tparticle>(spc.p.begin(), spc.p.end()));
                    if (rc==false)
                        std::cerr << "error saving xtc\n";
                }

            public:

                XTCtraj( const json &j, Tspace &s ) : xtc(1e6), spc(s)
            {
                from_json(j);
                name = "xtcfile";
                cite = "http://bit.ly/2A8lzpa";
            }
        };

        class VirtualVolume : public Analysisbase {
            private:
                double dV;
                Change c;
                Energy::Energybase& pot;
                std::function<double()> getVolume;
                std::function<void(double)> scaleVolume;
                Average<double> duexp; // < exp(-du/kT) >

                void _sample() override {
                    if (dV>0) {
                        double Vold = getVolume(),
                               Uold = pot.energy(c);
                        scaleVolume(Vold + dV);
                        double Unew = pot.energy(c);
                        scaleVolume(Vold);
                        duexp += std::exp(-(Uold - Unew));
                        assert(std::fabs(Uold - pot.energy(c)) < 1e-7);
                    }
                }

                void _from_json(const json &j) override { dV = j.at("dV"); }

                void _to_json(json &j) const override {
                    double pex = -std::log(duexp.avg()) / dV;
                    j["dV"] = dV;
                    j["Pex/mM"] = pex / 1.0_mM;
                    j["Pex/Pa"] = pex / 1.0_Pa;
                    j["Pex/kT/"+u8::angstrom+u8::cubed] = pex;
                    _roundjson(j,5);
                }

                void init(const json &j) {
                    from_json(j);
                    c.dV=true;
                    c.all=true;
                    name = "virtualvolume";
                    cite = "doi:10.1063/1.472721";
                }

            public:
                template<class Tspace>
                    VirtualVolume(const json &j, Tspace &spc, Energy::Energybase &pot) : pot(pot) {
                        init(j);
                        getVolume = [&spc](){ return spc.geo.getVolume(); };
                        scaleVolume = [&spc](double Vnew) { spc.scaleVolume(Vnew); };
                    }
        }; //!< Calculate excess pressure using virtual volume move

        class CombinedAnalysis : public BasePointerVector<Analysisbase> {
            public:
                template<class Tspace, class Tenergy>
                    CombinedAnalysis(const json &j, Tspace &spc, Tenergy &pot) {
                        for (auto &m : j.at("analysis")) // loop over move list
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                try {

                                    if (it.key()=="systemenergy")
                                        push_back<SystemEnergy>(it.value(), pot);

                                    if (it.key()=="savestate")
                                        push_back<SaveState>(it.value(), spc);

                                    if ( it.key()=="xtcfile")
                                        push_back<XTCtraj<Tspace>>(it.value(), spc);

                                    if ( it.key()=="atomprofile")
                                        push_back<AtomProfile<Tspace>>(it.value(), spc);

                                    if ( it.key()=="atomrdf")
                                        push_back<AtomRDF<Tspace>>(it.value(), spc);

                                    if ( it.key()=="density")
                                        push_back<Density<Tspace>>(it.value(), spc);

                                    if ( it.key()=="molrdf")
                                        push_back<MoleculeRDF<Tspace>>(it.value(), spc);

                                    if ( it.key()=="virtualvolume")
                                        push_back<VirtualVolume>(it.value(), spc, pot);

                                    if ( it.key()=="widom")
                                        push_back<WidomInsertion<Tspace>>(it.value(), spc, pot);

                                    // additional analysis go here...
                                } catch (std::exception &e) {
                                    throw std::runtime_error("Error adding analysis '"
                                            + it.key() + "': " + e.what() + "\n");
                                }
                            }
                    }

                inline void sample() {
                    for (auto i : this->vec)
                        i->sample();
                }

        }; //!< Aggregates analysis

        /** @brief Example analysis */
        template<class T, class Enable = void>
            struct _analyse {
                void sample(T &p) {
                    std::cout << "not a dipole!" << std::endl;
                } //!< Sample
            }; // primary template

        /** @brief Example analysis */
        template<class T>
            struct _analyse<T, typename std::enable_if<std::is_base_of<Dipole, T>::value>::type> {
                void sample(T &p) {
                    std::cout << "dipole!" << std::endl;
                } //!< Sample
            }; // specialized template

    }//namespace

}//namespace
