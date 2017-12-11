#pragma once

#include "space.h"
#include "io.h"

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
                        _j["relative time"] = timer.result();
                        _j["nstep"] = steps;
                        _j["samples"] = cnt;
                    }
                    if (!cite.empty())
                        _j["citation"] = cite;
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

                    for (auto &i : Nmol)
                        rho_mol[i.first] += i.second/V;

                    for (auto &i : Natom)
                        rho_atom[i.first] += i.second/V;
                }

                inline void _to_json(json &j) const override {
                    auto &_j = j["atomic"];
                    for (auto &i : rho_atom)
                        _j[ atoms<Tparticle>.at(i.first).name ] = json({{ "c/M", i.second.avg() / 1.0_molar }});

                    auto &_jj = j["molecular"];
                    for (auto &i : rho_mol)
                        _jj[ molecules<Tpvec>.at(i.first).name ] = json({{ "c/M", i.second.avg() / 1.0_molar }});
                }

            public:
                    Density( const json &j, Tspace &spc ) : spc(spc) {
                        from_json(j);
                        name = "density";
                    }
        };

        class SystemEnergy : public Analysisbase {
            private:
                std::string file;
                std::ofstream f;
                std::function<double()> energyFunc;
                Average<double> uavg; //!< mean energy
                double uinit;

                inline void _sample() override {
                    double u = energyFunc();
                    uavg+=u;
                    f << cnt*steps << " " << u << "\n";
                }

                inline void _to_json(json &j) const override {
                    j["file"] = file;
                    j["init"] = uinit;
                    j["final"] = energyFunc();
                    if (cnt>0)
                        j["mean"] = uavg.avg();
                }

                inline void _from_json(const json &j) override {
                    file = j.at("file");
                    if (f)
                        f.close();
                    f.open(file);
                    if (!f)
                        throw std::runtime_error(name + ": cannot open output file " + file);
                }

            public:
                template<class Tenergy>
                    SystemEnergy( const json &j, Tenergy &pot ) {
                        from_json(j);
                        name = "systemenergy";
                        energyFunc = [&pot]() {
                            Change change;
                            change.all = true;
                            return pot.energy(change);
                        };
                        uinit = energyFunc(); // initial energy
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
                                if (f)
                                    f << std::setw(4) << json(spc);
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
                    cout << hist.getMap().size() << endl;
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
                    file = j.at("file");
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

                                    if ( it.key()=="aromrdf")
                                        push_back<AtomRDF<Tspace>>(it.value(), spc);

                                    if ( it.key()=="density")
                                        push_back<Density<Tspace>>(it.value(), spc);

                                    if ( it.key()=="molrdf")
                                        push_back<MoleculeRDF<Tspace>>(it.value(), spc);

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
