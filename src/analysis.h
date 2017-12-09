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
                    if (cnt>0) {
                        _j["relative time"] = timer.result();
                        _j["nstep"] = steps;
                        _j["samples"] = cnt;
                    }
                    if (!cite.empty())
                        _j["citation"] = cite;
                    _to_json(_j);
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

                void _sample() override {
                    writeFunc(file);
                }

        };

        /**
         * @brief Base class for distribution functions etc.
         */
        class PairFunctionBase : public Analysisbase {

            protected:
                struct data {
                    int dim;
                    double dr;
                    Table2D<double,double> hist;
                    Table2D<double,Average<double>> hist2;
                    std::string name1, name2, file, file2;
                    double Rhypersphere; // Radius of 2D hypersphere
                };
                std::vector<data> datavec;        // vector of data sets
                Average<double> V;                // average volume (angstrom^3)
                virtual void normalize(data &);
            private:
                virtual void update(data &d)=0;   // called on each defined data set
                void sample() override;
                json _json();

            public:
                PairFunctionBase(json, std::string);
                //virtual ~PairFunctionBase();
        };

        void PairFunctionBase::sample()
        {
            for (auto &d : datavec)
                update(d);
        }

        json PairFunctionBase::_json()
        {
            json j;
            auto &_j = j[name];
            for (auto &d : datavec)
                _j[ d.name1+"-"+d.name2 ] = {
                    { "dr", d.dr },
                    { "file", d.file },
                    { "file2", d.file2 },
                    { "dim", d.dim },
                    { "Rhyper", d.Rhypersphere }
                };
            return j;
        }

        PairFunctionBase::PairFunctionBase( json j, std::string name ) {
            try {
                for (auto &i : j.at("pairs"))
                    if (i.is_object())
                    {
                        data d;
                        d.file = i.at("file");
                        d.file2 = i.value("file2",d.file+".avg");
                        d.name1 = i.at("name1");
                        d.name2 = i.at("name2");
                        d.dim = i.value("dim", 3);
                        d.dr = i.value("dr", 0.1);
                        d.hist.setResolution(d.dr);
                        d.hist2.setResolution(d.dr);
                        d.Rhypersphere = i.value("Rhyper", -1.0);
                        datavec.push_back( d );
                    }
            }
            catch(std::exception& e) {
                throw std::runtime_error(name + ": " + e.what());
            }

            if (datavec.empty())
                std::cerr << name + ": no sample sets defined for analysis\n";
        }

        void PairFunctionBase::normalize(data &d)
        {
            assert(V.cnt>0);
            double Vr=1, sum = d.hist.sumy();
            for (auto &i : d.hist.getMap()) {
                if (d.dim==3)
                    Vr = 4 * pc::pi * pow(i.first,2) * d.dr;
                if (d.dim==2) {
                    Vr = 2 * pc::pi * i.first * d.dr;
                    if (d.Rhypersphere > 0)
                        Vr = 2.0*pc::pi*d.Rhypersphere*sin(i.first/d.Rhypersphere) * d.dr;
                }
                if (d.dim==1)
                    Vr = d.dr;
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
                    typedef typename Tspace::Tparticle Tparticle;

                    void update(data &d) override
                    {
                        V += spc.geo.getVolume( d.dim );
                        int N = spc.p.size();

                        auto it = findName( atoms<Tparticle>, d.name1 ).id;
                        if (it == atoms<Tparticle>.end())
                            throw std::runtime_error("unknown atom '" + d.name1 + "'");
                        int id1 = it->id();

                        it = findName( atoms<Tparticle>, d.name2 ).id;
                        if (it == atoms<Tparticle>.end())
                            throw std::runtime_error("unknown atom '" + d.name2 + "'");
                        int id2 = it->id();

                        for (int i=0; i<N-1; i++)
                            for (int j=i+1; j<N; j++)
                                if (
                                        ( spc.p[i].id==id1 && spc.p[j].id==id2 ) ||
                                        ( spc.p[i].id==id2 && spc.p[j].id==id1 )
                                   )
                                {
                                    double r = spc.geo.dist( spc.p[i], spc.p[j] );
                                    d.hist(r)++;
                                }
                    }

                public:
                    AtomRDF( json j, Tspace &spc ) : PairFunctionBase(j,
                            "Atomic Pair Distribution Function"), spc(spc) {}

                    ~AtomRDF()
                    {
                        for (auto &d : datavec) {
                            normalize(d);
                            d.hist.save( d.file );
                        }
                    }
            };

        /** @brief Same as `AtomRDF` but for molecules. Identical input. */
        template<class Tspace>
            class MoleculeRDF : public PairFunctionBase {
                private:
                    Tspace &spc;
                    void update(data &d) override
                    {
                        V += spc.geo.getVolume( d.dim );
                        auto g1 = spc.findMolecules( d.name1 );
                        auto g2 = spc.findMolecules( d.name2 );

                        if (d.name1!=d.name2)
                            for (auto i : g1)
                                for (auto j : g2) {
                                    double r = spc.geo.dist( i->cm, j->cm );
                                    d.hist(r)++;
                                }
                        else {
                            for (int i=0; i<(int)g1.size()-1; i++)
                                for (int j=i+1; j<(int)g1.size(); j++) {
                                    double r = spc.geo.dist( g1[i]->cm, g1[j]->cm );
                                    d.hist(r)++;
                                }
                        }
                    }

                public:
                    MoleculeRDF( json j, Tspace &spc ) : PairFunctionBase(j,
                            "Molecular Pair Distribution Function"), spc(spc) {}

                    ~MoleculeRDF()
                    {
                        for (auto &d : datavec) {
                            normalize(d);
                            d.hist.save( d.file );
                        }
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

                                    // additional analysis go here...
                                } catch (std::exception &e) {
                                    throw std::runtime_error("Error adding analysis '" + it.key() + "': " + e.what());
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
