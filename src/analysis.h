#pragma once

#include <numeric>
#include <cerrno>
#include "move.h"
#include "space.h"
#include "io.h"
#include "energy.h"
#include "mpi.h"
#include "scatter.h"

namespace Faunus {

    namespace Analysis {

        class Analysisbase {
            virtual void _to_json(json &j) const;
            virtual void _from_json(const json &j);
            virtual void _sample()=0;
            int stepcnt=0;
            int totstepcnt=0;
            TimeRelativeOfTotal<std::chrono::microseconds> timer;

            protected:
            int steps=0; //!< Sample interval (do not modify)
            int nskip=0; //!< MC steps to skip before sampling
            int cnt=0;   //!< number of samples

            public:
            std::string name; //!< descriptive name
            std::string cite; //!< reference, url, doi etc. describing the analysis

            void to_json(json &j) const; //!< JSON report w. statistics, output etc.
            void from_json(const json &j); //!< configure from json object
            virtual void sample();
            virtual ~Analysisbase();
        };

        void to_json(json &j, const Analysisbase &base);

        /*
         * @brief Sample and save reaction coordinates to a file
         */
        class FileReactionCoordinate : public Analysisbase {
            private:
                Average<double> avg;
                std::string type, filename;
                std::ofstream file;
                std::shared_ptr<ReactionCoordinate::ReactionCoordinateBase> rc=nullptr;

                inline void _to_json(json &j) const override {
                    j = *rc;
                    j["type"] = type;
                    j["file"] = filename;
                    j.erase("range");     // these are for penalty function
                    j.erase("resolution");// use only, so no need to show
                    if (cnt>0)
                        j["average"] = avg.avg();
                }

                inline void _sample() override {
                    if (file) {
                        double val = (*rc)();
                        avg += val;
                        file << cnt*steps << " " << val << " " << avg.avg() << "\n";
                    }
                }

            public:
                template<class Tspace>
                    FileReactionCoordinate(const json &j, Tspace &spc) {
                        using namespace Faunus::ReactionCoordinate;
                        from_json(j);
                        name = "reactioncoordinate";
                        filename = MPI::prefix + j.at("file").get<std::string>();
                        file.open(filename); // output file
                        type = j.at("type").get<std::string>();
                        try {
                            if      (type=="atom")     rc = std::make_shared<AtomProperty>(j, spc);
                            else if (type=="molecule") rc = std::make_shared<MoleculeProperty>(j, spc);
                            else if (type=="system")   rc = std::make_shared<SystemProperty>(j, spc);
                            else if (type=="cmcm")     rc = std::make_shared<MassCenterSeparation>(j, spc);
                            if (rc==nullptr)
                                throw std::runtime_error("unknown coordinate type");

                        } catch (std::exception &e) {
                            throw std::runtime_error("error for reaction coordinate '"
                                    + type + "': " + e.what() + usageTip["coords=["+type+"]"]  );
                        }
                    }
        };

        /**
         * @brief Excess chemical potential of molecules
         */
        template<typename Tspace>
            class WidomInsertion : public Analysisbase {
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

                void _sample() override {
                    if (!change.empty()) {
                        Tpvec pin;
                        auto &g = spc.groups.at( change.groups.at(0).index );
                        assert(g.empty());
                        g.resize(g.capacity()); // active group
                        for ( int i = 0; i < ninsert; ++i )
                        {
                            pin = rins(spc.geo, spc.p, molecules<Tpvec>.at(molid));
                            if (!pin.empty()) {
                                if (absolute_z)
                                    for (auto &p : pin)
                                        p.pos.z() = std::fabs(p.pos.z());

                                assert(pin.size() == g.size());
                                spc.geo.randompos( pin[0].pos, random );
                                spc.geo.randompos( pin[1].pos, random );

                                std::copy(pin.begin(), pin.end(), g.begin()); // copy into ghost group
                                if (!g.atomic) // update molecular mass-center
                                    g.cm = Geometry::massCenter(g.begin(), g.end(),
                                            spc.geo.getBoundaryFunc(), -g.begin()->pos);

                                expu += exp( -pot->energy(change) ); // widom average
                            }
                        }
                        g.resize(0); // deactive molecule
                    }
                }

                void _to_json(json &j) const override {
                    double excess = -std::log(expu.avg());
                    j = {
                        { "dir", rins.dir }, { "molecule", molname },
                        { "insertions", expu.cnt }, { "absz", absolute_z },
                        { u8::mu+"/kT",
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
                                d.internal = m.begin()->atomic;
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
                    }
            };

            /**
             * @brief Single particle hard sphere Widom insertion with charge scaling
             *
             * This will calculate excess chemical potentials for single particles
             * in the primitive model of electrolytes. Use the `add()` functions
             * to add test or *ghost* particles and call `sample()` to perform single
             * particle insertions.
             * Inserted particles are *non-perturbing* and thus removed again after
             * sampling. Electroneutrality for insertions of charged species is
             * maintaing by charge re-scaling according to
             *
             * - [Svensson and Woodward, Mol. Phys. 1988, 64:247]
             *   (http://doi.org/ft9bv9)
             *
             * Currently this works **only** for the primitive model of electrolytes, i.e.
             * hard, charged spheres interacting with a Coulomb potential.
             *
             * JSON input:
             *
             *  Keyword    | Description
             *  ---------- | ---------------------------------------
             *  `lB`       | Bjerrum length (angstrom)
             *  `ninsert`  | Number of intertions per sampling event
             *  `nstep`    | Sample every n'th step
             *
             * @warning Works only for the primitive model
             * @note This is a conversion of the Widom routine found in the `bulk.f`
             *       fortran program by Bolhuis/Jonsson/Akesson at Lund University.
             * @author Martin Trulsson and Mikael Lund
             * @date Lund / Prague 2007-2008.
             */
            // template<class Tspace>
            //     class WidomScaled : public AnalysisBase
            // {
            //
            //     private:
            //
                    // typedef typename Tspace::Tparticle Tparticle;
                    // typedef typename Tspace::Tpvec Tpvec;
                    // typedef MoleculeData<Tpvec> TMoleculeData;
            //
            //         Tspace &spc;
            //         Energy::Energybase* pot;
            //         RandomInserter<TMoleculeData> rins;
            //
            //         Tvec chel;          //!< electrostatic
            //         Tvec chhc;          //!< hard collision
            //         Tvec chex;          //!< excess
            //         Tvec chexw;         //!< excess
            //         Tvec chtot;         //!< total
            //         vector <Tvec> ewden; //!< charging denominator
            //         vector <Tvec> ewnom; //!< charging nominator
            //         vector <Tvec> chint; //!< charging integrand
            //         Tvec chid;          //!< ideal term
            //         Tvec expuw;
            //         vector<int> ihc, irej;
            //         int ghostin;        //< ghost insertions
            //         double lB;          //!< Bjerrum length
            //
            //         void init()
            //         {
            //             int gspec = g.size();
            //             chel.resize(gspec);
            //             chhc.resize(gspec);
            //             chex.resize(gspec);
            //             chtot.resize(gspec);
            //             ewden.resize(gspec);
            //             ewnom.resize(gspec);
            //             chint.resize(gspec);
            //             expuw.resize(gspec);
            //             chexw.resize(gspec);
            //             ihc.resize(gspec);
            //             irej.resize(gspec);
            //
            //             for ( int i = 0; i < gspec; i++ )
            //             {
            //                 chel[i] = 0;
            //                 chhc[i] = 0;
            //                 chex[i] = 0;
            //                 chtot[i] = 0;
            //                 ihc[i] = 0;
            //                 ewden[i].resize(11);
            //                 ewnom[i].resize(11);
            //                 chint[i].resize(11);
            //                 for ( int j = 0; j < 11; j++ )
            //                     ewden[i][j] = ewnom[i][j] = chint[i][j] = 0;
            //             }
            //         }
            //
            //         template<class Tgeo>
            //             bool overlap( const Tparticle &a, const Tparticle &b, const Tgeo &geo )
            //             {
            //                 double s = a.radius + b.radius;
            //                 return (geo.sqdist(a, b) < s * s) ? true : false;
            //             }
            //
            //         string _info() override
            //         {
            //             using namespace textio;
            //             std::ostringstream o;
            //             double aint4, aint2, aint1;
            //             for ( size_t i = 0; i < g.size(); i++ )
            //             {
            //                 for ( int cint = 0; cint < 11; cint++ )
            //                 {
            //                     if ( ewden[i][cint] == 0 )
            //                         std::cerr << "# WARNING: Widom denominator equals zero" << endl;
            //                     else
            //                         chint[i][cint] = ewnom[i][cint] / ewden[i][cint];
            //                 }
            //                 aint4 = chint[i][1] + chint[i][3] + chint[i][5] + chint[i][7] + chint[i][9];
            //                 aint2 = chint[i][2] + chint[i][4] + chint[i][6] + chint[i][8];
            //                 aint1 = chint[i][0] + chint[i][10];
            //                 chel[i] = 1. / 30. * (aint1 + 2 * aint2 + 4 * aint4);
            //             }
            //
            //             int cnttot = cnt * ghostin;
            //             o << pad(SUB, w, "Number of insertions") << cnttot << endl
            //                 << pad(SUB, w, "Excess chemical potentials (kT)") << endl
            //                 << "             total    elec.   hs             z        r" << endl;
            //             char w = 10;
            //             for ( size_t i = 0; i < g.size(); i++ )
            //             {
            //                 chhc[i] = -log(double(cnttot - ihc[i]) / cnttot);
            //                 chexw[i] = -log(expuw[i]);
            //                 chex[i] = chhc[i] + chel[i];
            //                 o.unsetf(std::ios_base::floatfield);
            //                 o << "    [" << i << "] "
            //                     << std::setprecision(4)
            //                     << std::setw(w) << chex[i]
            //                     << std::setw(w) << chel[i]
            //                     << std::setw(w) << chhc[i]
            //                     << std::setprecision(2) << std::fixed
            //                     << std::setw(w) << g[i].charge
            //                     << std::setw(w) << g[i].radius << endl;
            //             }
            //             return o.str();
            //         }
            //
            //         Tmjson _json() override
            //         {
            //             if ( cnt * ghostin>0 )
            //                 return {
            //                     { name,
            //                         {
            //                             { "number of insertions", cnt*ghostin },
            //                             { "ninsert", ghostin },
            //                             { "lB", lB }
            //                         }
            //                     }
            //                 };
            //             return Tmjson();
            //         }
            //
            //
            //         void _sample() override
            //         {
            //             auto &geo = spc.geo;
            //             auto &p = spc.p;
            //             if ( !g.empty())
            //                 if ( !p.empty())
            //                 {
            //                     Tparticle ghost;
            //                     double u, cu;
            //                     for ( int i = 0; i < ghostin; i++ )
            //                     {
            //                         geo.randompos(ghost);
            //                         int goverlap = 0;
            //                         for ( size_t k = 0; k < g.size(); k++ )
            //                         {
            //                             ghost.radius = g[k].radius;
            //                             irej[k] = 0;
            //                             int j = 0;
            //                             while ( !overlap(ghost, p[j], geo) && j < (int) p.size())
            //                                 j++;
            //                             if ( j != (int) p.size())
            //                             {
            //                                 ihc[k]++;
            //                                 irej[k] = 1;
            //                                 goverlap++;
            //                             }
            //                         }
            //
            //                         if ( goverlap != (int) g.size())
            //                         {
            //                             cu = 0;
            //                             u = 0;  //elelectric potential (Coulomb only!)
            //                             for ( auto &i : p )
            //                             {
            //                                 double invdi = 1 / geo.dist(ghost, i);
            //                                 cu += invdi;
            //                                 u += invdi * i.charge;
            //                             }
            //                             cu = cu * lB;
            //                             u = u * lB;
            //                             double ew, ewla, ewd;
            //                             for ( size_t k = 0; k < g.size(); k++ )
            //                             {
            //                                 if ( irej[k] == 0 )
            //                                 {
            //                                     expuw[k] += exp(-u * g[k].charge);
            //                                     for ( int cint = 0; cint < 11; cint++ )
            //                                     {
            //                                         ew = g[k].charge *
            //                                             (u - double(cint) * 0.1 * g[k].charge * cu / double(p.size()));
            //                                         ewla = ew * double(cint) * 0.1;
            //                                         ewd = exp(-ewla);
            //                                         ewden[k][cint] += ewd;
            //                                         ewnom[k][cint] += ew * ewd;
            //                                     }
            //                                 }
            //                             }
            //                         }
            //                     }
            //                 }
            //         }
            //
            //     public:
            //
            //         WidomScaled( Tmjson &j, Tspace &spc ) : spc(spc), AnalysisBase(j)
            //     {
            //         lB = j.value("lB", 7.0);
            //         ghostin = j.value("ninsert", 10);
            //         name = "Single particle Widom insertion w. charge scaling";
            //         cite = "doi:10/ft9bv9 + doi:10/dkv4s6";
            //
            //         add(spc.p);
            //     }
            //
            //         /**
            //          * @brief Add ghost particle
            //          *
            //          * This will add particle `p` to the list of ghost particles
            //          * to insert.
            //          */
            //         void add( const Tparticle &p )
            //         {
            //             g.push_back(p);
            //             init();
            //         }
            //
            //         /**
            //          * @brief Add ghost particles
            //          *
            //          * This will scan the particle vector for particles and each unique type
            //          * will be added to the list a ghost particles to insert.
            //          */
            //         template<class Tpvec>
            //             void add( const Tpvec &p )
            //             {
            //                 std::set<typename Tparticle::Tid> ids;
            //                 for ( auto &i : p )
            //                     ids.insert(i.id);
            //                 for ( auto i : ids )
            //                 {
            //                     Tparticle a;
            //                     a = atom[i];
            //                     add(a);
            //                 }
            //             }
            //
            // }; // end of WidomScaled


        template<class Tspace>
            class AtomProfile : public Analysisbase {
                Tspace &spc;
                typedef typename Tspace::Tparticle Tparticle;
                Equidistant2DTable<double,double> tbl;
                std::vector<std::string> names; // atom names to analyse
                std::set<int> ids; // atom ids to analyse
                std::string file; // output filename
                Point ref={0,0,0};
                double dr; // radial resolution
                bool count_charge=false;
                bool Vnormalise=true;

                void _from_json(const json &j) override {
                    ref = j.value("origo", Point(0,0,0));
                    file = j.at("file").get<std::string>();
                    names = j.at("atoms").get<decltype(names)>(); // atom names
                    auto vec_of_ids = names2ids(Faunus::atoms, names);     // names --> molids
                    ids = std::set<int>(vec_of_ids.begin(), vec_of_ids.end()); // copy vector to set
                    dr = j.value("dr", 0.1);
                    tbl.setResolution(dr,0);
                    count_charge = j.value("charge", false);
                }

                void _to_json(json &j) const override {
                    j = {
                        {"origo", ref}, {"atoms", names}, {"file", file}, {"dr", dr},
                        {"charge", count_charge} };
                }

                void _sample() override {
                    for (auto &g : spc.groups)
                        for (auto &p : g)
                            if (ids.count(p.id)!=0) {
                                double r = spc.geo.vdist(p.pos, ref).norm();
                                if (count_charge)
                                    tbl(r) += p.charge; // count charge
                                else
                                    tbl(r) += 1; // count atoms
                            }
                }

                public:

                AtomProfile(const json &j, Tspace &spc) : spc(spc) {
                    name = "atomprofile";
                    from_json(j);
                }

                ~AtomProfile() {
                    std::ofstream f(MPI::prefix + file);
                    if (f) {
                        tbl.stream_decorator = [&](std::ostream &o, double r, double N) {
                            if (r>0) {
                                double V = 4*pc::pi*r*r*dr;
                                N = N/double(cnt);
                                o << r << " " << N << " " << N / V * 1e27 / pc::Nav << "\n";
                            }
                        };
                        f << "# r N rho/M\n" << tbl;
                    }
                }
            };

        /**
         * @brief Measures the density of atoms along z axis
         */
        template<class Tspace>
            class SlicedDensity : public Analysisbase {
                Tspace &spc;
                typedef typename Tspace::Tparticle Tparticle;
                Table2D<double,unsigned int> N; // N(z)
                std::vector<std::string> names;
                std::vector<int> ids;
                std::string file;
                double dz;

                void _from_json(const json &j) override {
                    file = j.at("file").get<std::string>();
                    names = j.at("atoms").get<decltype(names)>(); // molecule names
                    ids = names2ids(atoms, names);     // names --> molids
                    dz = j.value("dz", 0.1);
                    N.setResolution(dz);
                }

                void _to_json(json &j) const override {
                    j = {{"atoms", names}, {"file", file}, {"dz", dz}};
                }

                void _sample() override {
                    // count atoms in slices
                    for (auto &g : spc.groups) // loop over all groups
                        for (auto &i : g)      // loop over active particles
                            if (std::find(ids.begin(), ids.end(), i.id) not_eq ids.end())
                                N( i.pos.z() )++;
                }

                public:

                SlicedDensity(const json &j, Tspace &spc) : spc(spc) {
                    name = "sliceddensity";
                    from_json(j);
                }

                ~SlicedDensity() {
                    std::ofstream f(MPI::prefix + file);
                    if (f and cnt>0) {
                        f << "# z rho/M\n";
                        Point L = spc.geo.getLength();
                        double halfz = 0.5*L.z();
                        double volume = L.x() * L.y() * dz;
                        for (double z=-halfz; z<=halfz; z+=dz)
                            f << z << " " << N(z) / volume / cnt * 1e27 / pc::Nav << "\n";
                    }
                }
            };

        /**
         * @brief Analysis of particle densities
         */
        template<class Tspace>
            class Density : public Analysisbase {
                Tspace& spc;
                typedef typename Tspace::Tparticle Tparticle;
                typedef typename Tspace::Tpvec Tpvec;

                std::map<int, Table2D<double,double>> seldhist;  // Density histograms for selected atoms
                std::map<int, Table2D<double,double>> atmdhist;  // Density histograms for atomic molecules
                std::map<int, Table2D<double,double>> moldhist;  // Density histograms for molecules
                std::map<int, Average<double>> rho_mol, rho_atom;
                std::map<int,int> Nmol, Natom;
                Average<double> Lavg, Vavg, invVavg;

                int capacity_limit=10; // issue warning if capacity get lower than this

                void _sample() override {
                    // count atom and groups of individual id's
                    Nmol.clear();
                    Natom.clear();

                    // make sure all atom counts are initially zero
                    for (auto &g : spc.groups) {
                        if (g.atomic)
                            for ( auto p = g.begin(); p < g.trueend() ; ++p)
                                Natom[p->id] = 0;
                        else
                            Nmol[g.id]=0;
                    }

                    double V = spc.geo.getVolume();
                    Vavg += V;
                    Lavg += std::cbrt(V);
                    invVavg += 1/V;

                    for (auto &g : spc.groups)
                        if (g.atomic) {
                            for (auto &p : g)
                                Natom[p.id]++;
                            atmdhist[g.id]( g.size() )++;
                        }
                        else if (not g.empty())
                            Nmol[g.id]++;

                    for (auto &i : Nmol) {
                        rho_mol[i.first] += i.second/V;
                        moldhist[i.first](i.second)++;
                    }

                    for (auto &i : Natom)
                        rho_atom[i.first] += i.second/V;

                    if ( reactions<Tpvec>.size()>0 ) { // in case of reactions involving atoms (swap moves)
                        for (auto &rit : reactions<Tpvec> ) {
                            for (auto pid : rit._prodid_a) {
                                auto atomlist = spc.findAtoms(pid.first);
                                seldhist[ pid.first ]( size(atomlist) )++;
                            }
                            for (auto rid : rit._reacid_a) {
                                auto atomlist = spc.findAtoms(rid.first);
                                seldhist[ rid.first ]( size(atomlist) )++;
                            }
                        }
                    }
                }

                void _to_json(json &j) const override {
                    using namespace u8;
                    j[ bracket( "V" ) ] = Vavg.avg();
                    j[ bracket( "1/V" ) ] = invVavg.avg();
                    j[ bracket( cuberoot + "V" ) ] = Lavg.avg();
                    j[ cuberoot + bracket("V") ] = std::cbrt(Vavg.avg());

                    auto &_j = j["atomic"];
                    for (auto &i : rho_atom)
                        if (i.second.cnt>0)
                            _j[ atoms.at(i.first).name ] = json({{ "c/M", i.second.avg() / 1.0_molar }});

                    auto &_jj = j["molecular"];
                    for (auto &i : rho_mol)
                        if (i.second.cnt>0)
                            _jj[ molecules<Tpvec>.at(i.first).name ] = json({{ "c/M", i.second.avg() / 1.0_molar }});
                    _roundjson(j,4);
                }

                public:
                Density( const json &j, Tspace &spc ) : spc(spc) {
                    from_json(j);
                    name = "density";
                }
                virtual ~Density() {
                    normalize();
                    for ( auto &m: atmdhist)
                        m.second.save( "rho-"s + molecules<Tpvec>.at(m.first).name + ".dat" );
                    for ( auto &m: moldhist)
                        m.second.save( "rho-"s + molecules<Tpvec>.at(m.first).name + ".dat" );
                    if ( reactions<Tpvec>.size()>0 ) { // in case of reactions involving atoms (swap moves)
                        for (auto &rit : reactions<Tpvec> ) {
                            for (auto pid : rit._prodid_a) 
                                seldhist.at( pid.first ).save( "rho-"s + atoms.at( pid.first ).name + ".dat" );
                            for (auto rid : rit._reacid_a) 
                                seldhist.at( rid.first ).save( "rho-"s + atoms.at( rid.first ).name + ".dat" );
                        }
                    } 
                }
                void normalize() {
                    for (auto &hist: atmdhist) {
                        double sum = hist.second.sumy();
                        for (auto &i : hist.second.getMap())
                            i.second = i.second/sum ;
                    }
                    for (auto &hist: moldhist) {
                        double sum = hist.second.sumy();
                        for (auto &i : hist.second.getMap())
                            i.second = i.second/sum ;
                    }
                }
            };

        template<class Tspace>
            class Multipole : public Analysisbase {
                const Tspace& spc;
                struct data { Average<double> Z, Z2, mu, mu2; };
                std::map<int,data> _map; //!< Molecular moments and their fluctuations

                void _sample() override {
                    for (auto &g : spc.groups)
                        if (!g.atomic) {
                            auto &d = _map[g.id];
                            auto p = Geometry::toMultipole(g, spc.geo.getBoundaryFunc());
                            d.Z += p.charge;
                            d.mu += p.mulen;
                            d.Z2 += p.charge*p.charge;
                            d.mu2 += p.mulen*p.mulen;
                        }
                }

                void _to_json(json &j) const override {
                    json &k = j["molecules"];
                    for (auto &d : _map)
                        k[ molecules<typename Tspace::Tpvec>[d.first].name ] = {
                            {"Z", d.second.Z.avg()}, {"Z2", d.second.Z2.avg()},
                            {"C", std::pow(d.second.Z.avg(),2)-d.second.Z2.avg()},
                            {u8::mu, d.second.mu.avg()}, {u8::mu+u8::squared, d.second.mu2.avg()}
                        };
                }

                public:
                Multipole(const json &j, const Tspace &spc ) : spc(spc) {
                    from_json(j);
                    name = "multipole";
                }
            }; // Molecular multipoles and their fluctuations

        class SystemEnergy : public Analysisbase {
            std::string file, sep=" ";
            std::ofstream f;
            std::function<std::vector<double>()> energyFunc;
            Average<double> uavg, u2avg; //!< mean energy and mean squared energy
            std::vector<std::string> names;
            Table2D<double,double> ehist;  //Density histograms

            double uinit;

            void normalize();

            void _sample() override;

            void _to_json(json &j) const override;

            void _from_json(const json &j) override;

            public:
            template<class Tspace>
                SystemEnergy( const json &j, Energy::Hamiltonian<Tspace> &pot ) {
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
                            u.push_back( i->energy(change) );
                        return u;
                    };
                    ehist.setResolution(0.25);
                    auto u = energyFunc();
                    uinit = std::accumulate(u.begin(), u.end(), 0.0); // initial energy
                }
        }; //!< Save system energy to disk. Keywords: `nstep`, `file`.

        /**
         * @brief Checks if system is sane. If not, abort program.
         */
        template<class Tspace>
            class SanityCheck : public Analysisbase {
                private:
                    Tspace &spc;
                    void _sample() override {
                        // loop over all groups
                        for (auto &g : spc.groups) {
                            // check if particles are inside container
                            for (auto &i : g) // loop over active particles
                                if (spc.geo.collision(i.pos))
                                    throw std::runtime_error("particle outside container");

                            // check if molecular mass centers are correct
                            if (not g.atomic)
                                if (not g.empty()) {
                                    Point cm =  Geometry::massCenter( g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.cm);
                                    double sqd = spc.geo.sqdist(g.cm, cm);
                                    if (sqd>1e-6) {
                                        std::cerr
                                            << "dist:      " << sqrt(sqd) << endl
                                            << "g.cm:      " << g.cm.transpose() << endl
                                            << "actual cm: " << cm.transpose() << endl;
                                        throw std::runtime_error("mass center-out-of-sync");
                                    }
                                }
                        }
                    }
                public:
                    SanityCheck(const json &j, Tspace &spc) : spc(spc) {
                        from_json(j);
                        name = "sanity";
                        steps = j.value("nstep", -1);
                    }
            };

        class SaveState : public Analysisbase {
            private:
                std::function<void(std::string)> writeFunc = nullptr;
                std::string file;
                void _to_json(json &j) const override;
                void _sample() override;
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

                        else if ( suffix == "gro" )
                            writeFunc = std::bind(
                                    []( std::string file, Tspace &s ) { FormatGRO::save(file, s); },
                                    _1, std::ref(spc));

                        else if ( suffix == "pqr" )
                            writeFunc = std::bind(
                                    []( std::string file, Tspace &s ) { FormatPQR::save(file, s.p, s.geo.getLength()); },
                                    _1, std::ref(spc));

                        else if ( suffix == "xyz" )
                            writeFunc = std::bind(
                                    []( std::string file, Tspace &s ) { FormatXYZ::save(file, s.p, s.geo.getLength()); },
                                    _1, std::ref(spc));

                        else if ( suffix == "json" ) // JSON state file
                            writeFunc = [&spc](const std::string &file) {
                                std::ofstream f(file);
                                if (f) {
                                    json j;
                                    Faunus::to_json(j, spc);
                                    j["random-move"] = Move::Movebase::slump;
                                    j["random-global"] = Faunus::random;
                                    f << std::setw(2) << j;
                                }
                            };

                        else if ( suffix == "ubj" ) // Universal Binary JSON state file
                            writeFunc = [&spc](const std::string &file) {
                                std::ofstream f(file, std::ios::binary);
                                if (f) {
                                    json j;
                                    Faunus::to_json(j, spc);
                                    j["random-move"] = Move::Movebase::slump;
                                    j["random-global"] = Faunus::random;
                                    auto v = json::to_ubjson(j); // json --> binary
                                    f.write( (const char*)v.data(), v.size()*sizeof(decltype(v)::value_type));
                                }
                            };

                        if (writeFunc==nullptr)
                            throw std::runtime_error("unknown file extension for '"+file+"'");
                    }

                ~SaveState();
        };

        /**
         * @brief Base class for distribution functions etc.
         */
        class PairFunctionBase : public Analysisbase {
            protected:
                int dim=3;
                int id1=-1, id2=-1; // particle id (mol or atom)
                double dr=0;
                Equidistant2DTable<double,double> hist;
                std::string name1, name2, file;
                double Rhypersphere=-1; // Radius of 2D hypersphere
                Average<double> V;   // average volume (angstrom^3)

            private:
                void _from_json(const json &j) override;
                void _to_json(json &j) const override;

            public:
                PairFunctionBase(const json &j);
                virtual ~PairFunctionBase();
        };

        /** @brief Atomic radial distribution function, g(r) */
        template<class Tspace>
            class AtomRDF : public PairFunctionBase {
                Tspace &spc;

                void _sample() override {
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

                    auto it = findName( atoms, name1 );
                    if ( it == atoms.end() )
                        throw std::runtime_error("unknown atom '" + name1 + "'");
                    id1 = it->id();

                    it = findName( atoms, name2 );
                    if ( it == atoms.end() )
                        throw std::runtime_error("unknown atom '" + name2 + "'");
                    id2 = it->id();
                }
            };

        /** @brief Same as `AtomRDF` but for molecules. Identical input. */
        template<class Tspace>
            class MoleculeRDF : public PairFunctionBase {
                typedef typename Tspace::Tpvec Tpvec;
                Tspace &spc;

                void _sample() override {
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

        /** @brief Write XTC trajectory file */
        template<class Tspace>
            class XTCtraj : public Analysisbase {
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

                void _sample() override {
                    xtc.setbox( spc.geo.getLength() );
                    bool rc = xtc.save(file, spc.p, Group<Tparticle>(spc.p.begin(), spc.p.end()));
                    if (rc==false)
                        std::cerr << "error saving xtc\n";
                }

                public:

                XTCtraj( const json &j, Tspace &s ) : xtc(1e6), spc(s) {
                    from_json(j);
                    name = "xtcfile";
                    cite = "http://bit.ly/2A8lzpa";
                }
            };

        class VirtualVolume : public Analysisbase {
            double dV;
            Change c;
            Energy::Energybase& pot;
            std::function<double()> getVolume;
            std::function<void(double)> scaleVolume;
            Average<double> duexp; // < exp(-du/kT) >

            void _sample() override;
            void _from_json(const json &j) override;
            void _to_json(json &j) const override;

            public:
            template<class Tspace>
                VirtualVolume(const json &j, Tspace &spc, Energy::Energybase &pot) : pot(pot) {
                    from_json(j);
                    c.dV=true;
                    c.all=true;
                    name = "virtualvolume";
                    cite = "doi:10.1063/1.472721";
                    getVolume = [&spc](){ return spc.geo.getVolume(); };
                    scaleVolume = [&spc](double Vnew) { spc.scaleVolume(Vnew); };
                }
        }; //!< Excess pressure using virtual volume move

        /**
         * @brief Multipolar decomposition between groups as a function of separation
         * @date Malmo 2014
         * @todo Add option to use charge center instead of mass center
         */
        template<class Tspace>
            class MultipoleDistribution : public Analysisbase {
                typedef typename Tspace::Tgroup Tgroup;
                typedef typename Tspace::Tparticle Tparticle;

                struct data {
                    Average<double> tot, ii, id, iq, dd, mucorr;
                };

                std::vector<std::string> names; //!< Molecule names (len=2)
                std::vector<int> ids;           //!< Molecule ids (len=2)
                std::string filename;           //!< output file name
                int id1, id2;                   //!< pair of molecular id's to analyse
                double dr;                      //!< distance resolution
                std::map<int, data> m;          //!< Energy distributions
                Tspace &spc;

                double g2g(const Tgroup &g1, const Tgroup &g2) {
                    double u = 0;
                    for ( auto &i : g1 )
                        for ( auto &j : g2 )
                            u += i.charge * j.charge / spc.geo.vdist(i.pos, j.pos).norm();
                    return u;
                } //<! exact ion-ion energy between particles

                void save() const {
                    using std::setw;
                    if (cnt>0) {
                        std::ofstream f(MPI::prefix + filename.c_str());
                        if (f) {
                            char w = 12;
                            f.precision(4);
                            f << "# Multipolar energies (kT/lB)\n"
                                << std::left << setw(w) << "# R/AA" << std::right << setw(w) << "exact"
                                << setw(w) << "total" << setw(w) << "ionion" << setw(w) << "iondip"
                                << setw(w) << "dipdip" << setw(w) << "ionquad" << setw(w) << "mucorr\n";
                            for (auto &i : m)
                                f << std::left << setw(w) << i.first * dr << std::right
                                    << setw(w) << i.second.tot << setw(w)
                                    << i.second.ii.avg()+i.second.id.avg()+i.second.dd.avg()+i.second.iq.avg()
                                    << setw(w) << i.second.ii << setw(w) << i.second.id
                                    << setw(w) << i.second.dd << setw(w) << i.second.iq
                                    << setw(w) << i.second.mucorr << "\n";
                        }
                    }
                } //!< save to disk

                void _sample() override {
                    for (auto &gi : spc.findMolecules(ids[0]))
                        for (auto &gj : spc.findMolecules(ids[1]))
                            if (gi!=gj) {
                                auto a = Geometry::toMultipole(gi, spc.geo.getBoundaryFunc());
                                auto b = Geometry::toMultipole(gj, spc.geo.getBoundaryFunc());
                                Point R = spc.geo.vdist(gi.cm, gj.cm);
                                auto &d = m[ to_bin(R.norm(), dr) ];
                                d.tot += g2g(gi, gj);
                                d.ii  += a.charge * b.charge / R.norm();
                                d.id  += q2mu( a.charge * b.mulen, b.mu, b.charge * a.mulen, a.mu, R);
                                d.dd  += mu2mu( a.mu, b.mu, a.mulen * b.mulen, R );
                                d.iq  += q2quad( a.charge, b.Q, b.charge, a.Q, R );
                                d.mucorr += a.mu.dot(b.mu);
                            }
                }

                void _to_json(json &j) const override {
                    j = {{"molecules", names}, {"file", filename}, {"dr", dr}};
                }

                public:
                MultipoleDistribution(const json &j, Tspace &spc) : spc(spc) {
                    from_json(j);
                    name = "Multipole Distribution";
                    dr = j.value("dr", 0.2);
                    filename = j.at("file").get<std::string>();
                    names = j.at("molecules").get<decltype(names)>(); // molecule names
                    ids = names2ids(molecules<typename Tspace::Tpvec>, names);// names --> molids
                    if (ids.size()!=2)
                        throw std::runtime_error("specify exactly two molecules");
                }

                ~MultipoleDistribution() { save(); }

            }; // end of multipole distribution

        /** @brief Sample scattering intensity */
        template<class Tspace, class Tformfactor=Scatter::FormFactorUnity<double>>
            class ScatteringFunction : public Analysisbase {
                Tspace &spc;
                bool usecom;                    // scatter from mass center, only?
                std::string filename;           // output file name
                std::vector<Point> p;           // vector of scattering points
                std::vector<int> ids;           // Molecule ids
                std::vector<std::string> names; // Molecule names
                Scatter::DebyeFormula<Tformfactor> debye;

                void _sample() override {
                    p.clear();
                    for (int id : ids) { // loop over molecule names
                        auto groups = spc.findMolecules(id);
                        for (auto &g : groups) // loop over groups
                            if (usecom && !g.atomic)
                                p.push_back( g.cm );
                            else
                                for (auto &i : g) // loop over particle index in group
                                    p.push_back( i.pos );
                    }
                    debye.sample( p, spc.geo.getVolume() );
                }

                void _to_json(json &j) const override {
                    j = { { "molecules", names }, { "com", usecom } };
                }

                public:
                ScatteringFunction(const json &j, Tspace &spc) try : spc(spc), debye(j) {
                    from_json(j);
                    name = "scatter";
                    usecom = j.value("com", true);
                    filename = j.at("file").get<std::string>();
                    names = j.at("molecules").get<decltype(names)>(); // molecule names
                    ids = names2ids(molecules<typename Tspace::Tpvec>, names);// names --> molids
                }
                catch( std::exception &e ) {
                    std::cerr << "Debye Formula Scattering: ";
                    throw;
                }

                ~ScatteringFunction() {
                    debye.save( filename );
                }
            };

        /**
         * @brief Analysis of polymer shape - radius of gyration, shape factor etc.
         * @date November, 2011
         *
         * This will analyse polymer groups and calculate Rg, Re and the shape factor. If
         * sample() is called with different groups these will be distinguished by their
         * *name* and sampled individually.
         */
        template<class Tspace>
            class PolymerShape : public Analysisbase {
                typedef typename Tspace::Tparticle Tparticle;
                Tspace &spc;
                std::map<int, Average<double>> Rg2, Rg, Re2, Re, Rs, Rs2, Rg2x, Rg2y, Rg2z;
                std::vector<int> ids; // molecule id's to analyse

                void _to_json(json &j) const override {
                    using namespace u8;
                    json &k = j["molecules"];
                    for (int i : ids)
                        k[ molecules<typename Tspace::Tpvec>[i].name ] = {
                            { bracket("Rg"), Rg.at(i).avg() },
                            { bracket("Re"), Re.at(i).avg() },
                            { bracket("Rg" + squared), Rg2.at(i).avg() },
                            { bracket("Rg" + squared) + "-" + bracket("Rg") + squared, Rg2.at(i).avg() - std::pow(Rg.at(i).avg(), 2.0) },
                            { bracket("Re" + squared) + "/" + bracket("Rg" + squared), Re2.at(i).avg() / Rg2.at(i).avg() },
                            { rootof + bracket("Rg"  + squared), sqrt( Rg2.at(i).avg() ) },
                            { rootof + bracket("Re"  + squared), sqrt( Re2.at(i).avg() )  },
                            { rootof + bracket("Rgxyz" + squared),
                                { 
                                    sqrt( Rg2x.at(i).avg() ),
                                    sqrt( Rg2y.at(i).avg() ),
                                    sqrt( Rg2z.at(i).avg() )
                                }
                            }
                        };
                }

                Point vectorgyrationRadiusSquared(typename Tspace::Tgroup &g) const {
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

                void _sample() override {
                    for (int i : ids)
                        for (auto &g : spc.findMolecules(i))
                            if (g.size()>1) {
                                Point r2 = vectorgyrationRadiusSquared(g);
                                double rg2 = r2.sum();
                                double re2 = spc.geo.sqdist( g.begin()->pos, (g.end()-1)->pos );
                                Rg[i] += sqrt(rg2);
                                Re[i] += sqrt(re2);
                                Rg2[i] += rg2;
                                Re2[i] += re2; //end-2-end squared
                                Rg2x[i] += r2.x();
                                Rg2y[i] += r2.y();
                                Rg2z[i] += r2.z();
                                double rs = Re2[i].avg() / Rg2[i].avg(); // fluctuations in shape factor
                                Rs[i] += rs;
                                Rs2[i] += rs * rs;
                            }
                }

                public:

                PolymerShape(const json &j, Tspace &spc) : spc(spc) {
                    from_json(j);
                    name = "Polymer Shape";
                    auto names = j.at("molecules").get<std::vector<std::string>>(); // molecule names
                    ids = names2ids(molecules<typename Tspace::Tpvec>, names);// names --> molids
                }

            };

        /** 
         * @brief "Trajectory" with charge and radius, only, for all (active, inactive) particles
         *
         * For use w. VMD to visualize charge fluctuations and grand canonical ensembles
         */
        class QRtraj : public Analysisbase {
            private:
                std::string file;
                std::ofstream f;
                std::function<void()> write_to_file;
                void _sample() override;
                void _to_json(json &j) const override;

            public:
                template<class Tspace>
                    QRtraj(const json &j, Tspace &spc) {
                        from_json(j);
                        name = "qrfile";
                        file = j.value("file", "qrtraj.dat"s);
                        f.open(file);
                        if (not f)
                            throw std::runtime_error("error opening "s + file);
                        f.precision(6);
                        write_to_file = [&groups=spc.groups, &f=f]() {
                            for (auto &g : groups)
                                for (auto it=g.begin(); it!=g.trueend(); ++it) // loop over *all* particles
                                    if (it<g.end())
                                        f << it->charge << " " << atoms[it->id].sigma*0.5 << " "; 
                                    else
                                        f << "0 0 "; // zero charge and radii for inactive particles
                            f << "\n"; // newline for every frame
                        };
                    }
        };

        struct CombinedAnalysis : public BasePointerVector<Analysisbase> {
            template<class Tspace, class Tenergy>
                CombinedAnalysis(const json &j, Tspace &spc, Tenergy &pot) {
                    if (j.is_array())
                        for (auto &m : j)
                            for (auto it=m.begin(); it!=m.end(); ++it)
                                if (it->is_object())
                                    try {
                                        size_t oldsize = this->vec.size();
                                        if      (it.key()=="atomprofile") push_back<AtomProfile<Tspace>>(it.value(), spc);
                                        else if (it.key()=="atomrdf") push_back<AtomRDF<Tspace>>(it.value(), spc);
                                        else if (it.key()=="density") push_back<Density<Tspace>>(it.value(), spc);
                                        else if (it.key()=="molrdf") push_back<MoleculeRDF<Tspace>>(it.value(), spc);
                                        else if (it.key()=="multipole") push_back<Multipole<Tspace>>(it.value(), spc);
                                        else if (it.key()=="multipoledist") push_back<MultipoleDistribution<Tspace>>(it.value(), spc);
                                        else if (it.key()=="polymershape") push_back<PolymerShape<Tspace>>(it.value(), spc);
                                        else if (it.key()=="qrfile") push_back<QRtraj>(it.value(), spc);
                                        else if (it.key()=="reactioncoordinate") push_back<FileReactionCoordinate>(it.value(), spc);
                                        else if (it.key()=="sanity") push_back<SanityCheck<Tspace>>(it.value(), spc);
                                        else if (it.key()=="savestate") push_back<SaveState>(it.value(), spc);
                                        else if (it.key()=="scatter") push_back<ScatteringFunction<Tspace>>(it.value(), spc);
                                        else if (it.key()=="sliceddensity") push_back<SlicedDensity<Tspace>>(it.value(), spc);
                                        else if (it.key()=="systemenergy") push_back<SystemEnergy>(it.value(), pot);
                                        else if (it.key()=="virtualvolume") push_back<VirtualVolume>(it.value(), spc, pot);
                                        else if (it.key()=="widom") push_back<WidomInsertion<Tspace>>(it.value(), spc, pot);
                                        else if (it.key()=="xtcfile") push_back<XTCtraj<Tspace>>(it.value(), spc);
                                        // additional analysis go here...
 
                                        if (this->vec.size()==oldsize)
                                            throw std::runtime_error("unknown analysis");

                                    } catch (std::exception &e) {
                                        throw std::runtime_error("Error adding analysis,\n\n\"" + it.key() + "\": "
                                                + it->dump() + "\n\n: " + e.what() + usageTip[it.key()]);
                                    }
                }

            void sample();

        }; //!< Aggregates analysis

        /** @brief Example analysis */
        template<class T, class Enable = void>
            struct _analyse {
                void sample(T&) {
                    std::cout << "not a dipole!" << std::endl;
                } //!< Sample
            }; // primary template

        /** @brief Example analysis */
        template<class T>
            struct _analyse<T, typename std::enable_if<std::is_base_of<Dipole, T>::value>::type> {
                void sample(T&) {
                    std::cout << "dipole!" << std::endl;
                } //!< Sample
            }; // specialized template

    }//namespace

}//namespace
