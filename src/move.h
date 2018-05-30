#pragma once

#include "core.h"
#include "energy.h"
#include "average.h"
#include "analysis.h"
#include "potentials.h"
#include "mpi.h"

namespace Faunus {
    namespace Move {

        class Movebase {
            private:
                virtual void _move(Change&)=0; //!< Perform move and modify change object
                virtual void _accept(Change&) {}; //!< Call after move is accepted
                virtual void _reject(Change&) {}; //!< Call after move is rejected
                TimeRelativeOfTotal<std::chrono::microseconds> timer;
            protected:
                virtual void _to_json(json &j) const=0; //!< Extra info for report if needed
                virtual void _from_json(const json &j)=0; //!< Extra info for report if needed
                unsigned long cnt=0;
                unsigned long accepted=0;
                unsigned long rejected=0;
            public:
                static Random slump;   //!< Shared for all moves
                std::string name;      //!< Name of move
                std::string cite;      //!< Reference
                int repeat=1;          //!< How many times the move should be repeated per sweep

                inline void from_json(const json &j) {
                    auto it = j.find("repeat");
                    if (it!=j.end()) {
                        if (it->is_number())
                            repeat = it->get<double>();
                        else
                            if (it->is_string())
                                if (it->get<std::string>()=="N")
                                    repeat = -1;
                    }
                    _from_json(j);
                    if (repeat<0)
                        repeat=0;
                }

                inline void to_json(json &j) const {
                    _to_json(j);
                    j["relative time"] = timer.result();
                    j["acceptance"] = double(accepted)/cnt;
                    j["repeat"] = repeat;
                    j["moves"] = cnt;
                    if (!cite.empty())
                        j["cite"] = cite;
                    _roundjson(j, 3);
                } //!< JSON report w. statistics, output etc.

                inline void move(Change &change) {
                    timer.start();
                    cnt++;
                    change.clear();
                    _move(change);
                    if (change.empty())
                        timer.stop();
                } //!< Perform move and modify given change object

                inline void accept(Change &c) {
                    accepted++;
                    _accept(c);
                    timer.stop();
                }

                inline void reject(Change &c) {
                    rejected++;
                    _reject(c);
                    timer.stop();
                }

                inline virtual double bias(Change &c, double uold, double unew) {
                    return 0; // du
                } //!< adds extra energy change not captured by the Hamiltonian
        };

        Random Movebase::slump; // static instance of Random (shared for all moves)

        inline void from_json(const json &j, Movebase &m) {
            m.from_json( j );
        } //!< Configure any move via json

        inline void to_json(json &j, const Movebase &m) {
            assert( !m.name.empty() );
            m.to_json(j[m.name]);
        }

        /**
         * @brief Translate and rotate a molecular group
         */
        template<typename Tspace>
            class AtomicTranslateRotate : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tparticle Tparticle;
                    Tspace& spc; // Space to operate on
                    int molid=-1;
                    Point dir={1,1,1};
                    Average<double> msqd; // mean squared displacement
                    double _sqd; // squared displament
                    std::string molname; // name of molecule to operate on
                    Change::data cdata;

                    void _to_json(json &j) const override {
                        j = {
                            {"dir", dir},
                            {"molid", molid},
                            {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
                            {"molecule", molname}
                        };
                        _roundjson(j,3);
                    }

                    void _from_json(const json &j) override {
                        assert(!molecules<Tpvec>.empty());
                        try {
                            molname = j.at("molecule");
                            auto it = findName(molecules<Tpvec>, molname);
                            if (it == molecules<Tpvec>.end())
                                throw std::runtime_error("unknown molecule '" + molname + "'");
                            molid = it->id();
                            dir = j.value("dir", Point(1,1,1));
                            if (repeat<0) {
                                auto v = spc.findMolecules(molid);
                                repeat = std::distance(v.begin(), v.end()); // repeat for each molecule...
                                if (repeat>0)
                                    repeat = repeat * v.front().size();     // ...and for each atom
                            }
                        }
                        catch (std::exception &e) {
                            std::cerr << name << ": " << e.what();
                            throw;
                        }
                    } //!< Configure via json object

                    typename Tpvec::iterator randomAtom() {
                        assert(molid>=0);
                        auto mollist = spc.findMolecules( molid ); // all `molid` groups
                        if (size(mollist)>0) {
                            auto git = slump.sample( mollist.begin(), mollist.end() ); // random molecule iterator
                            if (!git->empty()) {
                                auto p = slump.sample( git->begin(), git->end() ); // random particle iterator  
                                cdata.index = Faunus::distance( spc.groups.begin(), git ); // integer *index* of moved group
                                cdata.atoms[0] = std::distance(git->begin(), p);  // index of particle rel. to group
                                return p; 
                            }
                        }
                        return spc.p.end();
                    }

                    void _move(Change &change) override {
                        auto p = randomAtom();
                        if (p!=spc.p.end()) {
                            double dp = atoms<Tparticle>.at(p->id).dp;
                            double dprot = atoms<Tparticle>.at(p->id).dprot;
                            auto& g = spc.groups[cdata.index];

                            if (dp>0) { // translate
                                Point oldpos = p->pos;
                                p->pos +=  0.5 * dp * ranunit(slump).cwiseProduct(dir);
                                spc.geo.boundaryFunc(p->pos);
                                _sqd = spc.geo.sqdist(oldpos, p->pos); // squared displacement
                                if (!g.atomic)
                                    g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.boundaryFunc, -g.cm);
                            }

                            if (dprot>0) { // rotate
                                Point u = ranunit(slump);
                                double angle = dprot * (slump()-0.5);
                                Eigen::Quaterniond Q( Eigen::AngleAxisd(angle, u) );
                                p->rotate(Q, Q.toRotationMatrix());
                            }

                            if (dp>0 || dprot>0)
                                change.groups.push_back( cdata ); // add to list of moved groups
                         }
                        else
                            std::cerr << name << ": no atoms found" << std::endl;
                    }

                    void _accept(Change &change) override { msqd += _sqd; }
                    void _reject(Change &change) override { msqd += 0; }

                public:
                    AtomicTranslateRotate(Tspace &spc) : spc(spc) {
                        name = "transrot";
                        repeat = -1; // meaning repeat N times
                        cdata.atoms.resize(1);
                        cdata.internal=true;
                    }
            };

        /**
         * @brief Translate and rotate a molecular group
         */
        template<typename Tspace>
            class TranslateRotate : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace& spc; // Space to operate on
                    int molid=-1;
                    double dptrans=0;
                    double dprot=0;
                    Point dir={1,1,1};
                    Average<double> msqd; // mean squared displacement
                    double _sqd; // squared displament

                    void _to_json(json &j) const override {
                        j = {
                            {"dir", dir}, {"dp", dptrans}, {"dprot", dprot},
                            {"molid", molid},
                            {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
                            {"molecule", molecules<Tpvec>[molid].name}
                        };
                        _roundjson(j,3);
                    }

                    void _from_json(const json &j) override {
                        assert(!molecules<Tpvec>.empty());
                        try {
                            std::string molname = j.at("molecule");
                            auto it = findName(molecules<Tpvec>, molname);
                            if (it == molecules<Tpvec>.end())
                                throw std::runtime_error("unknown molecule '" + molname + "'");
                            molid = it->id();
                            dir = j.value("dir", Point(1,1,1));
                            dprot = j.at("dprot");
                            dptrans = j.at("dp");
                            if (repeat<0) {
                                auto v = spc.findMolecules(molid);
                                repeat = std::distance(v.begin(), v.end());
                            }
                        }
                        catch (std::exception &e) {
                            std::cerr << name << ": " << e.what();
                            throw;
                        }
                    } //!< Configure via json object

                    void _move(Change &change) override {
                        assert(molid>=0);
                        assert(!spc.groups.empty());
                        assert(spc.geo.getVolume()>0);

                        // pick random group from the system matching molecule type
                        // TODO: This can be slow -- implement look-up-table in Space
                        auto mollist = spc.findMolecules( molid ); // list of molecules w. 'molid'
                        if (size(mollist)>0) {
                            auto it = slump.sample( mollist.begin(), mollist.end() );
                            if (!it->empty()) {
                                assert(it->id==molid);

                                if (dptrans>0) { // translate
                                    Point oldcm = it->cm;
                                    Point dp = 0.5*ranunit(slump).cwiseProduct(dir) * dptrans;
                                    it->translate( dp, spc.geo.boundaryFunc );
                                    _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
                                }

                                if (dprot>0) { // rotate
                                    Point u = ranunit(slump);
                                    double angle = dprot * (slump()-0.5);
                                    Eigen::Quaterniond Q( Eigen::AngleAxisd(angle, u) );
                                    it->rotate(Q, spc.geo.boundaryFunc);
                                }

                                if (dptrans>0|| dprot>0) { // define changes
                                    Change::data d;
                                    d.index = Faunus::distance( spc.groups.begin(), it ); // integer *index* of moved group
                                    d.all = true; // *all* atoms in group were moved
                                    change.groups.push_back( d ); // add to list of moved groups
                                }
                                assert( spc.geo.sqdist( it->cm,
                                            Geometry::massCenter(it->begin(),it->end(),spc.geo.boundaryFunc,-it->cm) ) < 1e-9 );
                            }
                        }
                        else std::cerr << name << ": no molecules found" << std::endl;
                    }

                    void _accept(Change &change) override { msqd += _sqd; }
                    void _reject(Change &change) override { msqd += 0; }

                public:
                    TranslateRotate(Tspace &spc) : spc(spc) {
                        name = "moltransrot";
                        repeat = -1; // meaning repeat N times
                    }
            };

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] TranslateRotate")
        {
            typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
            typedef Space<Geometry::Cuboid, Tparticle> Tspace;
            typedef typename Tspace::Tpvec Tpvec;

            CHECK( !atoms<Tparticle>.empty() ); // set in a previous test
            CHECK( !molecules<Tpvec>.empty() ); // set in a previous test

            Tspace spc;
            TranslateRotate<Tspace> mv(spc);
            json j = R"( {"molecule":"B", "dp":1.0, "dprot":0.5, "dir":[0,1,0], "repeat":2 })"_json;
            mv.from_json(j);

            j = json(mv).at(mv.name);
            CHECK( j.at("molecule")   == "B");
            CHECK( j.at("dir")   == Point(0,1,0) );
            CHECK( j.at("dp")    == 1.0 );
            CHECK( j.at("repeat")  == 2 );
            CHECK( j.at("dprot") == 0.5 );
        }
#endif

        template<typename Tspace>
            class VolumeMove : public Movebase {
                private:
                    const std::map<std::string, Geometry::VolumeMethod> methods = {
                        {"xy", Geometry::XY},
                        {"isotropic", Geometry::ISOTROPIC},
                        {"isochoric", Geometry::ISOCHORIC}
                    };
                    typename decltype(methods)::const_iterator method;
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace& spc;
                    Average<double> msqd; // mean squared displacement
                    double dV=0, deltaV=0, Vnew=0, Vold=0;

                    void _to_json(json &j) const override {
                        using namespace u8;
                        j = {
                            {"dV", dV}, {"method", method->first},
                            {rootof + bracket(Delta + "V" + squared), std::sqrt(msqd.avg())},
                            {cuberoot + rootof + bracket(Delta + "V" + squared),
                                std::cbrt(std::sqrt(msqd.avg()))}
                        };
                        _roundjson(j,3);
                    }

                    void _from_json(const json &j) override {
                        method = methods.find( j.value("method", "isotropic") );
                        if (method==methods.end())
                            std::runtime_error("unknown volume change method");
                        dV = j.at("dV");
                    }

                    void _move(Change &change) override {
                        if (dV>0) {
                            change.dV=true;
                            change.all=true;
                            Vold = spc.geo.getVolume();
                            if (method->second == Geometry::ISOCHORIC)
                                Vold = std::pow(Vold,1.0/3.0); // volume is constant
                            Vnew = std::exp(std::log(Vold) + (slump()-0.5) * dV);
                            deltaV = Vnew-Vold;
                            spc.scaleVolume(Vnew, method->second);
                        } else deltaV=0;
                    }

                    void _accept(Change &change) override { msqd += deltaV*deltaV; }
                    void _reject(Change &change) override { msqd += 0; }

                public:
                    VolumeMove(Tspace &spc) : spc(spc) {
                        name = "volume";
                        repeat = 1;
                    }
            };

        template<typename Tspace>
            class Cluster : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace& spc;
                    Average<double> msqd, msqd_angle, N;
                    double thresholdsq=0, dptrans=0, dprot=0, angle=0, _bias=0;
                    Point dir={1,1,1}, dp;
                    std::vector<std::string> names;
                    std::vector<int> ids;
                    std::vector<size_t> index; // all possible molecules to move

                    void _to_json(json &j) const override {
                        using namespace u8;
                        j = {
                            {"threshold", std::sqrt(thresholdsq)}, {"dir", dir}, {"dp", dptrans}, {"dprot", dprot},
                            {rootof + bracket("r" + squared), std::sqrt(msqd.avg())},
                            {rootof + bracket(theta + squared) + "/" + degrees, std::sqrt(msqd_angle.avg()) / 1.0_deg},
                            {bracket("N"), N.avg()}
                        };
                        _roundjson(j,3);
                    }

                    void _from_json(const json &j) override {
                        dptrans = j.at("dp");
                        dir = j.value("dir", Point(1,1,1));
                        dprot = j.at("dprot");
                        thresholdsq = std::pow(j.at("threshold").get<double>(), 2);
                        names = j.at("molecules").get<decltype(names)>(); // molecule names
                        ids = names2ids(molecules<Tpvec>, names);     // names --> molids
                        index.clear();
                        for (auto &g : spc.groups)
                            if (!g.atomic)
                                if (std::find(ids.begin(), ids.end(), g.id)!=ids.end() )
                                    index.push_back( &g-&spc.groups.front() );
                        if (repeat<0)
                            repeat = index.size();
                    }

                    void findCluster(Tspace &spc, size_t first, std::set<size_t>& cluster) {
                        std::set<size_t> pool(index.begin(), index.end());
                        cluster.clear();
                        cluster.insert(first);
                        pool.erase(first);
                        size_t n;
                        do { // find cluster (not very clever...)
                            n = cluster.size();
                            for (size_t i : cluster)
                                if (!spc.groups[i].empty()) // check if group is inactive
                                    for (size_t j : pool)
                                        if (!spc.groups[j].empty()) // check if group is inactive
                                            if (i!=j)
                                                if (spc.geo.sqdist(spc.groups[i].cm, spc.groups[j].cm)<=thresholdsq) {
                                                    cluster.insert(j);
                                                    pool.erase(j);
                                                }
                        } while (cluster.size()!=n);

                        // check if cluster is too large
                        double max = spc.geo.getLength().minCoeff()/2;
                        for (auto i : cluster)
                            for (auto j : cluster)
                                if (j>i)
                                    if (spc.geo.sqdist(spc.groups[i].cm, spc.groups[j].cm)>=max*max)
                                        throw std::runtime_error(name+": cluster larger than half box length");
                    }

                    void _move(Change &change) override {
                        if (thresholdsq>0 && !index.empty()) {
                            std::set<size_t> cluster; // all group index in cluster
                            size_t first = *slump.sample(index.begin(), index.end()); // random molecule (nuclei)
                            findCluster(spc, first, cluster); // find cluster around first

                            N += cluster.size(); // average cluster size
                            Change::data d;
                            d.all=true;
                            dp = 0.5*ranunit(slump).cwiseProduct(dir) * dptrans;
                            angle = dprot * (slump()-0.5);

                            Point COM = Geometry::trigoCom(spc, cluster); // cluster center
                            Eigen::Quaterniond Q;
                            Q = Eigen::AngleAxisd(angle, ranunit(slump)); // quaternion

                            for (auto i : cluster) { // loop over molecules in cluster
                                auto &g = spc.groups[i];

                                Geometry::rotate(g.begin(), g.end(), Q, spc.geo.boundaryFunc, -COM);
                                g.cm = g.cm-COM;
                                spc.geo.boundary(g.cm);
                                g.cm = Q*g.cm+COM;
                                spc.geo.boundary(g.cm);

                                g.translate( dp, spc.geo.boundaryFunc );
                                d.index=i;
                                change.groups.push_back(d);
                            }
                            _bias += 0; // one may add bias here...
#ifndef NDEBUG
                            Point newCOM = Geometry::trigoCom(spc, cluster);
                            double _zero = std::sqrt( spc.geo.sqdist(COM,newCOM) ) - dp.norm();
                            if (fabs(_zero)>1)
                                std::cerr << _zero << " ";
#endif
                        }
                    }

                    double bias(Change &change, double uold, double unew) override {
                        return _bias;
                    } //!< adds extra energy change not captured by the Hamiltonian

                    void _reject(Change &change) override { msqd += 0; msqd_angle += 0; }

                    void _accept(Change &change) override {
                        msqd += dp.squaredNorm();
                        msqd_angle += angle*angle;
                    }

                public:
                    Cluster(Tspace &spc) : spc(spc) {
                        cite = "doi:10/cj9gnn";
                        name = "cluster";
                        repeat = -1; // meaning repeat N times
                    }
            };

        template<typename Tspace>
            class Pivot : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    std::vector<std::reference_wrapper<const Potential::BondData>> bonds;
                    std::vector<int> index; // atom index to rotate
                    Tspace& spc;
                    std::string molname;
                    int molid;
                    double dprot;
                    double d2; // cm movement, squared
                    Average<double> msqd; // cm mean squared displacement

                    void _to_json(json &j) const override {
                        using namespace u8;
                        j = {
                            {"molecule", molname}, {"dprot", dprot},
                            {u8::rootof + u8::bracket("r_cm" + u8::squared), std::sqrt(msqd.avg())}
                        };
                        _roundjson(j,3);
                    }

                    void _from_json(const json &j) override {
                        dprot = j.at("dprot");
                        molname = j.at("molecule");
                        auto it = findName(molecules<Tpvec>, molname);
                        if (it == molecules<Tpvec>.end())
                            throw std::runtime_error("unknown molecule '" + molname + "'");
                        molid = it->id();
                        bonds = Potential::filterBonds(
                                molecules<Tpvec>[molid].bonds, Potential::BondData::harmonic);
                        if (repeat<0) {
                            auto v = spc.findMolecules(molid);
                            repeat = std::distance(v.begin(), v.end()); // repeat for each molecule...
                            if (repeat>0)
                                repeat *= bonds.size();
                        }
                    }

                    void _move(Change &change) override {
                        d2=0;
                        if (std::fabs(dprot)>1e-9) {
                            auto it = spc.randomMolecule(molid, slump);
                            if (it!=spc.groups.end())
                                if (it->size()>2) {
                                    auto b = slump.sample(bonds.begin(), bonds.end()); // random harmonic bond
                                    if (b != bonds.end()) {
                                        int i1 = b->get().index.at(0);
                                        int i2 = b->get().index.at(1);
                                        int offset = std::distance( spc.p.begin(), it->begin() );

                                        index.clear();
                                        if (slump()>0.0)
                                            for (size_t i=i2+1; i<it->size(); i++)
                                                index.push_back(i+offset);
                                        else
                                            for (int i=0; i<i1; i++)
                                                index.push_back(i+offset);
                                        i1+=offset;
                                        i2+=offset;

                                        if (!index.empty()) {
                                            Point oldcm = it->cm;
                                            it->unwrap(spc.geo.distanceFunc); // remove pbc
                                            Point u = (spc.p[i1].pos - spc.p[i2].pos).normalized();
                                            double angle = dprot * (slump()-0.5);
                                            Eigen::Quaterniond Q( Eigen::AngleAxisd(angle, u) );
                                            auto M = Q.toRotationMatrix();
                                            for (auto i : index) {
                                                spc.p[i].rotate(Q, M); // internal rot.
                                                spc.p[i].pos = Q * ( spc.p[i].pos - spc.p[i1].pos)
                                                    + spc.p[i1].pos; // positional rot.
                                            }
                                            it->cm = Geometry::massCenter(it->begin(), it->end());
                                            it->wrap(spc.geo.boundaryFunc); // re-apply pbc

                                            d2 = spc.geo.sqdist(it->cm, oldcm); // CM movement

                                            Change::data d;
                                            d.index = Faunus::distance( spc.groups.begin(), it ); // integer *index* of moved group
                                            d.all = true; // *all* atoms in group were moved
                                            change.groups.push_back( d ); // add to list of moved groups
                                        }
                                    }
                                }
                        }
                    }

                    void _accept(Change &change) override { msqd += d2; }
                    void _reject(Change &change) override { msqd += 0; }

                public:
                    Pivot(Tspace &spc) : spc(spc) {
                        name = "pivot";
                        repeat = -1; // --> repeat=N
                    }
            }; //!< Pivot move around random harmonic bond axis

#ifdef ENABLE_MPI
        /**
         * @brief Class for parallel tempering (aka replica exchange) using MPI
         *
         * Although not completely correct, the recommended way of performing a temper move
         * is to do `N` Monte Carlo passes with regular moves and then do a tempering move.
         * This is because the MPI nodes must be in sync and if you have a system where
         * the random number generator calls are influenced by the Hamiltonian we could
         * end up in a deadlock.
         *
         * @date Lund 2012, 2018
         */
        template<class Tspace>
            class ParallelTempering : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tparticle Tparticle;

                    Tspace& spc; // Space to operate on
                    MPI::MPIController& mpi;

                    int partner;                  //!< Exchange replica (partner)
                    enum extradata {VOLUME=0};    //!< Structure of extra data to send
                    std::map<std::string, Average<double>> accmap;

                    MPI::FloatTransmitter ft;   //!< Class for transmitting floats over MPI
                    MPI::ParticleTransmitter<Tpvec> pt;//!< Class for transmitting particles over MPI

                    void findPartner() {
                        int dr=0;
                        partner = mpi.rank();
                        (mpi.random()>0.5) ? dr++ : dr--;
                        (mpi.rank() % 2 == 0) ? partner+=dr : partner-=dr;
                    } //!< Find replica to exchange with

                    bool goodPartner() {
                        assert(partner!=mpi.rank() && "Selfpartner! This is not supposed to happen.");
                        if (partner>=0)
                            if ( partner<mpi.nproc() )
                                if ( partner!=mpi.rank() )
                                    return true;
                        return false;
                    } //!< Is partner valid?

                    void _to_json(json &j) const override {
                        j = {
                            { "replicas", mpi.nproc() },
                            { "datasize", pt.getFormat() }
                        };
                        json &_j = j["exchange"];
                        _j = json::object();
                        for (auto &m : accmap)
                            _j[m.first] = {
                                {"attempts", m.second.cnt},
                                {"acceptance", m.second.avg()}
                            };
                    }

                    void _move(Change &change) override {
                        double Vold = spc.geo.getVolume();
                        findPartner();
                        Tpvec p; // temperary storage
                        p.resize(spc.p.size());
                        if (goodPartner()) {
                            change.all=true;
                            pt.sendExtra[VOLUME]=Vold;  // copy current volume for sending
                            pt.recv(mpi, partner, p); // receive particles
                            pt.send(mpi, spc.p, partner);     // send everything
                            pt.waitrecv();
                            pt.waitsend();

                            double Vnew = pt.recvExtra[VOLUME];
                            if (Vnew<1e-9 || spc.p.size() != p.size())
                                MPI_Abort(mpi.comm, 1);

                            if (std::fabs(Vnew-Vold)>1e-9)
                                change.dV=true; 

                            spc.p = p;
                            spc.geo.setVolume(Vnew);

                            // update mass centers
                            for (auto& g : spc.groups)
                                if (g.atomic==false)
                                    g.cm = Geometry::massCenter(g.begin(), g.end(),
                                            spc.geo.boundaryFunc, -g.begin()->pos);
                        }
                    }

                    double exchangeEnergy(double mydu) {
                        std::vector<MPI::FloatTransmitter::floatp> duSelf(1), duPartner;
                        duSelf[0]=mydu;
                        duPartner = ft.swapf(mpi, duSelf, partner);
                        return duPartner.at(0);               // return partner energy change
                    } //!< Exchange energy with partner

                    double bias(Change &change, double uold, double unew) override {
                        return exchangeEnergy(unew-uold); // Exchange dU with partner (MPI) 
                    }

                    std::string id() {
                        std::ostringstream o;
                        if (mpi.rank() < partner)
                            o << mpi.rank() << " <-> " << partner;
                        else
                            o << partner << " <-> " << mpi.rank();
                        return o.str();
                    } //!< Unique string to identify set of partners

                    void _accept(Change &change) override {
                        if ( goodPartner() )
                            accmap[ id() ] += 1;
                    }

                    void _reject(Change &change) override {
                        if ( goodPartner() )
                            accmap[ id() ] += 0;
                    }

                    void _from_json(const json &j) override {
                        pt.setFormat( j.value("format", std::string("XYZQI") ) );
                    }

                public:
                    ParallelTempering(Tspace &spc, MPI::MPIController &mpi ) : spc(spc), mpi(mpi) {
                        name="temper";
                        partner=-1;
                        pt.recvExtra.resize(1);
                        pt.sendExtra.resize(1);
                    }
            };
#endif

        template<typename Tspace>
            class Propagator : public BasePointerVector<Movebase> {
                private:
                    int _repeat;
                    std::discrete_distribution<> dist;
                    std::vector<double> w; // list of weights for each move

                    void addWeight(double weight=1) {
                        w.push_back(weight);
                        dist = std::discrete_distribution<>(w.begin(), w.end());
                        _repeat = int(std::accumulate(w.begin(), w.end(), 0.0));
                    }

                public:
                    using BasePointerVector<Movebase>::vec;
                    inline Propagator() {}
                    inline Propagator(const json &j, Tspace &spc, MPI::MPIController &mpi) {

                        if (j.count("random")==1)
                            Movebase::slump = j["random"]; // slump is static --> shared for all moves

                        for (auto &m : j.at("moves")) {// loop over move list
                            size_t oldsize = vec.size();
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                try {
#ifdef ENABLE_MPI
                                    if (it.key()=="temper") this->template push_back<Move::ParallelTempering<Tspace>>(spc, mpi);
#endif
                                    if (it.key()=="moltransrot") this->template push_back<Move::TranslateRotate<Tspace>>(spc);
                                    if (it.key()=="transrot") this->template push_back<Move::AtomicTranslateRotate<Tspace>>(spc);
                                    if (it.key()=="pivot") this->template push_back<Move::Pivot<Tspace>>(spc);
                                    if (it.key()=="volume") this->template push_back<Move::VolumeMove<Tspace>>(spc);
                                    if (it.key()=="cluster") this->template push_back<Move::Cluster<Tspace>>(spc);

                                    if (vec.size()==oldsize+1) {
                                        vec.back()->from_json( it.value() );
                                        addWeight(vec.back()->repeat);
                                    } else
                                        std::cerr << "warning: ignoring unknown move '" << it.key() << "'" << endl;
                                } catch (std::exception &e) {
                                    throw std::runtime_error("Error adding move '" + it.key() + "': " + e.what());
                                }
                            }
                        }
                    }

                    int repeat() { return _repeat; }

                    auto sample() {
                        if (!vec.empty()) {
                            assert(w.size() == vec.size());
                            return vec.begin() + dist( Move::Movebase::slump.engine );
                        }
                        return vec.end();
                    } //!< Pick move from a weighted, random distribution
            };

    }//Move namespace

    template<class Tgeometry, class Tparticle>
        class MCSimulation {
            private:
                typedef Space<Tgeometry, Tparticle> Tspace;
                typedef typename Tspace::Tpvec Tpvec;

                bool metropolis(double du) const {
                    if (std::isnan(du))
                        return false;
                    if (du<0)
                        return true;
                    return ( Move::Movebase::slump() > std::exp(-du)) ? false : true;
                } //!< Metropolis criterion (true=accept)

                struct State {
                    Tspace spc;
                    Energy::Hamiltonian<Tspace> pot;
                    State(const json &j) : spc(j), pot(spc,j) {}

                    void sync(State &other, Change &change) {
                        spc.sync( other.spc, change );
                        pot.sync( &other.pot, change );
                    }
                }; //!< Contains everything to describe a state

                State state1, // old state
                      state2; // new state (trial);
                double uinit=0, dusum=0;
                Average<double> uavg;

                void init() {
                    state1.pot.key = Energy::Energybase::OLD; // this is the old energy (current)
                    state2.pot.key = Energy::Energybase::NEW; // this is the new energy (trial)
                    dusum=0;
                    Change c; c.all=true;
                    uinit = state1.pot.energy(c);
                    state2.sync(state1, c);
                    assert(state1.pot.energy(c) == state2.pot.energy(c));
                }

            public:
                Move::Propagator<Tspace> moves;

                auto& pot() { return state1.pot; }
                auto& space() { return state1.spc; }
                const auto& pot() const { return state1.pot; }
                const auto& space() const { return state1.spc; }
                const auto& geometry() const { return state1.spc.geo; }
                const auto& particles() const { return state1.spc.p; }

                double drift() {
                    Change c; c.all=true;
                    double ufinal = state1.pot.energy(c);
                    return ( ufinal-(uinit+dusum) ) / uinit; 
                } //!< Calculates the relative energy drift from initial configuration

                MCSimulation(const json &j, MPI::MPIController &mpi) : state1(j), state2(j), moves(j, state2.spc, mpi) {
                    init();
                }

                void restore(const json &j) {
                    state1.spc = j;
                    state2.spc = j;
                    init();
                } //!< restore system from previously saved state

                void move() {
                    Change change;
                    for (int i=0; i<moves.repeat(); i++) {
                        auto mv = moves.sample(); // pick random move
                        if (mv != moves.end() ) {

                            change.clear();
                            (**mv).move(change);

                            if (!change.empty()) {
                                double unew, uold, du;
#pragma omp parallel sections
                                {
#pragma omp section
                                    { unew = state2.pot.energy(change); }
#pragma omp section
                                    { uold = state1.pot.energy(change); }
                                }
                                du = unew - uold;
                                double bias = (**mv).bias(change, uold, unew);
                                if ( metropolis(du + bias) ) { // accept move
                                    state1.sync( state2, change );
                                    (**mv).accept(change);
                                }
                                else { // reject move
                                    state2.sync( state1, change );
                                    (**mv).reject(change);
                                    du=0;
                                }
                                dusum+=du; // sum of all energy changes
                            }
                        }
                    }
                }

                void to_json(json &j) {
                    j = state1.spc.info();
                    j["temperature"] = pc::temperature / 1.0_K;
                    j["moves"] = moves;
                    j["energy"].push_back(state1.pot);
                }
        };

    template<class Tgeometry, class Tparticle>
        void to_json(json &j, MCSimulation<Tgeometry,Tparticle> &mc) { mc.to_json(j); }

}//namespace
