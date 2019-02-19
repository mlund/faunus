#pragma once

#include "core.h"
#include "energy.h"
#include "average.h"
//#include "analysis.h"
#include "potentials.h"
#include "mpi.h"

namespace Faunus {
    namespace Move {

        class Movebase {
            private:
                virtual void _move(Change&)=0; //!< Perform move and modify change object
                virtual void _accept(Change&); //!< Call after move is accepted
                virtual void _reject(Change&); //!< Call after move is rejected
                virtual void _to_json(json &j) const=0; //!< Extra info for report if needed
                virtual void _from_json(const json &j)=0; //!< Extra info for report if needed
                TimeRelativeOfTotal<std::chrono::microseconds> timer; //!< Timer for whole move
                TimeRelativeOfTotal<std::chrono::microseconds> timer_move; //!< Timer for _move() only
            protected:
                unsigned long cnt=0;
                unsigned long accepted=0;
                unsigned long rejected=0;
            public:
                static Random slump;   //!< Shared for all moves
                std::string name;      //!< Name of move
                std::string cite;      //!< Reference
                int repeat=1;          //!< How many times the move should be repeated per sweep

                void from_json(const json &j);
                void to_json(json &j) const; //!< JSON report w. statistics, output etc.
                void move(Change &change); //!< Perform move and modify given change object
                void accept(Change &c);
                void reject(Change &c);
                virtual double bias(Change&, double uold, double unew); //!< adds extra energy change not captured by the Hamiltonian
                inline virtual ~Movebase() {};
        };

        void from_json(const json &j, Movebase &m); //!< Configure any move via json
        void to_json(json &j, const Movebase &m);

        /**
         * @brief Swap the charge of a single atom
         */
        template<typename Tspace>
            class AtomicSwapCharge : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tparticle Tparticle;
                    Tspace& spc; // Space to operate on
                    int molid=-1;
                    double ln10 = log(10);
                    double pKa, pH;
                    Average<double> msqd; // mean squared displacement
                    double _sqd, _bias; // squared displament
                    std::string molname; // name of molecule to operate on
                    Change::data cdata;

                    void _to_json(json &j) const override {
                        j = {
                            {"pH", pH},
                            {"pka", pKa},
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
                            pH = j.at("pH").get<double>();
                            pKa = j.at("pKa").get<double>();
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
                            auto& g = spc.groups[cdata.index];
                            double oldcharge = p->charge;
                            p->charge = fabs(oldcharge - 1);
                            _sqd = fabs(oldcharge - 1) - oldcharge;
                            change.groups.push_back( cdata ); // add to list of moved groups
                            _bias = _sqd*(pH-pKa)*ln10; // one may add bias here...
                        }
                    }

                    double bias(Change&, double, double) override {
                        return _bias;
                    } //!< adds extra energy change not captured by the Hamiltonian

                    void _accept(Change&) override { msqd += _sqd; }
                    void _reject(Change&) override { msqd += 0; }

                public:
                    AtomicSwapCharge(Tspace &spc) : spc(spc) {
                        name = "swapcharge";
                        repeat = -1; // meaning repeat N times
                        cdata.atoms.resize(1);
                        cdata.internal=true;
                    }
            };

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
                            assertKeys(j, {"molecule", "dir", "repeat"});
                            molname = j.at("molecule");
                            auto it = findName(molecules<Tpvec>, molname);
                            if (it == molecules<Tpvec>.end())
                                throw std::runtime_error("unknown molecule '" + molname + "'");
                            molid = it->id();
                            dir = j.value("dir", Point(1,1,1));
                            if (repeat<0) {
                                auto v = spc.findMolecules(molid, Tspace::ALL );
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
                        auto mollist = spc.findMolecules( molid, Tspace::ALL  ); // all `molid` groups
                        if (size(mollist)>0) {
                            auto git = slump.sample( mollist.begin(), mollist.end() ); // random molecule iterator
                            if (not git->empty()) {
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
                        if (p not_eq spc.p.end()) {
                            double dp = atoms.at(p->id).dp;
                            double dprot = atoms.at(p->id).dprot;
                            auto& g = spc.groups[cdata.index];

                            if (dp>0) { // translate
                                Point oldpos = p->pos;
                                p->pos += ranunit(slump,dir) * dp * slump();

                                spc.geo.boundary(p->pos);
                                _sqd = spc.geo.sqdist(oldpos, p->pos); // squared displacement
                                if (not g.atomic) // recalc mass-center for non-molecular groups
                                    g.cm = Geometry::massCenter(g.begin(), g.end(), spc.geo.getBoundaryFunc(), -g.cm);
                            }

                            if (dprot>0) { // rotate
                                Point u = ranunit(slump);
                                double angle = dprot * (slump()-0.5);
                                Eigen::Quaterniond Q( Eigen::AngleAxisd(angle, u) );
                                p->rotate(Q, Q.toRotationMatrix());
                            }

                            if (dp>0 or dprot>0)
                                change.groups.push_back( cdata ); // add to list of moved groups
                        }
                        // else
                        //    std::cerr << name << ": no atoms found" << std::endl;
                    }

                    void _accept(Change&) override { msqd += _sqd; }
                    void _reject(Change&) override { msqd += 0; }

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
                protected:
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace& spc; // Space to operate on
                    int molid=-1;
                    double dptrans=0;
                    double dprot=0;
                    double _sqd; // squared displacement
                    Average<double> msqd; // mean squared displacement

                    Point dir={1,1,1};

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
                            throw std::runtime_error(name+": " + e.what());
                        }
                    } //!< Configure via json object

                    template<class Tmollist>
                      auto selectMolecule(Tmollist &mollist) {
                        //typename Tspace::Tgvec::iterator it;
                         
                        return slump.sample( mollist.begin(), mollist.end() );
                        //return it;
                     }
                    
                    void _move(Change &change) override {
                        assert(molid>=0);
                        assert(!spc.groups.empty());
                        assert(spc.geo.getVolume()>0);
                        // pick random group from the system matching molecule type
                        // TODO: This can be slow -- implement look-up-table in Space
                        auto mollist = spc.findMolecules( molid, Tspace::ACTIVE ); // list of molecules w. 'molid'
                        if (size(mollist)>0) {
                            auto it = selectMolecule( mollist );

                            if (not it->empty()) {
                                

                                    
                                    
                                assert(it->id==molid);

                                if (dptrans>0) { // translate
                                    Point oldcm = it->cm;
                                    Point dp = ranunit(slump,dir) * dptrans * slump();

                                    it->translate( dp, spc.geo.getBoundaryFunc() );
                                    _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
                                    }

                                if (dprot>0) { // rotate
                                    Point u = ranunit(slump);
                                    double angle = dprot * (slump()-0.5);
                                    Eigen::Quaterniond Q( Eigen::AngleAxisd(angle, u) );
                                    it->rotate(Q, spc.geo.getBoundaryFunc());
                                }

                                if (dptrans>0||dprot>0) { // define changes
                                    Change::data d;
                                    d.index = Faunus::distance( spc.groups.begin(), it ); // integer *index* of moved group
                                    d.all = true; // *all* atoms in group were moved
                                    change.groups.push_back( d ); // add to list of moved groups
                                }
                                assert( spc.geo.sqdist( it->cm,
                                            Geometry::massCenter(it->begin(), it->end(),
                                                spc.geo.getBoundaryFunc(),-it->cm) ) < 1e-6 );
                                    

                            }
                        }
                    }
                
                    void _accept(Change&) override { msqd += _sqd; }
                    void _reject(Change&) override { msqd += 0; }

                public:
                    TranslateRotate(Tspace &spc) : spc(spc) {
                        name = "moltransrot";
                        repeat = -1; // m eaning repeat N times
                    }
            };         

        /**
         * @brief Translate and rotate a molecular group
         */
        template<typename Tspace>
            class BiasedTranslateRotate : public Movebase {
                protected:
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace& spc; // Space to operate on
                    
                    int molid=-1, refid1=-1, refid2=-1, cntAcross=0, cnt;
                    double dptrans=0, dprot=0;
                    double p=1;
                    double apad=0, a, b=0, c=0;
                    double _sqd; // squared displacement
                    Average<double> msqd; // mean squared displacement
                    
                    double cosTheta, theta;
                    double x, y;
                    double coord, coordNew, coordTemp;
                    double randNbr;
                    double countNavg=0;
                    double _bias=0, var = 0, varSum = 0, variance = 0.05, countN [1000000],  cntInner=0, countNtot=0, countNsum=0, countNtotsum=0;

                    Point dir={1,1,1};
                    Point cylAxis={0,0,0};
                    Point origo={0,0,0};
                    Point molV={0,0,0};
                    
                    bool findBias = true;

                    void _to_json(json &j) const override {
                        j = {
                            {"dir", dir}, {"dp", dptrans}, {"dprot", dprot}, 
                            {"p", p}, {"origo", origo}, {"apad", apad}, {"b", b}, {"c", c},
                            {"molid", molid}, {"refid", refid1}, {"refid", refid2}, 
                            {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
                            {"# of water in ellipsoid in last loop", countN[cntAcross]},
                            {"Average # of water in ellipsoid calculated in loop", countNavg},
                            //{"# of water in ellipsoid specified", countN},
                            {"Fraction inner/total moves", cntInner/cnt},
                            {"Number of inner moves", cntInner},
                            {"Total number of attempted moves", cnt},
                            {"molecule", molecules<Tpvec>[molid].name}, {"ref1", atoms[refid1].name}, {"ref2", atoms[refid2].name}
                        };
                        _roundjson(j,3);
                    }

                    void _from_json(const json &j) override {
                        assert(!molecules<Tpvec>.empty());
                        try {
                            std::string molname = j.at("molecule");
                            std::string refname1 = j.at("ref1");
                            std::string refname2 = j.at("ref2");
                            auto it = findName(molecules<Tpvec>, molname);
                            auto ref1 = findName(atoms, refname1);
                            auto ref2 = findName(atoms, refname2);
                            if (it == molecules<Tpvec>.end())
                                throw std::runtime_error("unknown molecule '" + molname + "'");
                            if (ref1 == atoms.end())
                                throw std::runtime_error("unknown reference atom '" + refname1 + "'");
                            if (ref2 == atoms.end())
                                throw std::runtime_error("unknown reference atom '" + refname2 + "'");
                            molid = it->id();
                            refid1 = ref1->id();
                            refid2 = ref2->id();
                            dir = j.value("dir", Point(1,1,1));
                            dprot = j.at("dprot");
                            dptrans = j.at("dp");
                            p = j.at("p");
                            apad = j.at("apad");
                            b = j.at("b");
                            c = j.at("c");
                            variance = j.at("var");
                            //countN = j.at("n");
                            //countNtot = j.at("Ntot");
                            
                            if (repeat<0) {
                                auto v = spc.findMolecules(molid);
                                repeat = std::distance(v.begin(), v.end());
                            }
                        }
                        catch (std::exception &e) {
                            throw std::runtime_error(name+": " + e.what());
                        }
                    } //!< Configure via json object

                    template<class Tmollist>
                      auto selectMolecule(Tmollist &mollist) {
                        //typename Tspace::Tgvec::iterator it;
                         
                        return slump.sample( mollist.begin(), mollist.end() );
                        //return it;
                     }
                    
                    void _move(Change &change) override {
                        assert(molid>=0);
                        assert(!spc.groups.empty());
                        assert(spc.geo.getVolume()>0);
                        _bias = 0.0;
                        // pick random group from the system matching molecule type
                        // TODO: This can be slow -- implement look-up-table in Space
                        auto mollist = spc.findMolecules( molid, Tspace::ACTIVE ); // list of molecules w. 'molid'
                        auto reflist1 = spc.findAtoms( refid1 ); 
                        auto reflist2 = spc.findAtoms( refid2 );
                        if (size(mollist)>0) {
                            auto it = selectMolecule( mollist );
                            auto ref1 = selectMolecule( reflist1 );
                            auto ref2 = selectMolecule( reflist2 );
                            cylAxis = spc.geo.vdist( ref2->pos, ref1->pos )*0.50000;  // half vector between reference points
                            origo = ref2->pos - cylAxis;  // Coordinates of middle point between reference points: new origo
                             // kolla (ref1 + ref2)/2 med pbc och att det blir lika med origo
                            a = cylAxis.norm()+apad;


                            if (not it->empty()) {
                                
                                assert(it->id==molid);
                                
                                randNbr = slump(); 
                                // if ((spc.geo.vdist(it->cm, it0->cm).norm() < dist) || ((spc.geo.vdist(it->cm, it0->cm).norm() > dist) && (p > slump()))) {
                                molV = spc.geo.vdist( it->cm, origo);
                                cosTheta = molV.dot(cylAxis)/molV.norm()/cylAxis.norm();
                                theta = acos(cosTheta);
                                x = cosTheta*molV.norm();
                                y = sin(theta)*molV.norm();
                                coord = x*x/(a*a)+y*y/(b*b);
                                //cntTot++; 

                                if ( not ( coord > 1.00000 && p < randNbr ) ) {
                                    
                                    if (coord <= 1.00000) 
                                        cntInner += 1.00000;
                                    
                                    
                                    cnt +=1;
                                    var = 0.00; 
                                    varSum = 0.00;
                                    countNtot = 0.00;
                                    //cout << "moves across boundary: " << cntacross << "\n"; 
                                    countN[cnt] = 0.00;
                                    for ( auto &g : mollist) {
                                        molV = spc.geo.vdist( g.cm, origo );
                                        cosTheta = molV.dot(cylAxis)/molV.norm()/cylAxis.norm();
                                        theta = acos(cosTheta);
                                        x = cosTheta*molV.norm();
                                        y = sin(theta)*molV.norm();
                                        coordTemp = x*x/(a*a)+y*y/(b*b);
                                        if (coordTemp <= 1.00000) {
                                            countN[cnt] += 1.00;
                                        }
                                    }
                                    countNsum += countN[cnt];
                                    countNavg = countNsum/cnt;
                                    for (int i=1; i<=cnt; i++) {
                                        varSum += (countN[i]-countNavg)*(countN[i]-countNavg);
                                        //cout << "variance sum in loop, not normalized: " << varsum << "\n\n";
                                    }
                                    if (cnt > 500) {
                                        var = varSum/cnt;
                                        //cout << "variance sum, not normalized: " << varsum << "\n";
                                        //cout << "variance: " << var << "\n\n";
                                        /*
                                        if (var <= variance) {
                                             
                                            findBias = false;
                                            cout << "bias found, with average number of particles inside ellipsoid = " << countnavg << " after " << cntacross << " steps across boundary\n\n";
                                            cout << "this steps number of water molecules inside: " << countn[cntacross] << "\n";
                                            cout << "this steps average number of water molecules inside: " << countnavg << "\n";
                                            cout << "this steps bias when going in->out: " << _bias << "\n";
                                            cout << "moves across boundary: " << cntacross << "\n"; 
                                            cout << "variance: " << var << "\n\n";

                                        }
                                        */
                                    }
                                     

                                    if (cnt%50000 == 0) {
                                        cout << "this steps number of water molecules inside: " << countN[cnt] << "\n";
                                        cout << "this steps average number of water molecules inside: " << countNavg << "\n";
                                        cout << "this steps bias: " << _bias << "\n";
                                        cout << "moves: " << cnt << "\n"; 
                                        cout << "variance: " << var << "\n\n";
                                    }

                                    if (dptrans>0) { // translate
                                        Point oldcm = it->cm;
                                        Point dp = ranunit(slump,dir) * dptrans * slump();

                                        it->translate( dp, spc.geo.getBoundaryFunc() );
                                        _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
                                    }

                                    if (dprot>0) { // rotate
                                        Point u = ranunit(slump);
                                        double angle = dprot * (slump()-0.5);
                                        Eigen::Quaterniond Q( Eigen::AngleAxisd(angle, u) );
                                        it->rotate(Q, spc.geo.getBoundaryFunc());
                                    }

                                    if (dptrans>0||dprot>0) { // define changes
                                        Change::data d;
                                        d.index = Faunus::distance( spc.groups.begin(), it ); // integer *index* of moved group
                                        d.all = true; // *all* atoms in group were moved
                                        change.groups.push_back( d ); // add to list of moved groups
                                    }
                                    assert( spc.geo.sqdist( it->cm,
                                                Geometry::massCenter(it->begin(), it->end(),
                                                    spc.geo.getBoundaryFunc(),-it->cm) ) < 1e-6 );
                                      
                                    /*
                                    if (findBias == true) {  
                                        countN[cnt-1] = 0;
                                        countNtot = 0; 
                                        for ( auto &g : mollist) {
                                            molV = spc.geo.vdist( g.cm, origo );
                                            cosTheta = molV.dot(cylAxis)/molV.norm()/cylAxis.norm();
                                            theta = acos(cosTheta);
                                            x = cosTheta*molV.norm();
                                            y = sin(theta)*molV.norm();
                                            coordTemp = x*x/(a*a)+y*y/(b*b);
                                            countNtot++;
                                            if (coordTemp < 1.0) {
                                                countN[cnt-1]++;
                                            }
                                        }
                                        
                                        countNsum += countN[cnt-1];
                                        countNavg = countNsum/cnt;
                                        for (int i=0; i<cnt; i++) {
                                            varSum = (countN[i]-countNavg)*(countN[i]-countNavg);
                                        }
                                        if (cntInner > 1) {
                                            var = varSum/cnt;

                                            if (var < variance) {
                                      
                                                findBias = false;
                                                cout << "Bias found, with average number of particles inside ellipsoid = " << countNavg << "\n\n";
                                            }
                                        } 
                                    }



                                        //countNtotsum += countNtot; 
                                      
                                    
                                    if (countN == 0 && countNtot == 0) {
                                        for ( auto &g : mollist) {
                                            molV = spc.geo.vdist( g.cm, origo );
                                            cosTheta = molV.dot(cylAxis)/molV.norm()/cylAxis.norm();
                                            theta = acos(cosTheta);
                                            x = cosTheta*molV.norm();
                                            y = sin(theta)*molV.norm();
                                            coordTemp = x*x/(a*a)+y*y/(b*b);
                                            countNtot++;
                                            if (coordTemp < 1.0) {
                                                countN++;
                                            }
                                        }
                                        countNsum += countN;
                                        countNtotsum += countNtot; 
                                        cout << "Particles within ellipsoid: " << countN << "\n";
                                        cout << "Sum of particles within ellipsoid: " << countNsum << "\n";
                                        cout << "Sum of particles in box: " << countNtotsum << "\n\n";
                                        countN = 0;
                                        countNtot = 0;

                                    }
                                    
                                    if (cnt == 1) {
                                        for ( auto &g : mollist) {
                                            molV = spc.geo.vdist( g.cm, origo );
                                            cosTheta = molV.dot(cylAxis)/molV.norm()/cylAxis.norm();
                                            theta = acos(cosTheta);
                                            x = cosTheta*molV.norm();
                                            y = sin(theta)*molV.norm();
                                            coordTemp = x*x/(a*a)+y*y/(b*b);
                                            initNtot++;
                                            if (coordTemp < 1.0) {
                                                initN++;
                                            }
                                        }    
                                    }
                                    */
                                    
                                    molV = spc.geo.vdist( it->cm, origo);
                                    cosTheta = molV.dot(cylAxis)/molV.norm()/cylAxis.norm();
                                    theta = acos(cosTheta);
                                    x = cosTheta*molV.norm();
                                    y = sin(theta)*molV.norm();
                                    coordNew = x*x/(a*a)+y*y/(b*b);

                                    if (coord <= 1.00000 && coordNew > 1.00000) {
                                        _bias = -log(p/(1-(1-p)/(p*countNtot+(1-p)*countN[cnt])));
                                    }
                                    else if (coord > 1.0000 && coordNew <= 1.00000) {
                                        _bias = -log(1/(1+(1-p)/(p*countNtot+(1-p)*countN[cnt])));
                                    }

                                    /*
                                    if (coord <= 1.00000 && coordNew > 1.00000) {
                                        if (findBias == true) {
                                            varsum = 0.00;
                                            countntot = 0.00;
                                            cntacross += 1.00;
                                            //cout << "moves across boundary: " << cntacross << "\n"; 
                                            countn[cntacross] = 0.00;
                                            for ( auto &g : mollist) {
                                                molv = spc.geo.vdist( g.cm, origo );
                                                costheta = molv.dot(cylaxis)/molv.norm()/cylaxis.norm();
                                                theta = acos(costheta);
                                                x = costheta*molv.norm();
                                                y = sin(theta)*molv.norm();
                                                coordtemp = x*x/(a*a)+y*y/(b*b);
                                                if (coordtemp <= 1.00000) {
                                                    countn[cntacross] += 1.00;
                                                }
                                            }
                                            countnsum += countn[cntacross];
                                            countnavg = countnsum/cntacross;
                                            for (int i=1; i<=cntacross; i++) {
                                                varsum += (countn[i]-countnavg)*(countn[i]-countnavg);
                                                //cout << "variance sum in loop, not normalized: " << varsum << "\n\n";
                                            }
                                            if (cntacross > 4) {
                                                var = varsum/cntacross;
                                                //cout << "variance sum, not normalized: " << varsum << "\n";
                                                //cout << "variance: " << var << "\n\n";
                                                if (var <= variance) {
                                                     
                                                    findbias = false;
                                                    cout << "bias found, with average number of particles inside ellipsoid = " << countnavg << " after " << cntacross << " steps across boundary\n\n";
                                                    cout << "this steps number of water molecules inside: " << countn[cntacross] << "\n";
                                                    cout << "this steps average number of water molecules inside: " << countnavg << "\n";
                                                    cout << "this steps bias when going in->out: " << _bias << "\n";
                                                    cout << "moves across boundary: " << cntacross << "\n"; 
                                                    cout << "variance: " << var << "\n\n";

                                                }
                                            }
                                            _bias = -log(p/(1-(1-p)/(p*countntot+(1-p)*countn[cntacross])));
                                            if (cntacross%1000) {
                                                cout << "this steps number of water molecules inside: " << countn[cntacross] << "\n";
                                                cout << "this steps average number of water molecules inside: " << countnavg << "\n";
                                                cout << "this steps bias when going in->out: " << _bias << "\n";
                                                cout << "moves across boundary: " << cntacross << "\n"; 
                                                cout << "variance: " << var << "\n\n";
                                            }
                                        }
                                        else {
                                            _bias = -log(p/(1-(1-p)/(p*countNtot+(1-p)*countNavg)));
                                            //cout << "Average number of water molecules inside: " << countNavg << "\n";
                                            //cout << "Bias when going in->out: " << _bias << "\n\n";
                                        }
                                    }
                                    else if (coord > 1.00000 && coordNew <= 1.00000) {
                                        if (findBias == true) {
                                            varSum = 0.00;
                                            countNtot = 0.00;
                                            cntAcross += 1.00;
                                            countN[cntAcross] = 0.00;
                                            for ( auto &g : mollist) {
                                                molV = spc.geo.vdist( g.cm, origo );
                                                cosTheta = molV.dot(cylAxis)/molV.norm()/cylAxis.norm();
                                                theta = acos(cosTheta);
                                                x = cosTheta*molV.norm();
                                                y = sin(theta)*molV.norm();
                                                coordTemp = x*x/(a*a)+y*y/(b*b);
                                                if (coordTemp <= 1.00000) {
                                                    countN[cntAcross] += 1.00;
                                                }
                                            }
                                            countNsum += countN[cntAcross];
                                            countNavg = countNsum/cntAcross;
                                            for (int i=1; i<=cntAcross; i++) {
                                                varSum += (countN[i]-countNavg)*(countN[i]-countNavg);
                                                //cout << "Variance sum in loop, not normalized: " << varSum << "\n\n";
                                            }
                                            if (cntAcross > 4) {
                                                var = varSum/cntAcross;
                                                //cout << "Variance sum, not normalized: " << varSum << "\n";
                                                if (var < variance) {
                                      
                                                    findBias = false;
                                                    cout << "Bias found, with average number of particles inside ellipsoid = " << countNavg << " after " << cntAcross << " steps across boundary\n\n";
                                                    cout << "This steps number of water molecules inside: " << countN[cntAcross] << "\n";
                                                    cout << "This steps average number of water molecules inside: " << countNavg << "\n";
                                                    cout << "This steps bias when going in->out: " << _bias << "\n\n";
                                                    cout << "Moves across boundary: " << cntAcross << "\n"; 
                                                    cout << "Variance: " << var << "\n\n";

                                                }
                                            }
                                            _bias = -log(1/(p*(1+(1-p)/(p*countNtot+(1-p)*countN[cntAcross]))));
                                            if (cntAcross%1000) {
                                                cout << "This steps number of water molecules inside: " << countN[cntAcross] << "\n";
                                                cout << "This steps average number of water molecules inside: " << countNavg << "\n";
                                                cout << "This steps bias when going out->in: " << _bias << "\n\n";
                                                cout << "Moves across boundary: " << cntAcross << "\n"; 
                                                cout << "Variance: " << var << "\n\n";
                                            }
                                        }                                            
                                        else { 
                                            _bias = -log(1/(p*(1+(1-p)/(p*countNtot+(1-p)*countNavg))));
                                            //cout << "Average number of water molecules inside: " << countNavg << "\n";
                                            //cout << "Bias when going out->in: " << _bias << "\n\n";
                                        }

                                    }
                                    else {
                                        _bias = 0.00;
                                        countN[cntAcross] = 0.00;
                                        for ( auto &g : mollist) {
                                            molV = spc.geo.vdist( g.cm, origo );
                                            cosTheta = molV.dot(cylAxis)/molV.norm()/cylAxis.norm();
                                            theta = acos(cosTheta);
                                            x = cosTheta*molV.norm();
                                            y = sin(theta)*molV.norm();
                                            coordTemp = x*x/(a*a)+y*y/(b*b);
                                            if (coordTemp <= 1.00000) {
                                                countN[cntAcross] += 1.00;
                                            }
                                        }
                                        //cout << "This steps number of water molecules inside: " << countN[cntAcross] << "\n";
                                        //cout << "This steps average number of water molecules inside: " << countNavg << "\n";
                                        //cout << "This steps bias when going out->in: " << _bias << "\n\n";
                                        //cout << "Moves across boundary: " << cntAcross << "\n"; 
                                        //cout << "Variance: " << var << "\n\n";
                                    }
                                    */     
                                }
                            }
                        }
                    }
                    double bias(Change&, double, double) override { 
                        return _bias;
                    }

                    void _accept(Change&) override { 
                        msqd += _sqd; 
                        /*
                        cout << "This steps number of water molecules inside: " << countN[cntAcross] << "\n";
                        cout << "This steps average number of water molecules inside: " << countNavg << "\n";
                        cout << "This steps bias when going in->out: " << _bias << "\n";
                        cout << "Moves across boundary: " << cntAcross << "\n"; 
                        cout << "Variance: " << var << "\n\n";
                        */
                    }
                    void _reject(Change&) override { msqd += 0; }

                public:
                    BiasedTranslateRotate(Tspace &spc) : spc(spc) {
                        name = "biasedmoltransrot";
                        repeat = -1; // m eaning repeat N times
                    }
            };

        /**
         * @brief Move that will swap conformation of a molecule
         *
         * This will swap between different molecular conformations
         * as defined in `MoleculeData` with `traj` and `weight`.
         * If defined, the weight
         * distribution is respected, otherwise all conformations
         * have equal intrinsic weight. Upon insertion, the new conformation
         * is randomly oriented and placed on top of the mass-center of
         * an exising molecule. That is, there is no mass center movement.
         *
         * @todo Add feature to align molecule on top of an exiting one
         * @todo Expand `_info()` to show number of conformations
         * @warning Weighted distributions untested and not verified for correctness
         * @date Malmo, November 2016
         */
        template<class Tspace>
            class ConformationSwap : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef MoleculeData<Tpvec> Tmoldata;
                    RandomInserter<Tmoldata> inserter;
                    Tspace& spc; // Space to operate on
                    int molid=-1;
                    int newconfid=-1;

                    void _to_json(json &j) const override {
                        j = {
                            {"molid", molid},
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
                            if ( molecules<Tpvec>[molid].conformations.size()<2)
                                throw std::runtime_error("minimum two conformations required");
                            if (repeat<0) {
                                auto v = spc.findMolecules(molid);
                                repeat = std::distance(v.begin(), v.end());
                            }
                        }
                        catch (std::exception &e) {
                            throw std::runtime_error(name+": " + e.what());
                        }
                    } //!< Configure via json object

                    void _move(Change &change) override {
                        assert(molid>=0);
                        assert(change.empty());

                        auto mollist = spc.findMolecules( molid, Tspace::ACTIVE ); // list of molecules w. 'molid'
                        if ( size(mollist)>0 ) {
                            auto g = slump.sample( mollist.begin(), mollist.end() );
                            if (not g->empty()) {
                                inserter.offset = g->cm;

                                // Get a new conformation that should be properly wrapped around the boundaries
                                // (if applicable) and have the same mass-center as "g->cm".
                                Tpvec p = inserter(spc.geo, spc.p, molecules<Tpvec>[molid]);
                                if (p.size() not_eq g->size())
                                    throw std::runtime_error(name + ": conformation atom count mismatch");

                                newconfid = molecules<Tpvec>[molid].conformations.index;

                                std::copy( p.begin(), p.end(), g->begin() ); // override w. new conformation
#ifndef NDEBUG
                                // this move shouldn't move mass centers, so let's check if this is true:
                                Point newcm = Geometry::massCenter(p.begin(), p.end(), spc.geo.getBoundaryFunc(), -g->cm);
                                if ( (newcm - g->cm).norm()>1e-6 )
                                    throw std::runtime_error(name + ": unexpected mass center movement");
#endif
                                Change::data d;
                                d.index = Faunus::distance(spc.groups.begin(), g); // integer *index* of moved group
                                d.all = true; // *all* atoms in group were moved
                                d.internal = false; // we *don't* want to calculate the internal energy
                                change.groups.push_back( d ); // add to list of moved groups
                            }
                        }
                    }

                    void _accept(Change &change) override {
                        assert(change.groups.size()==1);
                        spc.groups[ change.groups.front().index ].confid = newconfid;
                    }

                public:

                    ConformationSwap(Tspace &spc) : spc(spc) {
                        name = "conformationswap";
                        repeat = -1; // meaning repeat n times
                        inserter.dir = {0,0,0};
                        inserter.rotate = true;
                    }

            }; // end of conformation swap move


        /**
         * @brief Sketch for MD move
         */
        template<typename Tspace>
            class ForceMove : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    void _to_json(json&) const {};
                    void _from_json(const json&) {};
                    std::vector<Point> forces, velocities;
                public:
                    ForceMove() {
                        // resize forces and velocities to mathc spc.p
                    }
            }; // end of forcemove

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] TranslateRotate")
        {
            typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
            typedef Space<Geometry::Chameleon, Tparticle> Tspace;
            typedef typename Tspace::Tpvec Tpvec;

            CHECK( !atoms.empty() ); // set in a previous test
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
                    Average<double> msqd, Vavg; // mean squared displacement
                    double dV=0, deltaV=0, Vnew=0, Vold=0;

                    void _to_json(json &j) const override {
                        using namespace u8;
                        if (cnt>0) {
                            j = {
                                {"dV", dV}, {"method", method->first},
                                {bracket("V"), Vavg.avg()},
                                {rootof + bracket(Delta + "V" + squared), std::sqrt(msqd.avg())},
                                {cuberoot + rootof + bracket(Delta + "V" + squared),
                                    std::cbrt(std::sqrt(msqd.avg()))}
                            };
                            _roundjson(j,3);
                        }
                    }

                    void _from_json(const json &j) override {
                        try {
                            method = methods.find( j.value("method", "isotropic") );
                            if (method==methods.end())
                                std::runtime_error("unknown volume change method");
                            dV = j.at("dV");
                        } catch (std::exception &e) {
                            throw std::runtime_error(e.what());
                        }
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

                    void _accept(Change&) override {
                        msqd += deltaV*deltaV;
                        Vavg += spc.geo.getVolume();
                    }
                    void _reject(Change&) override {
                        msqd += 0;
                        Vavg += spc.geo.getVolume();
                    }

                public:
                    VolumeMove(Tspace &spc) : spc(spc) {
                        name = "volume";
                        repeat = 1;
                    }
            }; // end of VolumeMove

        /**
         * @brief Displaces charge on a single atom
         */
        template<typename Tspace>
            class ChargeMove : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace& spc; // Space to operate on
                    Average<double> msqd; // mean squared displacement
                    double dq=0, deltaq=0;
                    int atomIndex;
                    Change::data cdata;

                    void _to_json(json &j) const override {
                        using namespace u8;
                        j = {
                            {"index", atomIndex},
                            {"dq", dq},
                            {rootof + bracket(Delta + "q" + squared), std::sqrt(msqd.avg())},
                            {cuberoot + rootof + bracket(Delta + "q" + squared),
                                std::cbrt(std::sqrt(msqd.avg()))}
                        };
                        _roundjson(j,3);
                    }

                    void _from_json(const json &j) override {
                        dq = j.at("dq").get<double>();
                        atomIndex = j.at("index").get<int>();
                        auto git = spc.findGroupContaining( spc.p[atomIndex] ); // group containing atomIndex
                        cdata.index = std::distance( spc.groups.begin(), git ); // integer *index* of moved group
                        cdata.atoms[0] = std::distance(git->begin(), spc.p.begin()+atomIndex );  // index of particle rel. to group
                    }

                    void _move(Change &change) override {
                        if (dq>0) {
                            auto &p = spc.p[atomIndex]; // refence to particle
                            double qold = p.charge;
                            p.charge +=  dq * (slump()-0.5);
                            deltaq = p.charge-qold;
                            change.groups.push_back( cdata ); // add to list of moved groups
                        } else deltaq=0;
                    }

                    void _accept(Change&) override { msqd += deltaq*deltaq; }
                    void _reject(Change&) override { msqd += 0; }

                public:
                    ChargeMove(Tspace &spc) : spc(spc) {
                        name = "charge";
                        repeat = 1;
                        cdata.internal=true; // the group is internally changed
                        cdata.atoms.resize(1); // we change exactly one atom
                    }
            };

        /*
         * @brief Establishes equilibrium of matter
         * Establishes equilibrium of matter between all species
         *
         * Consider the dissociation process AX=A+X. This class will locate
         * all species of type AX and A and make a MC swap move between them.
         * X can be implicit, meaning that it enters only with its chemical potential
         * (activity). The reacting species, the equilibrium constant,
         * and the activities are read from the JSON input file.
         *
         */
        template<typename Tspace>
            class SpeciationMove : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tparticle Tparticle;

                    Tspace& spc;
                    Tspace *otherspc;
                    ReactionData<Tpvec> *trialprocess;
                    std::map<std::string, Average<double>> accmap;

                    double lnK;
                    bool forward;
                    std::vector<int> molDel;                  // index of groups to delete
                    std::vector<int> atomDel;                 // atom index to delete
                    std::map<int, int> molcnt_ins, atomcnt_ins,
                        molcnt_del, atomcnt_del,
                        molcnt, atomcnt;   // id's and number of inserted/deleted mols and atoms
                    std::multimap<int, Tpvec> pmap;      // coordinates of mols and atoms to be inserted
                    unsigned int Ndeleted, Ninserted;    // Number of accepted deletions and insertions

                    void _to_json(json &j) const override {
                        j = {
                            // { "replicas", mpi.nproc() },
                            // { "datasize", pt.getFormat() }
                        };
                        json &_j = j["reactions"];
                        _j = json::object();
                        for (auto &m : accmap)
                            _j[m.first] = {
                                {"attempts", m.second.cnt},
                                {"acceptance", m.second.avg()}
                            };
                    }

                    void _from_json(const json&) override {};

                public:

                    SpeciationMove(Tspace &spc) : spc(spc) {
                        name = "speciation";
                        repeat = 1;
                    }

                    void setOther(Tspace &ospc) {
                        otherspc = &ospc;
                    }

                    double energy(); //!< Returns intrinsic energy of the process

                    void _move(Change &change) override {
                        if ( reactions<Tpvec>.size()>0 ) {
                            auto rit = slump.sample( reactions<Tpvec>.begin(), reactions<Tpvec>.end() );
                            lnK = rit->lnK;
                            forward = (bool)slump.range(0,1); // random boolean
                            trialprocess = &(*rit);
                            if ( rit->empty(forward) ) // Enforce canonic constraint if invoked
                                return; // Out of material, slip out the back door


                            for (auto &m : rit->Molecules2Add( !forward )) { // Delete checks
                                auto mollist = spc.findMolecules( m.first, Tspace::ALL);
                                if ( molecules<Tpvec>[m.first].atomic ) {
                                    if( size(mollist)!=1 ) // There can be only one
                                        throw std::runtime_error("Bad definition: One group per atomic molecule!");
                                    auto git = mollist.begin();
                                    if ( git->size() < m.second ) // Assure that there are atoms enough in the group
                                        return; // Slip out the back door
                                } else {
                                    mollist = spc.findMolecules( m.first, Tspace::ACTIVE);
                                    if ( size(mollist) <  m.second ) 
                                        return; // Not possible to perform change, escape through the back door
                                }
                            }
                            for (auto &m : rit->Molecules2Add( forward )) { // Addition checks
                                auto mollist = spc.findMolecules( m.first, Tspace::ALL);
                                if ( molecules<Tpvec>[m.first].atomic ) {
                                    if( size(mollist)!=1 ) // There can be only one
                                        throw std::runtime_error("Bad definition: One group per atomic molecule!");
                                    auto git = mollist.begin();
                                    if ( (git->size() + m.second) > git->capacity() ) {  // Assure that there are atoms enough in the group
                                        return; // Slip out the back door
                                    }
                                } else {
                                    mollist = spc.findMolecules( m.first, Tspace::INACTIVE);
                                    if ( size(mollist) <  m.second ) {
                                        return; // Not possible to perform change, escape through the back door
                                    }
                                }
                            }

                            if (rit->swap==true) {
                                auto m1 = rit->Atoms2Add( not forward ); // Swap checks
                                auto m2 = rit->Atoms2Add( forward );
                                assert( (m1.size() == 1) and (m2.size() == 1) && "Bad definition: Only 2 explicit atoms per reaction!"); // Swap A = B
                                auto atomlist = spc.findAtoms( m1.begin()->first );
                                if ( size(atomlist) < 1 ) // Make sure that there are any active atoms to swap
                                    return; // Slip out the back door
                                auto ait = slump.sample( atomlist.begin(), atomlist.end() ); // Random particle iterator
                                auto git = spc.findGroupContaining( *ait );
                                Change::data d;
                                d.atoms.push_back( Faunus::distance(git->begin(), ait) ); // Index of particle rel. to group
                                d.index = Faunus::distance( spc.groups.begin(), git );
                                d.internal = true;
                                d.dNswap = true;
                                change.groups.push_back( d ); // Add to list of moved groups
                                Tparticle p = atoms.at( m2.begin()->first );
                                p.pos = ait->pos;
                                *ait = p;
                                assert(ait->id == m2.begin()->first);
                            }

                            change.dN=true; // Attempting to change the number of atoms / molecules
                            for (auto &m : rit->Molecules2Add( not forward )) { // Delete
                                auto mollist = spc.findMolecules( m.first, Tspace::ALL);
                                if ( molecules<Tpvec>[m.first].atomic ) {
                                    auto git = mollist.begin();
                                    auto othermollist = otherspc->findMolecules(m.first, Tspace::ALL);  // implies that new and old are in sync
                                    auto othergit=othermollist.begin();
                                    Change::data d;
                                    d.index = Faunus::distance( spc.groups.begin(), git ); // integer *index* of moved group
                                    d.internal = true;
                                    d.dNatomic = true;
                                    for ( int N=0; N<m.second; N++ ) {  // deactivate m.second m.first atoms
                                        auto ait = slump.sample( git->begin(), git->end()); // iterator to random atom
                                        // Shuffle back to end, both in trial and new
                                        auto nait = git->end()-1; //iterator to last atom
                                        int dist = Faunus::distance( ait, git->end() ); // distance to random atom from end
                                        if ( Faunus::distance( ait, nait) > 1 ) {
                                            std::iter_swap(ait, nait);
                                            std::iter_swap(othergit->end()-dist-N, othergit->end() - (1+N) );
                                        }
                                        d.atoms.push_back ( Faunus::distance(git->begin(), nait) );
                                        git->deactivate( nait, git->end());
                                    }
                                    std::sort( d.atoms.begin(), d.atoms.end() );
                                    change.groups.push_back( d ); // add to list of moved groups
                                } else {
                                    mollist = spc.findMolecules( m.first, Tspace::ACTIVE);
                                    for ( int N=0; N <m.second; N++ ) {
                                        auto git = slump.sample(mollist.begin(), mollist.end());
                                        git->deactivate( git->begin(), git->end());
                                        Change::data d;
                                        d.index = Faunus::distance( spc.groups.begin(), git ); // integer *index* of moved group
                                        d.all = true; // *all* atoms in group were moved
                                        d.internal = true;
                                        for (int i=0; i<git->capacity(); i++)
                                            d.atoms.push_back(i);
                                        change.groups.push_back( d ); // add to list of moved groups
                                    }
                                }
                            }

                            for (auto &m : rit->Molecules2Add( forward )) { // Add
                                auto mollist = spc.findMolecules( m.first, Tspace::ALL);
                                if ( molecules<Tpvec>[m.first].atomic ) {
                                    auto git = mollist.begin();
                                    Change::data d;
                                    d.index = Faunus::distance( spc.groups.begin(), git);
                                    d.internal = true;
                                    d.dNatomic = true;
                                    for ( int N=0; N<m.second; N++ ) {  // Activate m.second m.first atoms
                                        git->activate( git->end(), git->end() + 1);
                                        auto ait = git->end()-1;
                                        spc.geo.randompos(ait->pos, slump);
                                        spc.geo.getBoundaryFunc()(ait->pos);
                                        d.atoms.push_back( Faunus::distance(git->begin(), ait) );  // Index of particle rel. to group
                                    }
                                    std::sort( d.atoms.begin(), d.atoms.end());
                                    change.groups.push_back( d ); // Add to list of moved groups
                                } else {
                                    mollist = spc.findMolecules( m.first, Tspace::INACTIVE);
                                    for ( int N=0; N <m.second; N++ ) {
                                        auto git = slump.sample(mollist.begin(), mollist.end());
                                        git->activate( git->inactive().begin(), git->inactive().end());
                                        Point newpoint; // = git->cm;
                                        spc.geo.randompos(newpoint, random);
                                        git->translate( -git->cm, spc.geo.getBoundaryFunc() );
                                        git->translate( newpoint, spc.geo.getBoundaryFunc() );
                                        Point u = ranunit(slump);
                                        Eigen::Quaterniond Q( Eigen::AngleAxisd(2*pc::pi*random(), u) );
                                        git->rotate(Q, spc.geo.getBoundaryFunc());
                                        Change::data d;
                                        d.index = Faunus::distance( spc.groups.begin(), git ); // Integer *index* of moved group
                                        d.all = true; // All atoms in group were moved
                                        d.internal = true;
                                        for (int i=0; i<git->capacity(); i++)
                                            d.atoms.push_back(i);
                                        change.groups.push_back( d ); // Add to list of moved groups
                                        assert( spc.geo.sqdist( git->cm, Geometry::massCenter(git->begin(),git->end(),spc.geo.getBoundaryFunc(), -git->cm ) ) < 1e-9 );
                                    }

                                }
                            }

                            std::sort(change.groups.begin(), change.groups.end() );
                        } else
                            throw std::runtime_error("No reactions in list, disable speciation or add reactions");
                    }

                    double bias(Change&, double, double) override {
                        if (forward)
                            return -lnK;
                        return lnK;
                    } //!< adds extra energy change not captured by the Hamiltonian

                    void _accept(Change&) override {
                        accmap[ trialprocess->name ] += 1;
                        trialprocess->N_reservoir += (forward == true) ? -1 : 1;
                        if( trialprocess->N_reservoir < 0 && trialprocess->canonic == true )
                            throw std::runtime_error("There are no negative number of molecules");
                    }

                    void _reject(Change&) override {
                        accmap[ trialprocess->name ] += 0;
                    }

            }; // End of class SpeciationMove

        /**
         * @brief QuadrantJump translates a molecule to another quadrant
         * considering as the origin the center of the box or the center of mass 
         * of a range of atomic indexes specified by "index": [start:stop].
         */
        template<typename Tspace>
            class QuadrantJump : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tparticle Tparticle;
                    Tspace& spc; // Space to operate on
                    int molid=-1;
                    Point dir={1,1,1};
                    std::vector<size_t> index;
                    double _sqd; // squared displacement
                    Average<double> msqd; // mean squared displacement

                    void _to_json(json &j) const override {
                        j = {
                            {"dir", dir},
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
                            index = j.value("index", decltype(index)());
                            if (repeat<0) {
                                auto v = spc.findMolecules(molid);
                                repeat = std::distance(v.begin(), v.end());
                            }
                        }
                        catch (std::exception &e) {
                            throw std::runtime_error(name+": " + e.what());
                        }
                    } //!< Configure via json object


                    void _move(Change &change) override {
                        assert(molid>=0);
                        assert(!spc.groups.empty());
                        assert(spc.geo.getVolume()>0);

                        // pick random group from the system matching molecule type
                        // TODO: This can be slow -- implement look-up-table in Space
                        auto mollist = spc.findMolecules( molid, Tspace::ACTIVE ); // list of molecules w. 'molid'
                        if (size(mollist)>0) {
                            auto it = slump.sample( mollist.begin(), mollist.end() );
                            if (not it->empty()) {
                                assert(it->id==molid);
                                Point oldcm = it->cm;
                                if (index.size()==2) {
                                    auto cm_O = Geometry::massCenter(spc.p.begin()+index[0], spc.p.begin()+index[1]+1, spc.geo.getBoundaryFunc() );
                                    it->translate( -2*spc.geo.vdist(oldcm, cm_O).cwiseProduct(dir.cast<double>()), spc.geo.getBoundaryFunc() );
                                } else {
                                    it->translate( -2*oldcm.cwiseProduct(dir.cast<double>()), spc.geo.getBoundaryFunc() );
                                }
                                _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement
                                Change::data d;
                                d.index = Faunus::distance( spc.groups.begin(), it ); // integer *index* of moved group
                                d.all = true; // *all* atoms in group were moved
                                change.groups.push_back( d ); // add to list of moved groups

                                assert( spc.geo.sqdist( it->cm,
                                            Geometry::massCenter(it->begin(),it->end(),spc.geo.getBoundaryFunc(),-it->cm) ) < 1e-9 );
                            }
                        }
                        else std::cerr << name << ": no molecules found" << std::endl;
                    }

                    void _accept(Change&) override { msqd += _sqd; }
                    void _reject(Change&) override { msqd += 0; }


                public:
                    QuadrantJump(Tspace &spc) : spc(spc) {
                        name = "quadrantjump";
                        repeat = -1; // meaning repeat N times
                    }
            };

        /**
         * @brief Molecular cluster move
         *
         * @todo fix so that it works w. GC molecules (index are calculating before simulation)
         */
        template<typename Tspace>
            class Cluster : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tgroup Tgroup;
                    Tspace& spc;
                    Average<double> msqd, msqd_angle, N;
                    double thresholdsq=0, dptrans=0, dprot=0, angle=0, _bias=0;
                    size_t bias_rejected=0;
                    bool rotate; // true if cluster should be rotated
                    Point dir={1,1,1}, dp;
                    std::vector<std::string> names; // names of molecules to be considered
                    std::vector<int> ids; // molecule id's of molecules to be considered
                    std::vector<size_t> index; // index of all possible molecules to be considered

                    virtual double clusterProbability(const Tgroup &g1, const Tgroup &g2) const {
                        if (spc.geo.sqdist(g1.cm, g2.cm)<=thresholdsq)
                            return 1.0;
                        return 0.0;
                    }

                    void _to_json(json &j) const override {
                        using namespace u8;
                        j = {
                            {"threshold", std::sqrt(thresholdsq)}, {"dir", dir}, {"dp", dptrans}, {"dprot", dprot},
                            {rootof + bracket("r" + squared), std::sqrt(msqd.avg())},
                            {rootof + bracket(theta + squared) + "/" + degrees, std::sqrt(msqd_angle.avg()) / 1.0_deg},
                            {bracket("N"), N.avg()},
                            {"bias rejection rate", double(bias_rejected) / cnt}
                        };
                        _roundjson(j,3);
                    }

                    void _from_json(const json &j) override {
                        assertKeys(j, {"dp", "dprot", "dir", "threshold", "molecules", "repeat"});
                        dptrans = j.at("dp");
                        dir = j.value("dir", Point(1,1,1));
                        dprot = j.at("dprot");
                        thresholdsq = std::pow(j.at("threshold").get<double>(), 2);
                        names = j.at("molecules").get<decltype(names)>(); // molecule names
                        ids = names2ids(molecules<Tpvec>, names);     // names --> molids
                        index.clear();
                        for (auto &g : spc.groups) // loop over all groups
                            if (not g.atomic) // only molecular groups
                                if (g.size()==g.capacity()) // only active particles
                                    if (std::find(ids.begin(), ids.end(), g.id)!=ids.end() )
                                        index.push_back( &g-&spc.groups.front() );
                        if (repeat<0)
                            repeat = index.size();
                    }

                    /**
                     * @param spc Space
                     * @param first Index of initial molecule (randomly selected)
                     * @param index w. all molecules clustered around first (first included)
                     */
                    void findCluster(Tspace &spc, size_t first, std::set<size_t>& cluster) {
                        assert(first < spc.p.size());
                        std::set<size_t> pool(index.begin(), index.end());
                        assert(pool.count(first)>0);

                        cluster.clear();
                        cluster.insert(first);
                        pool.erase(first);

                        size_t n;
                        do { // find cluster (not very clever...)
start:
                            n = cluster.size();
                            for (size_t i : cluster)
                                if (not spc.groups.at(i).empty()) // check if group is inactive
                                    for (size_t j : pool)
                                        if (i!=j)
                                            if (not spc.groups.at(j).empty()) { // check if group is inactive
                                                // probability to cluster
                                                double P = clusterProbability(spc.groups.at(i), spc.groups.at(j));
                                                if ( Movebase::slump() <= P ) {
                                                    cluster.insert(j);
                                                    pool.erase(j);
                                                    goto start; // wow, first goto ever!
                                                }
                                            }
                        } while (cluster.size() != n);

                        // check if cluster is too large
                        double max = spc.geo.getLength().minCoeff()/2;
                        for (auto i : cluster)
                            for (auto j : cluster)
                                if (j>i)
                                    if (spc.geo.sqdist(spc.groups.at(i).cm, spc.groups.at(j).cm)>=max*max)
                                        rotate=false; // skip rotation if cluster larger than half the box length
                    }

                    void _move(Change &change) override {
                        _bias=0;
                        rotate=true;
                        if (thresholdsq>0 and not index.empty()) {
                            std::set<size_t> cluster; // all group index in cluster
                            size_t first = *slump.sample(index.begin(), index.end()); // random molecule (nuclei)
                            findCluster(spc, first, cluster); // find cluster around first

                            N += cluster.size(); // average cluster size
                            Change::data d;
                            d.all=true;
                            dp = ranunit(slump, dir) * dptrans * slump();

                            if (rotate)
                                angle = dprot * (slump()-0.5);
                            else
                                angle = 0;

                            auto boundary = spc.geo.getBoundaryFunc();

                            // lambda function to calculate cluster COM
                            auto clusterCOM = [&]() {
                                double m, sum_m=0;
                                Point cm(0,0,0);
                                Point O = spc.groups[*cluster.begin()].cm;
                                for (auto i : cluster) {
                                    auto &g = spc.groups[i];
                                    Point t = g.cm-O;
                                    boundary(t);
                                    m = g.mass();
                                    cm += m*t;
                                    sum_m += m;
                                }
                                cm = cm/sum_m + O;
                                boundary(cm);
                                return cm;
                            };

                            Point COM = clusterCOM(); // org. cluster center
                            Eigen::Quaterniond Q;
                            Q = Eigen::AngleAxisd(angle, ranunit(slump)); // quaternion

                            for (auto i : cluster) { // loop over molecules in cluster
                                auto &g = spc.groups[i];
                                if (rotate) {
                                    Geometry::rotate(g.begin(), g.end(), Q, boundary, -COM);
                                    g.cm = g.cm-COM;
                                    boundary(g.cm);
                                    g.cm = Q*g.cm+COM;
                                    boundary(g.cm);
                                }
                                g.translate(dp, boundary);
                                d.index=i;
                                change.groups.push_back(d);
                            }

                            // Reject if cluster composition changes during move
                            // Note: this only works for the binary 0/1 probability function
                            // currently implemented in `findCluster()`.

                            std::set<size_t> aftercluster; // all group index in cluster _after_move
                            findCluster(spc, first, aftercluster); // find cluster around first
                            if (aftercluster == cluster)
                                _bias = 0;
                            else {
                                _bias = pc::infty; // bias is infinite --> reject
                                bias_rejected++; // count how many time we reject due to bias
                            }
#ifndef NDEBUG
                            // check if cluster mass center movement matches displacement
                            if (_bias==0) {
                                Point newCOM = clusterCOM(); // org. cluster center
                                Point d = spc.geo.vdist(COM, newCOM); // distance between new and old COM
                                double _zero = (d+dp).norm(); // |d+dp| should ideally be zero...
                                assert(std::fabs(_zero) < 1e-6 && "cluster likely too large");
                            }
#endif
                        }
                    }

                    double bias(Change&, double, double) override {
                        return _bias;
                    } //!< adds extra energy change not captured by the Hamiltonian

                    void _reject(Change&) override { msqd += 0; msqd_angle += 0; }

                    void _accept(Change&) override {
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
                    std::vector<std::shared_ptr<Potential::BondData>> bonds;
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
                                molecules<Tpvec>[molid].bonds, Potential::BondData::HARMONIC);
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
                            auto g = spc.randomMolecule(molid, slump); // look for random group
                            if (g!=spc.groups.end())
                                if (g->size()>2) { // must at least have three atoms
                                    auto b = slump.sample(bonds.begin(), bonds.end()); // random harmonic bond
                                    if (b != bonds.end()) {
                                        // index in `bonds are relative to the group
                                        int i1 = (*b)->index.at(0);
                                        int i2 = (*b)->index.at(1);
                                        int offset = std::distance( spc.p.begin(), g->begin() );

                                        index.clear();
                                        if (slump()>0.5)
                                            for (size_t i=i2+1; i<g->size(); i++)
                                                index.push_back(i+offset);
                                        else
                                            for (int i=0; i<i1; i++)
                                                index.push_back(i+offset);
                                        i1+=offset;
                                        i2+=offset;

                                        if (not index.empty()) {
                                            Point oldcm = g->cm;
                                            g->unwrap(spc.geo.getDistanceFunc()); // remove pbc
                                            Point u = (spc.p[i1].pos - spc.p[i2].pos).normalized();
                                            double angle = dprot * (slump()-0.5);
                                            Eigen::Quaterniond Q( Eigen::AngleAxisd(angle, u) );
                                            auto M = Q.toRotationMatrix();
                                            for (auto i : index) {
                                                spc.p[i].rotate(Q, M); // internal rot.
                                                spc.p[i].pos = Q * ( spc.p[i].pos - spc.p[i1].pos)
                                                    + spc.p[i1].pos; // positional rot.
                                            }
                                            g->cm = Geometry::massCenter(g->begin(), g->end());
                                            g->wrap(spc.geo.getBoundaryFunc()); // re-apply pbc

                                            d2 = spc.geo.sqdist(g->cm, oldcm); // CM movement

                                            Change::data d;
                                            for (int i : index)
                                                d.atoms.push_back(i-offset); // `atoms` index are relative to group
                                            d.index = Faunus::distance( spc.groups.begin(), g ); // integer *index* of moved group
                                            d.all = false;
                                            d.internal = true;    // trigger internal interactions
                                            change.groups.push_back( d ); // add to list of moved groups
                                        }
                                    }
                                }
                        }
                    }

                    void _accept(Change&) override { msqd += d2; }
                    void _reject(Change&) override { msqd += 0; }

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
                                            spc.geo.getBoundaryFunc(), -g.begin()->pos);
                        }
                    }

                    double exchangeEnergy(double mydu) {
                        std::vector<MPI::FloatTransmitter::floatp> duSelf(1), duPartner;
                        duSelf[0]=mydu;
                        duPartner = ft.swapf(mpi, duSelf, partner);
                        return duPartner.at(0);               // return partner energy change
                    } //!< Exchange energy with partner

                    double bias(Change&, double uold, double unew) override {
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

                    void _accept(Change&) override {
                        if ( goodPartner() )
                            accmap[ id() ] += 1;
                    }

                    void _reject(Change&) override {
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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
                    inline Propagator(const json &j, Tspace &spc, MPI::MPIController &mpi) {
#pragma GCC diagnostic pop
                        if (j.count("random")==1)
                            Movebase::slump = j["random"]; // slump is static --> shared for all moves

                        for (auto &m : j.at("moves")) {// loop over move list
                            size_t oldsize = vec.size();
                            for (auto it : m.items()) {
                                try {
                                    if (it.key()=="moltransrot") this->template push_back<Move::TranslateRotate<Tspace>>(spc); 
                                    else if (it.key()=="biasedmoltransrot") this->template push_back<Move::BiasedTranslateRotate<Tspace>>(spc);
                                    else if (it.key()=="conformationswap") this->template push_back<Move::ConformationSwap<Tspace>>(spc);
                                    else if (it.key()=="transrot") this->template push_back<Move::AtomicTranslateRotate<Tspace>>(spc);
                                    else if (it.key()=="pivot") this->template push_back<Move::Pivot<Tspace>>(spc);
                                    else if (it.key()=="volume") this->template push_back<Move::VolumeMove<Tspace>>(spc);
                                    else if (it.key()=="charge") this->template push_back<Move::ChargeMove<Tspace>>(spc);
                                    else if (it.key()=="speciation") this->template push_back<Move::SpeciationMove<Tspace>>(spc);
                                    else if (it.key()=="quadrantjump") this->template push_back<Move::QuadrantJump<Tspace>>(spc);
                                    else if (it.key()=="cluster") this->template push_back<Move::Cluster<Tspace>>(spc);
                                    // new moves go here...
#ifdef ENABLE_MPI
                                    else if (it.key()=="temper") this->template push_back<Move::ParallelTempering<Tspace>>(spc, mpi);
                                    // new moves requiring MPI go here...
#endif
                                    if (vec.size()==oldsize+1) {
                                        vec.back()->from_json( it.value() );
                                        addWeight(vec.back()->repeat);
                                    } else
                                        std::cerr << "warning: ignoring unknown move '" << it.key() << "'" << endl;
                                } catch (std::exception &e) {
                                    throw std::runtime_error("Error adding move '" + it.key() + "': " + e.what() + usageTip[it.key()]);
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
                        throw std::runtime_error("Metropolis error: energy cannot be NaN");
                    if (du<0)
                        return true;
                    if (-du > pc::max_exp_argument)
                        std::cerr << "warning: large metropolis energy" << std::endl;
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

                State state1, // old state (accepted)
                      state2; // new state (trial)
                double uinit=0, dusum=0;
                Average<double> uavg;

                void init() {
                    dusum=0;
                    Change c; c.all=true;

                    state1.pot.key = Energy::Energybase::OLD; // this is the old energy (current, accepted)
                    state2.pot.key = Energy::Energybase::NEW; // this is the new energy (trial)

                    state1.pot.init();
                    double u1 = state1.pot.energy(c);
                    uinit = u1;

                    state2.sync(state1, c); // copy all information from state1 into state2
                    state2.pot.init();
                    double u2 = state2.pot.energy(c);

                    // check that the energies in state1 and state2 are *identical*
                    if (u1 not_eq u2) {
                        std::cerr << "u1 = " << u1 << "  u2 = " << u2 << endl;
                        throw std::runtime_error("error aligning energies - this could be a bug...");
                    }

                    // inject reference to state1 in SpeciationMove (needed to calc. *differences*
                    // in ideal excess chem. potentials)
                    for (auto base : moves.vec) {
                        auto derived = std::dynamic_pointer_cast<Move::SpeciationMove<Tspace>>(base);
                        if (derived)
                            derived->setOther(state1.spc);
                    }
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
                    if (uinit!=0)
                        return ( ufinal-(uinit+dusum) ) / uinit;
                    else
                        return ufinal-(uinit+dusum); // in case uinit==0, report absolute drift
                } //!< Calculates the relative or absolute energy drift from initial configuration

                MCSimulation(const json &j, MPI::MPIController &mpi) : state1(j), state2(j), moves(j, state2.spc, mpi) {
                    init();
                }

                void store(json &j) const {
                    j = state1.spc;
                    j["random-move"] = Move::Movebase::slump;
                    j["random-global"] = Faunus::random;
                } // store system to json object

                void restore(const json &j) {
                    try {
                        state1.spc = j; // old/accepted state
                        state2.spc = j; // trial state
                        Move::Movebase::slump = j["random-move"]; // restore move random number generator
                        Faunus::random = j["random-global"];      // restore global random number generator
                        reactions<Tpvec> = j.at("reactionlist").get<decltype(reactions<Tpvec>)>(); // should be handled by space
                        init();
                    } catch (std::exception &e) {
                        throw std::runtime_error("error initialising simulation: "s + e.what());
                    }
                } //!< restore system from previously store json object

                void move() {
                    Change change;
                    for (int i=0; i<moves.repeat(); i++) {
                        auto mv = moves.sample(); // pick random move
                        if (mv != moves.end() ) {
                            change.clear();
                            (**mv).move(change);

                            if (change) {
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

                                double bias = (**mv).bias(change, uold, unew) + IdealTerm( state2.spc, state1.spc , change);

                                if ( metropolis(du + bias) ) { // accept move
                                    state1.sync( state2, change );
                                    (**mv).accept(change);
                                } else { // reject move
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
        void to_json(json &j, MCSimulation<Tgeometry,Tparticle> &mc) {
            mc.to_json(j);
        }

    /**
     * @brief Ideal energy contribution of a speciation move
     * This funciton calculates the contribution to the energy change arising from the
     * change in concentration of reactant and products in the current and in the trial state.
     *
     * @f[
     *     \beta \Delta U = - \sum \ln ( N_o!/N_n! V^{N_n - N_o} )
     * @f]
     *
     * where the sum runs over all products and reactants.
     *
     * @todo
     * - use exception message to suggest how to fix the problem
     */
    template<typename Tspace>
        double IdealTerm( Tspace &spc_n, Tspace &spc_o, const Change &change) {
            using Tpvec = typename Tspace::Tpvec;
            double NoverO=0;
            if ( change.dN ) { // Has the number of any molecules changed?
                for ( auto &m : change.groups ) {
                    int N_o = 0;
                    int N_n = 0;
                    if ( m.dNswap ) {
                        assert( m.atoms.size()==1 );
                        auto &g_n = spc_n.groups.at(m.index);
                        int id1 = (g_n.begin()+m.atoms.front())->id;
                        auto &g_o = spc_o.groups.at(m.index);
                        int id2 = (g_o.begin()+m.atoms.front())->id;
                        for (int id : {id1, id2}) {
                            auto atomlist_n = spc_n.findAtoms(id);
                            auto atomlist_o = spc_o.findAtoms(id);
                            N_n = size(atomlist_n);
                            N_o = size(atomlist_o);
                            int dN = N_n - N_o;
                            double V_n = spc_n.geo.getVolume();
                            if (dN>0)
                                for (int n=0; n < dN; n++)
                                    NoverO += std::log( (N_o + 1 + n) / ( V_n * 1.0_molar ));
                            else 
                                for (int n=0; n < (-dN); n++)
                                    NoverO -= std::log( (N_o - n) / ( V_n * 1.0_molar ));
                        }
                    } else {
                        if ( m.dNatomic ) {
                            auto mollist_n = spc_n.findMolecules(spc_n.groups[m.index].id, Tspace::ALL);
                            auto mollist_o = spc_o.findMolecules(spc_o.groups[m.index].id, Tspace::ALL);
                            if ( size(mollist_n) > 1 || size(mollist_o) > 1 )
                                throw std::runtime_error("Bad definition: One group per atomic molecule!");
                            if ( not molecules<Tpvec>[ spc_n.groups[m.index].id ].atomic)
                                throw std::runtime_error("Only atomic molecules!");
                            // Below is safe due to the catches above
                            // add consistency criteria with m.atoms.size() == N
                            N_n =  mollist_n.begin()->size();
                            N_o =  mollist_o.begin()->size();
                        } else {
                            if ( not molecules<Tpvec>[ spc_n.groups[m.index].id ].atomic ) { // Molecular species
                                auto mollist_n = spc_n.findMolecules(m.index, Tspace::ACTIVE);
                                auto mollist_o = spc_o.findMolecules(m.index, Tspace::ACTIVE);
                                N_n=size(mollist_n);
                                N_o=size(mollist_o);
                            }
                        }
                        int dN = N_n - N_o;
                        if (dN!=0) {
                            double V_n = spc_n.geo.getVolume();
                            //double V_o = spc_o.geo.getVolume();
                            if (dN>0)
                                for (int n=0; n < dN; n++)
                                    NoverO += std::log( (N_o + 1 + n) / ( V_n * 1.0_molar ));
                            else 
                                for (int n=0; n < (-dN); n++)
                                    NoverO -= std::log( (N_o - n) / ( V_n * 1.0_molar ));
                        }
                    }
                }
            }
            return NoverO; // negative sign since Pref exp{-beta(dU)} = exp{-beta(dU -ln(Pref)}
        }

}//Faunus namespace
