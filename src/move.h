#pragma once

#include "core.h"
#include "energy.h"
#include "average.h"
#include "analysis.h"

namespace Faunus {
    namespace Move {

        class Movebase {
            private:
                virtual void _move(Change&)=0; //!< Perform move and modify change object
                virtual void _accept(Change&) {}; //!< Call after move is accepted
                virtual void _reject(Change&) {}; //!< Call after move is rejected
            protected:
                virtual void _to_json(json &j) const=0; //!< Extra info for report if needed
                virtual void _from_json(const json &j)=0; //!< Extra info for report if needed
                double prob=1;
                unsigned long cnt=0;
                unsigned long accepted=0;
                unsigned long rejected=0;
            public:
                static Random slump;   //!< Shared for all moves
                std::string name;      //!< Name of move
                unsigned int repeat=1; //!< How many times the move should be repeated
                bool permol=true;

                inline void from_json(const json &j) {
                    _from_json(j);
                }

                inline void to_json(json &j) const {
                    assert( !name.empty() );
                    auto& _j = j[name];
                    _to_json(_j);
                    _j["acceptance"] = double(accepted)/cnt;
                    _j["N attempted"] = cnt;
                    _j["N accepted"] = accepted;
                    _j["N rejected"] = cnt-accepted;
                } //!< JSON report w. statistics, output etc.

                inline void move(Change &change) {
                    cnt++;
                    change.clear();
                    _move(change);
                } //!< Perform move and modify given change object

                inline void accept(Change &c) {
                    accepted++;
                    _accept(c);
                }

                inline void reject(Change &c) {
                    rejected++;
                    _reject(c);
                }
        };

        Random Movebase::slump; // static instance of Random (shared for all moves)

        inline void from_json(const json &j, Movebase &m) {
            m.from_json(j);
        } //!< Configure any move via json

        inline void to_json(json &j, const Movebase &m) {
            m.to_json(j);
        }

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
                            {"molid", molid}, {"prob", prob},
                            {u8::rootof + u8::bracket("r" + u8::squared), std::sqrt(msqd.avg())},
                            {"mol", molecules<Tpvec>[molid].name}
                        };
                    }

                    void _from_json(const json &j) override {
                        assert(!molecules<Tpvec>.empty());
                        try {
                            std::string molname = j.at("mol");
                            auto it = findName(molecules<Tpvec>, molname);
                            if (it == molecules<Tpvec>.end())
                                throw std::runtime_error("unknown molecule '" + molname + "'");
                            molid = it->id();
                            prob = j.value("prob", 1.0);
                            dir = j.value("dir", Point(1,1,1));
                            dprot = j.at("dprot");
                            dptrans = j.at("dp");
                            permol = j.value("permol", true);
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
                        auto mollist = spc.findMolecules( molid ); // list of molecules w. 'molid'
                        if (size(mollist)>0) {
                            auto it  = (mollist | ranges::view::sample(1, slump.engine) ).begin();
                            assert(it->id==molid);

                            // translate group
                            Point oldcm = it->cm;
                            Point dp = 0.5*ranunit(slump).cwiseProduct(dir) * dptrans;
                            it->translate( dp, spc.geo.boundaryFunc );
                            _sqd = spc.geo.sqdist(oldcm, it->cm); // squared displacement

                            // add rotation here...

                            Change::data d; // object that describes what was moved in a single group
                            d.index = Faunus::distance( spc.groups.begin(), it ); // integer *index* of moved group
                            d.all = true; // *all* atoms in group were moved
                            change.groups.push_back( d ); // add to list of moved groups
                        }
                        else std::cerr << name << ": no molecules found" << std::endl;
                    }

                    void _accept(Change &change) override { msqd += _sqd; }
                    void _reject(Change &change) override { msqd += 0; }

                public:
                    TranslateRotate(Tspace &spc) : spc(spc) {
                        name = "Molecular Translation and Rotation";
                    }
                    void dummy() {
                        std::cout << "dummy\n";
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
            json j = R"( {"mol":"B", "dp":1.0, "dprot":0.5, "dir":[0,1,0], "prob":0.2 })"_json;
            mv.from_json(j);

            j = json(mv).at(mv.name);
            CHECK( j.at("mol")   == "B");
            CHECK( j.at("dir")   == Point(0,1,0) );
            CHECK( j.at("dp")    == 1.0 );
            CHECK( j.at("prob")  == 0.2 );
            CHECK( j.at("dprot") == 0.5 );
            //CHECK( j.at("molid") == 1 );
            //CHECK( j["B"].at("accepted") == 0 );
            //CHECK( j["B"].at("rejected") == 0 );
        }
#endif

        template<typename Tspace>
            class Propagator : public BasePointerVector<Movebase> {
                public:
                    using BasePointerVector<Movebase>::vec;
                    inline Propagator() {}
                    inline Propagator(const json &j, Tspace &spc) {
                        for (auto &m : j.at("moves")) {// loop over move list
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                if (it.key()=="moltransrot") {
                                    this->template push_back<Move::TranslateRotate<Tspace>>(spc);
                                    vec.back()->from_json( it.value() );
                                }
                                // additional moves go here...
                            }
                        }
                    }
            };

        template<typename Tspace>
            void to_json(json &j, const Propagator<Tspace> &m) {
                for (auto i : m.vec)
                    j.push_back(*i);
            }

    }//Move namespace

    template<class Tgeometry, class Tparticle>
        class MCSimulation {
            private:
                bool metropolis(double du) {
                    if (std::isnan(du))
                        return false;
                    if (du<0)
                        return true;
                    return ( Move::Movebase::slump() > std::exp(-du)) ? false : true;
                } //!< Metropolis criterion (true=accept)

            public:
                typedef Space<Tgeometry, Tparticle> Tspace;
                typedef typename Tspace::Tpvec Tpvec;

                struct State {
                    Tspace spc;
                    Energy::Hamiltonian<Tspace> pot;
                    State(const json &j) : spc(j), pot(spc,j) {
                    }
                }; //!< Contains everything to describe a state

                State state1, state2;
                Move::Propagator<Tspace> moves;
                double uinit=0, dusum=0;
                Analysis::SystemEnergy Usys;

                const auto& geo() { return state1.spc.geo; }
                const Tpvec& p() { return state1.spc.p; }

                MCSimulation(const json &j) : state1(j), state2(j), moves(j, state2.spc), Usys(j.at("analysis"), state1.pot) {
                    assert( state1.spc.groups.front().begin() == state1.spc.p.begin());
                    assert( state2.spc.groups.front().begin() == state2.spc.p.begin());

                    Change c;
                    c.all=true;
                    state2.spc.sync(state1.spc, c);
                    assert(state1.spc.p.begin() != state2.spc.p.begin());
                    assert(state1.spc.groups.front().begin() == state1.spc.p.begin());
                    assert(state2.spc.groups.front().begin() == state2.spc.p.begin());

                    c.all=true;
                    double uinit1 = state1.pot.energy(c);
                    double uinit2 = state2.pot.energy(c);
                    uinit = uinit1;

                    cout << "uinit         = " << uinit1 << endl;
                    cout << "uinit (state2) = " << uinit2 << endl;
                }

                ~MCSimulation() {
                    Change c;
                    c.all=true;
                    double ufinal = state1.pot.energy(c);

                    cout << "initial energy = " << uinit << endl;
                    cout << "final energy   = " << ufinal << endl;
                    cout << "sum of changes = " << dusum << endl;
                    cout << "drift          = " << ufinal-(uinit+dusum) << endl;
                    std::ofstream f("out.json");
                    json j;
                    j["moves"] = json(moves);
                    j["space"] = state1.spc;
                    j["energy"] = json(state1.pot);
                    if (f) {
                        f << std::setw(4) << j << endl;
                        f.close();
                    }
                }

                void move() {

                    Change change;

                    // pick move
                    auto mv = Move::Movebase::slump.sample( moves.begin(), moves.end() );
                    if (mv == moves.end() ) {
                        std::cerr << "no mc moves found" << endl;
                    } else {
                        (**mv).move(change);

                        if (!change.empty()) {

                            double unew = state2.pot.energy(change),
                                   uold = state1.pot.energy(change),
                                   du = unew - uold;

                            if (metropolis(du) ) // accept
                            {
                                state1.spc.sync( state2.spc, change ); // sync newspc -> state1.pc
                                (**mv).accept(change);
                            }
                            else // reject
                            {
                                state2.spc.sync( state1.spc, change ); // copy state1.pc -> newspc
                                (**mv).reject(change);
                                du=0;
                            }
                            dusum+=du;
                        }
                    }
                    Usys.sample();
                }

        };

}//namespace
