#pragma once

#include "core.h"
#include "energy.h"
#include "average.h"

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
                Random slump;     //!< Temporarily here
                double prob=1;
                unsigned long cnt=0;
                unsigned long accepted=0;
                unsigned long rejected=0;
            public:
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
                    _j["cnt"] = accepted;
                    _j["accepted"] = accepted;
                    _j["rejected"] = cnt-accepted;
                } //!< JSON report w. statistics, output etc.

                inline void move(Change &change) {
                    change.clear();
                    _move(change);
                    cnt++;
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
                    Tspace* oldspc = nullptr;  // pointer to old space
                    Tspace* newspc = nullptr;  // pointer to new (trial) space

                    int molid=-1;
                    double dptrans=0;
                    double dprot=0;
                    Point dir={1,1,1};
                    Average<double> msd; // mean square displacement

                    void _to_json(json &j) const override {
                        j = {
                            {"dir", dir}, {"dp", dptrans}, {"dprot", dprot},
                                {"molid", molid}, {"prob", prob},
                                {"group", molecules<Tpvec>[molid].name}
                        };
                    }

                    void _from_json(const json &j) override {
                        assert(oldspc!=nullptr);
                        assert(newspc!=nullptr);
                        assert(!molecules<Tpvec>.empty());
                        try {
                            std::string molname = j.at("group");
                            auto it = findName(molecules<Tpvec>, molname);
                            if (it == molecules<Tpvec>.end())
                                throw std::runtime_error("unknown molecule '" + molname + "'");
                            //name += ": [" + molname + "]";
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

                        // pick random group from the system matching molecule type
                        auto mollist =  newspc->findMolecules( molid ); // list of molecules w. 'molid'
                        auto g_iter  = (mollist | ranges::view::sample(1, slump.engine) ).begin(); // iterator to random group

                        // translate group
                        Point dp = ranunit(slump).cwiseProduct(dir) * dptrans;
                        g_iter->translate( dp, newspc->geo.boundaryFunc );
                        // add rotation here...

                        Change::data d; // object that describes what was moved in a single group
                        d.index = Faunus::distance( newspc->groups.begin(), g_iter ); // integer *index* of moved group
                        d.all = true; // *all* atoms in group were moved
                        change.groups.push_back( d ); // add to list of moved groups
                    }

                public:
                    TranslateRotate(Tspace &oldspc, Tspace &newspc) : oldspc(&oldspc), newspc(&newspc) {
                        name = "Molecular Translation and Rotation";
                    }
            };

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] TranslateRotate")
        {
            typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
            typedef Space<Geometry::Cuboid, Tparticle> Tspace;
            typedef typename Tspace::Tpvec Tpvec;
            Tspace spc, trial;

            CHECK( !atoms<Tparticle>.empty() ); // set in a previous test
            CHECK( !molecules<Tpvec>.empty() ); // set in a previous test

            Tspace oldspc, newspc;
            TranslateRotate<Tspace> mv(oldspc, newspc);
            json j = R"( {"group":"B", "dp":1.0, "dprot":0.5, "dir":[0,1,0], "prob":0.2 })"_json;
            mv.from_json(j);

            j = json(mv).at(mv.name);
            CHECK( j.at("group") == "B");
            CHECK( j.at("dir")   == Point(0,1,0) );
            CHECK( j.at("dp")    == 1.0 );
            CHECK( j.at("prob")  == 0.2 );
            CHECK( j.at("dprot") == 0.5 );
            //CHECK( j.at("molid") == 1 );
            //CHECK( j["B"].at("accepted") == 0 );
            //CHECK( j["B"].at("rejected") == 0 );
        }
#endif

        class Propagator : public BasePointerVector<Movebase> {
            public:
                inline double move(Change &change) {
                    for (auto i : this->vec)
                        i->move(change);
                    return 0;
                }
        };

        inline void to_json(json &j, const Propagator &m) {
            for (auto i : m.vec)
               j = *i;
        }

        inline void from_json(const json &j, Propagator &m) {
        }


    }//Move namespace

    template<class Tspace>
        class MCSimulation {
            public:
                struct state {
                    Tspace spc;
                    Energy::Hamiltonian<Tspace> pot;
                }; //!< Contains everything to describe a state

                state old, trial;

            private:
                Random slump;
                //Energy::Nonbonded<Tspace, Potential::Dummy> pot;

                bool metropolis(double du) {
                    if (du<0)
                        return true;
                    return ( slump() > std::exp(-du)) ? false : true;
                } //!< Metropolis criterion (true=accept)

                void addMoves(const json &j) {
                    for (auto &m : j.at("moves")) // loop over all moves
                        for (auto it=m.begin(); it!=m.end(); ++it)
                            if (it.key()=="moltransrot") {
                                moves.push_back<Move::TranslateRotate<Tspace>>(old.spc, trial.spc);
                                moves.vec.back()->from_json( it.value() );
                            }
                } // add MC moves

            public:
                typedef typename Tspace::Tparticle Tparticle;
                typedef typename Tspace::Tpvec Tpvec;

                const Tpvec& p() { return old.spc.p; }
                const auto& geo() { return old.spc.geo; }
                Move::Propagator moves;

                MCSimulation(const json &j) {
                    atoms<Tparticle> = j.at("atomlist").get<decltype(atoms<Tparticle>)>();
                    molecules<Tpvec> = j.at("moleculelist").get<decltype(molecules<Tpvec>)>();
                    old.spc.geo = j.at("system").at("geometry");
                    insertMolecules(old.spc);
                    Change c;
                    c.dV=1;
                    trial.spc.sync(old.spc, c);
                    assert(old.spc.p.size() == trial.spc.p.size());
                    addMoves(j);

                    old.pot.template push_back<Energy::Nonbonded<Tspace, Potential::HardSphere>>(old.spc);
                    trial.pot.template push_back<Energy::Nonbonded<Tspace, Potential::HardSphere>>(trial.spc);
                }

                ~MCSimulation() {
                    std::ofstream f("out.json");
                    json j0, j1, j2;
                    moves.vec.front()->to_json(j0);
                    if (f) {
                        f << std::setw(4) << j0 << endl;
                        //f << std::setw(4) << moves.front()->to_json() << endl;
                        //f << std::setw(4) << j1 << endl;
                        //f << std::setw(4) << j2 << endl;
                        f.close();
                    }
                }

                void move() {
                    typename Tspace::Tchange change;

                    auto mv = *slump.sample( moves.vec.begin(), moves.vec.end() ); // pick random move
                    mv->move(change);

                    if (!change.empty()) {

                        trial.pot.update(change); // update trial potential
                        double unew = trial.pot.energy(change); // new energy
                        double uold = old.pot.energy(change); // old energy
                        double du = unew - uold;

                        if (metropolis(du) ) // accept
                        {
                            old.spc.sync( trial.spc, change ); // sync newspc -> oldspc
                            // todo: sync hamiltonian
                            mv->accept(change);
                        }
                        else // reject
                        {
                            trial.spc.sync( old.spc, change ); // copy oldspc -> newspc
                            // todo: sync hamiltonian
                            mv->reject(change);
                        }
                    }
                }

        };

}//namespace
