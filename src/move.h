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
                Random slump;     //!< Temporarily here
                inline virtual json _to_json() const { return json(); } //!< Extra info for report if needed
                double prob=1;
                unsigned long cnt=0;
                unsigned long accepted=0;
                unsigned long rejected=0;
            public:
                std::string name;      //!< Name of move
                unsigned int repeat=1; //!< How many times the move should be repeated
                bool permol=true;

                inline json to_json() const {
                    assert( !name.empty() );
                    json j = {{ name, _to_json() }};
                    j[name]["acceptance"] = double(accepted)/cnt;
                    j[name]["cnt"] = accepted;
                    j[name]["accepted"] = accepted;
                    j[name]["rejected"] = cnt-accepted;
                    return j;
                } //!< JSON report w. statistics, output etc.

                inline void move(Change &change) {
                    change.clear();
                    _move(change);
                    cnt++;
                } //!< Perform move and modify given change object

                void accept(Change &c) {
                    accepted++;
                    _accept(c);
                }

                void reject(Change &c) {
                    rejected++;
                    _reject(c);
                }
        };

        inline void to_json(json &j, const Movebase &m) {
            j = m.to_json();
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

                    json _to_json() const override {
                        return {
                            {"dir", dir}, {"dp", dptrans}, {"dprot", dprot},
                            {"molid", molid}, {"prob", prob},
                            {"group", molecules<Tpvec>[molid].name}
                        };
                    }

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

                    void from_json(const json &j) {
                        assert(!molecules<Tpvec>.empty());
                        try {
                            std::string molname = j.at("group");
                            auto it = findName(molecules<Tpvec>, molname);
                            if (it == molecules<Tpvec>.end())
                                throw std::runtime_error("unknown molecule '" + molname + "'");
                            name += ": [" + molname + "]";
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
                    } //!< Configure via json object and check input syntax
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

            j = json(mv)[mv.name];
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
    }//namespace

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
                std::vector<std::shared_ptr<Move::Movebase>> moves;
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
                                auto ptr = std::make_shared<Move::TranslateRotate<Tspace>>(old.spc, trial.spc);
                                ptr->from_json( it.value() );
                                moves.push_back(ptr);
                            }
                } // add MC moves

            public:
                typedef typename Tspace::Tparticle Tparticle;
                typedef typename Tspace::Tpvec Tpvec;

                const Tpvec& p() { return old.spc.p; }
                const auto& geo() { return old.spc.geo; }

                MCSimulation(const json &j) {
                    atoms<Tparticle> = j.at("atomlist").get<decltype(atoms<Tparticle>)>();
                    molecules<Tpvec> = j.at("moleculelist").get<decltype(molecules<Tpvec>)>();
                    old.spc.geo = j.at("system").at("geometry");
                    insertMolecules(old.spc);
                    cout << old.spc.p.size() << endl;
                    Change c;
                    c.dV=1;
                    trial.spc.sync(old.spc, c);
                    assert(old.spc.p.size() == trial.spc.p.size());
                    addMoves(j);
                }

                ~MCSimulation() {
                    std::ofstream f("out.json");
                    json j0, j1, j2;
                    j0 = moves.front()->to_json();
                    j1["atomlist"] = atoms<Tparticle>;
                    j2["moleculelist"] = molecules<Tpvec>;
                    j0 = merge(j0, j1);
                    j0 = merge(j0, j2);
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

                    auto mv = *slump.sample( moves.begin(), moves.end() ); // pick random move
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
