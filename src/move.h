#pragma once

#include "energy.h"

namespace Faunus {
    namespace Move {

        class Movebase {
            private:
                virtual void _move(Change&)=0; //!< Perform move and modify change object
                virtual void _accept(Change&)=0; //!< Call after move is accepted
                virtual void _reject(Change&)=0; //!< Call after move is rejected
            protected:
                Random slump;     //!< Temporarily here
                json config;      //!< JSON object containing
                inline virtual json _to_json() const { return json(); } //!< Extra info for report if needed
                double prob=1;
                unsigned long cnt=0;
                unsigned long accepted=0;
            public:
                std::string name; //!< Name of move

                inline json to_json() const {
                    assert( !name.empty() );

                    //config["accepted"] = accepted;
                    //config["rejected"] = cnt-accepted;

                    json j1 = {{ name, _to_json() }};
                    json j2 = {{ name, config }};
                    return merge( j2, j1 );
                } //!< JSON report w. statistics, output etc.

                inline void operator()(Change &change) {
                    change.clear();
                    _move(change);
                } //!< Perform move and modify given change object
        };

        inline void to_json(json &j, const Movebase &m) {
            j = m.to_json();
        }

        /**
         * @brief Translate and rotate a molecular group
         *
         * Input format:
         *
         * ```.js
         *     { "mollist" :
         *         {
         *             "water" : { "dp" : 1.0, "dprot" : 0.5, "dir" : [1,1,1], "prob" : 1.0 }
         *         }
         *     }
         * ```
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
                    Point dir;

                    json _to_json() const override {
                        //config["molid"] = molid;
                        return json();
                    }

                    void _move(Change &change) override {
                        assert(molid>=0);

                        // pick random group from the system matching molecule type
                        auto g_iter = ( newspc->findMolecules( molid )
                                | ranges::view::sample(1, slump.engine) ).begin(); // iterator to random group

                        // translate group
                        Point dp = ranunit(slump).cwiseProduct(dir) * dptrans;
                        g_iter->translate( dp, newspc->geo.boundaryFunc );
                        // add rotation here...

                        Change::data d; // object that describes what was moved in a single group
                        d.index = Faunus::distance( newspc->groups.begin(), g_iter ); // integer *index* of moved group
                        d.all = true; // *all* atoms in group were moved
                        change.groups.push_back( d ); // add to list of moved groups
                    }

                    void _accept(Change &change) override {
                    }

                    void _reject(Change &change) override {
                    }

                public:
                    TranslateRotate(Tspace &oldspc, Tspace &newspc) : oldspc(&oldspc), newspc(&newspc) {
                        name = "Molecular Translation and Rotation";
                        //assert( &oldspc!=nullptr );
                        //assert( &newspc!=nullptr );
                    }

                    void from_json(const json &j) {
                        try {
                            config = j;
                            std::string molname = j.at("group");
                            name += ": " + molname;
                            auto it = findName(molecules<Tpvec>, molname);
                            if (it == molecules<Tpvec>.end())
                                throw std::runtime_error("unknown molecule '" + molname + "'");
                            molid = it->id();
                            prob = j.value("prob", 1.0);
                            dir = j.value("dir", Point(1,1,1));
                            dprot = j.at("dprot");
                            dptrans = j.at("dp");
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
            private:
                Random slump;
                Tspace oldspc, newspc;
                std::vector<std::shared_ptr<Move::Movebase>> moves;
                Energy::Nonbonded<Tspace, Potential::Dummy> pot;

                bool metropolis(double du) {
                    if (du<0)
                        return true;
                    return ( slump() > std::exp(-du)) ? false : true;
                } //!< Metropolis criterion (true=accept)

            public:
                typedef typename Tspace::Tparticle Tparticle;
                typedef typename Tspace::Tpvec Tpvec;

                MCSimulation(const json &j) : oldspc(j), newspc(oldspc) {
                    for (auto &m : j.at("moves")) // loop over all moves
                        for (auto it=m.begin(); it!=m.end(); ++it)
                            if (it.key()=="moltransrot") {
                                auto ptr = std::make_shared<Move::TranslateRotate<Tspace>>(oldspc, newspc);
                                ptr->from_json( it.value() );
                                moves.push_back(ptr);
                            }
                }

                ~MCSimulation() {
                    //std::cout << std::setw(4) << mv.to_json() << endl;
                }

                void move() {
                    typename Tspace::Tchange change;
                    //mv(change);
                    return;
                    if (!change.empty()) {
                        cout << "change!" << endl;
                        auto u = pot.energy(oldspc, newspc, change);
                        double du = u.second - u.first;
                        if (metropolis(du) ) // accept
                        {
                            oldspc.sync( newspc, change ); // copy newspc -> oldspc
                        }
                        else // reject
                        {
                            newspc.sync( oldspc, change ); // copy oldspc -> newspc
                        }
                    }
                }

        };

}//namespace
