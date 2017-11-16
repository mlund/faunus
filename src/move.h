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
            public:
                std::string name; //!< Name of move

                inline json to_json() const {
                    assert( !name.empty() );
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
                    json* mollist = nullptr;   // json list w. molecule types to move
                    json::iterator currentmol; // iterator to json entry of selected molecule type

                    void _move(Change &change) override {
                        assert( mollist!=nullptr );

                        // pick random molecule type from `mollist`
                        currentmol = slump.sample(mollist->begin(), mollist->end());
                        json& j = currentmol.value(); // json entry w. displacement parameters etc.

                        // pick random group from the system matching molecule type
                        auto g_iter = ( newspc->findMolecules( j["molid"].get<int>() )
                                | ranges::view::sample(1, slump.engine) ).begin(); // iterator to random group

                        // translate group
                        Point dp = ranunit(slump).cwiseProduct( j["dir"].get<Point>() ) * j["dp"].get<double>();
                        g_iter->translate( dp, newspc->geo.boundaryFunc );
                        // add rotation here...

                        Change::data d; // object that describes what was moved in a single group
                        d.index = Faunus::distance( newspc->groups.begin(), g_iter ); // integer *index* of moved group
                        d.all = true; // *all* atoms in group were moved
                        change.groups.push_back( d ); // add to list of moved groups
                    }

                    void _accept(Change &change) override {
                        currentmol.value()["accepted"]+=1;
                    }

                    void _reject(Change &change) override {
                        currentmol.value()["rejected"]+=1;
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
                            mollist = &config.at("mollist");
                            if (mollist->is_object()) {
                                for (auto it=mollist->begin(); it!=mollist->end(); ++it) {
                                    auto mol = findName(molecules<Tpvec>, it.key());
                                    if (mol == molecules<Tpvec>.end())
                                        throw std::runtime_error("unknown molecule '" + it.key() + "'");
                                    else {
                                        auto &v = it.value();
                                        v["molid"] = mol->id();
                                        v["dir"] = v.value("dir", Point(1,1,1));
                                        v["prob"] = v.value("prob", 1.0);
                                        v["accepted"] = v.value("accepted", int());
                                        v["rejected"] = v.value("rejected", int());
                                        if (!v.at("dp").is_number() ||
                                                !v.at("dprot").is_number() )
                                            throw std::runtime_error("'dp' and 'dprot' must be numbers");
                                    }
                                }
                            } else
                                throw std::runtime_error("'mollist' must be of type object");
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
            json j = R"( {"mollist" : { "B":{"dp":1.0, "dprot":0.5, "dir":[0,1,0], "prob":0.2 } }})"_json;
            mv.from_json(j);

            j = json(mv)[mv.name]["mollist"];
            CHECK( j["B"].at("dir")   == Point(0,1,0) );
            CHECK( j["B"].at("dp")    == 1.0 );
            CHECK( j["B"].at("prob")  == 0.2 );
            CHECK( j["B"].at("dprot") == 0.5 );
            CHECK( j["B"].at("molid") == 1 );
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
                Move::TranslateRotate<Tspace> mv;
                Energy::Nonbonded<Tspace, Potential::Coulomb> pot;

                bool metropolis(double du) {
                    if (du<0)
                        return true;
                    return ( slump() > std::exp(-du)) ? false : true;
                } //!< Metropolis criterion (true=accept)

            public:
                MCSimulation() : mv(oldspc, newspc) {
                }

                void move() {
                    typename Tspace::Tchange change;
                    mv(change);
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
