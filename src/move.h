#pragma once

#include "energy.h"

namespace Faunus {
    namespace Move {

        class Movebase {
            private:
                virtual void move(Change&)=0; //!< Perform move and modify change object
            protected:
                Random slump; //!< Temporarily here
            public:
                json config;
                std::string name;            //!< Name of move

                auto randomGroup() {
                    auto &j = config.at("groups");
                    return slump.sample( j.begin(), j.end() );
                } //!< iterator to random group to move

                inline void operator()(Change &change) {
                    change.clear();
                    move(change);
                } //!< Perform move and return change object
        };

        template<typename Tspace>
            class TranslateRotate : public Movebase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace* trial;
                    json::iterator currentmol; // iterator to json entry of random molecule type

                    void move(Change &change) override {
                        currentmol = randomGroup(); // iterator to json entry of random molecule id
                        json& d = currentmol.value();

                        auto g_iter = ( trial->findMolecules( d["molid"].get<int>() )
                                | ranges::view::sample(1, slump.engine) ).begin(); // iterator to random group

                        Point dp = ranunit(slump).cwiseProduct( d["dir"].get<Point>() ) * d["dp"].get<double>();
                        g_iter->translate( dp, trial->geo.boundaryFunc );

                        Change::data cdata; // object that describes what was moved
                        cdata.index = Faunus::distance( trial->groups.begin(), g_iter ); // index of moved group
                        cdata.all = true; // *all* atoms were moved
                        change.groups.push_back( cdata ); // add to list of moved groups
                    }

                public:
                    TranslateRotate(Tspace &trial) : trial(&trial) {
                        name = "Translate-Rotate";
                    }

                    json to_json() const {
                        return config;
                    }

                    void from_json(const json &j) {
                        try {
                            config = j;
                            auto& g = config.at("mollist");
                            if (g.is_object()) {
                                for (auto it=g.begin(); it!=g.end(); ++it) {
                                    auto mol = findName(molecules<Tpvec>, it.key());
                                    if (mol == molecules<Tpvec>.end())
                                        throw std::runtime_error("unknown molecule '" + it.key() + "'");
                                    else {
                                        auto &v = it.value();
                                        v["molid"] = mol->id();
                                        v["dir"] = v.value("dir", Point(1,1,1));
                                        v["prob"] = v.value("prob", 1.0);
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
                    } //!< Configure via json object
            };

        template<class T>
            void to_json(json &j, const TranslateRotate<T> mv) {
                j = mv.config;
            }

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] TranslateRotate")
        {
            typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
            typedef Space<Geometry::Cuboid, Tparticle> Tspace;
            typedef typename Tspace::Tpvec Tpvec;
            Tspace spc, trial;

            CHECK( !atoms<Tparticle>.empty() ); // set in a previous test
            CHECK( !molecules<Tpvec>.empty() ); // set in a previous test

            TranslateRotate<Tspace> mv(trial);
            json j = R"( {"mollist" : { "B":{"dp":1.0, "dprot":0.5, "dir":[0,1,0]}  }})"_json;
            mv.from_json(j);

            CHECK( mv.config["mollist"]["B"].at("dp")    == 1.0 );
            CHECK( mv.config["mollist"]["B"].at("dir")   == Point(0,1,0) );
            CHECK( mv.config["mollist"]["B"].at("dprot") == 0.5 );
            CHECK( mv.config["mollist"]["B"].at("molid") == 1 );
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
                    if ( slump() > std::exp(-du)) // core of MC!
                        return false;
                    return true;
                } //!< Metropolis criterion


            public:
                MCSimulation() : mv(newspc) {
                }

                void move() {
                    typename Tspace::Tchange change;
                    mv(change);
                    if (!change.empty()) {
                        cout << "change!" << endl;
                        auto u = pot.energy(oldspc, newspc, change);
                        double du = u.second - u.first;
                        if (metropolis(du) )
                            oldspc.sync( newspc, change );
                        else
                            newspc.sync( oldspc, change );
                    }
                    // make change to trial
                    //std::pair<double> u = pot.energyChange( spc, trial, change );
                    //if (metropolis( u.second-u.first ) )
                    //    spc.sync( trial, change );
                    //else
                    //    trial.sync( spc, change );
                }

        };

}//namespace
