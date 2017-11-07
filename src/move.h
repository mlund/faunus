#pragma once

#include "energy.h"

namespace Faunus {
    namespace Move {

        struct MoleculeMoveParameters {
            int molid; //!< Molecule type to move
            float prob=1; //!< Probability
            unsigned int trials=0;   //!< Number of trial moves
            unsigned int accepted=0; //!< Number of accepted moves
            Point dir={1,1,1}; //!< Direction of move if appropriate
        };

        class Movebase {
            private:
                virtual void move(Change&)=0; //!< Perform move and modify change object
            protected:
                Random slump; //!< Temporarily here
            public:
                std::string name;            //!< Name of move
                struct data {
                    int molid;               //!< Molecule id
                    float prob=1;            //!< Probability
                    unsigned int trials=0;   //!< Number of trial moves
                    unsigned int accepted=0; //!< Number of accepted moves
                    double dp1=0;
                    double dp2=0;
                    Point dir={1,1,1};
                };
                std::vector<data> mollist;   //!< Vector of molecule id's to operate on

                auto randomMolecule() {
                    return slump.sample( mollist.begin(), mollist.end() );
                } //!< Iterator to random molecule (data object)

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

                    void move(Change &change) override {
                        auto d = randomMolecule(); // random registered molid
                        auto it = ( trial->findMolecules(d->molid)
                                | ranges::view::sample(1, slump.engine) ).begin(); // random group

                        Point dp = ranunit(slump).cwiseProduct( d->dir ) * d->dp1;
                        it->translate( dp, trial->geo.boundaryFunc );

                        Change::data cdata; // object that describes what was moved
                        cdata.index = Faunus::distance( trial->groups.begin(), it ); // index of moved group
                        cdata.all = true; // *all* atoms were moved
                        change.groups.push_back( cdata ); // add to list of moved groups
                    }

                public:
                    TranslateRotate(Tspace &trial) : trial(&trial) {
                        name = "Translate-Rotate";
                    }

                    json to_json() const {
                        json j;
                        auto& _j = j[name] = json::array();
                        for (auto &i : mollist) {
                            auto& molname = molecules<Tpvec>.at(i.molid).name;
                            _j.push_back({{ molname, {
                                    { "trials", i.trials },
                                    { "dp", i.dp1 },
                                    { "dprot", i.dp2 },
                                    { "dir", i.dir } }}});
                            if (i.trials>0)
                                _j.push_back( {{ molname, {
                                        { "acceptance", i.accepted / double(i.trials) }
                                        }}});
                        }
                        return j;
                    }

                    void from_json(const json &j) {
                        try {
                            if (j.is_array()) {
                                mollist.clear();
                                for (auto &x : j) {
                                    for (auto it=x.begin(); it!=x.end(); ++it) {
                                        auto mol = findName(molecules<Tpvec>, it.key());
                                        if (mol == molecules<Tpvec>.end())
                                            throw std::runtime_error("unknown molecule '" + it.key() + "'");
                                        else {
                                            data d;
                                            d.molid = mol->id();
                                            d.dp1 = it.value().at("dp").get<double>();
                                            d.dp2 = it.value().at("dprot").get<double>();
                                            d.dir = it.value().value("dir", d.dir);
                                            mollist.push_back(d);
                                        }
                                    }
                                }
                            } else
                                throw std::runtime_error("json object must be of type array.");
                        }
                        catch (std::exception &e) {
                            std::cerr << name << ": " << e.what();
                            throw;
                        }
                    } //!< Configure via json object
            };

        template<class T>
            void to_json(json &j, const TranslateRotate<T> mv) {
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
            json j = R"( [ { "A" : { "dp":1.0, "dprot":0.5, "dir":[0,1,0] } } ] )"_json;
            mv.from_json(j);
            cout << std::setw(4) << mv.to_json() << endl;
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
