#pragma once

//#include "core.h"
//#include "potentials.h"

namespace Faunus {
    namespace Move {

        class Movebase {
            private:
                virtual void move(Change&)=0; //!< Perform move and return change object
            protected:
                Random slump; //!< Temporarily here
            public:
                std::string name;            //!< Name of move
                struct data {
                    int molid;
                    double weight=1;
                    unsigned int trials=0;   //!< Number of trial moves
                    unsigned int accepted=0; //!< Number of accepted moves
                    std::vector<double> prop_f;
                    std::vector<Point> prop_v;
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
                    Tspace* trial;

                    enum keys {DP=0, DPROT=1, DIR=0};

                    void move(Change &change) override {
                        auto d = randomMolecule(); // random registered molid
                        auto it = ( trial->findMolecules(d->molid)
                            | ranges::view::sample(1, slump.engine) ).begin(); // random group

                        Point dp = ranunit(slump).cwiseProduct( d->prop_v.at(DIR) ) * d->prop_f.at(DP);
                        it->translate( dp, trial->geo.boundaryFunc );

                        Change::data cdata;
                        cdata.index = Faunus::distance( trial->groups.begin(), it );
                        cdata.all = true;
                        change.groups.push_back( cdata );
                    }

                public:
                    TranslateRotate(Tspace &trial) : trial(&trial) {
                        name = "Translate-Rotate";
                    } 

                    void from_json(const json &j) {
                        typedef typename Tspace::Tpvec Tpvec;
                        for ( auto it = j.begin(); it != j.end(); ++it ) {
                            auto mol = findName(molecules<Tpvec>, it.key());
                            if (mol!=molecules<Tpvec>.end()) {
                            }
                            else
                                throw std::runtime_error("Unknown molecule '" + it.key() + "'");
                        }
                    } //!< Configure via json object

            };

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] TranslateRotate")
        {
            typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
            typedef Space<Geometry::Cuboid, Tparticle> Tspace;
            Tspace spc, trial;

            TranslateRotate<Tspace> mv(trial);
            json j;
            mv.from_json(j);
        }
#endif

        template<class Tspace>
            class MCSimulation {
                private:
                    Tspace spc, backup;
                    TranslateRotate<Tspace> mv;
                    //Propagator mv;
                    //Reporter analyse;

                public:
                    MCSimulation() {
                    }

                    void move() {
                        typename Tspace::Tchange change;
                        // make change to trial
                        //std::pair<double> u = pot.energyChange( spc, trial, change );
                        //if (metropolis( u.second-u.first ) )
                        //    spc.sync( trial, change );
                        //else
                        //    trial.sync( spc, change );
                    }

            };

    }//namespace
}//namespace
