#pragma once

//#include "core.h"
//#include "potentials.h"

namespace Faunus {
    namespace Move {

        struct MoleculeMoveParameters {
            int molid; //!< Molecule type to move
            float prob=1; //!< Probability
            unsigned int trials=0;   //!< Number of trial moves
            unsigned int accepted=0; //!< Number of accepted moves
            Point dir={1,1,1};
        };

        class Movebase {
            private:
                virtual void move(Change&)=0; //!< Perform move and modify change object
            protected:
                Random slump; //!< Temporarily here
            public:
                std::string name;            //!< Name of move
                struct data {
                    int molid;
                    float prob=1;
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

                    void from_json(const json &j) {
                        assert( j.is_object() );
                        for ( auto it = j.begin(); it != j.end(); ++it ) {
                            typedef typename Tspace::Tpvec Tpvec;
                            auto mol = findName(molecules<Tpvec>, it.key());
                            if (mol!=molecules<Tpvec>.end()) {
                                data d;
                                d.molid = mol->id();
                                d.dp1 = j.at("dp").get<double>();
                                d.dp2 = j.at("dprot").get<double>();
                                d.dir = j.value("dir",  d.dir);
                                mollist.push_back(d);
                            }
                            else
                                throw std::runtime_error("Unknown molecule '" + it.key() + "'");
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
            Tspace spc, trial;

            TranslateRotate<Tspace> mv(trial);
            json j;
            mv.from_json(j);
            j = mv;
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
