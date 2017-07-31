#pragma once

//#include "core.h"
//#include "potentials.h"

namespace Faunus {
    namespace Move {

        class Movebase {
            private:
                virtual void move(Change&)=0; //!< Perform move and return change object
            protected:
                static Random slump; //!< Temporarily here
            public:
                std::string name;            //!< Name of move
                struct data {
                    int molid;
                    double weight=1;
                    unsigned int trials=0;   //!< Number of trial moves
                    unsigned int accepted=0; //!< Number of accepted moves
                };
                std::vector<data> mollist;   //!< Vector of molecule id's to operate on

                inline auto randomMolecule() {
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
                    Tspace& trial;

                    void move(Change &change) override {
                        auto it = randomMolecule(); // random registered molid
                        auto g = trial.findMolecules(it->molid)
                            | ranges::view::sample(1, slump.engine); // random group
                        //Point cm = g.cm;// = Point(0,0,0);
                    }

                public:
                    TranslateRotate(Tspace &trial) : trial(trial) {
                        name = "Translate-Rotate";
                    } 
            };
#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] TranslateRotate")
    {
        typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
        typedef Space<Geometry::Cuboid, Tparticle> Tspace;
        Tspace spc, trial;

        //TranslateRotate<Tspace> mv(trial);
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
