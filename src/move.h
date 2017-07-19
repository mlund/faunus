#pragma once

#include "core.h"
#include "potentials.h"

namespace Faunus {
    namespace Move {

        template<typename Tspace>
        struct TranslateMove {
            typename Tspace::Tchange change;
            Random slump;
            std::vector<int> mollist; // list of molecule id to pick from
            Tspace spc, trial;

            void move(Tspace &spc) {
                assert(!mollist.empty());

                change.clear();

                int molid = *slump.sample(mollist.begin(), mollist.end()); // random molecule id
                int igroup = trial.random(trial.groups, molid, slump); // random group index
                if (igroup>-1)
                    if (!trial.g[igroup].empty()) {
                        int iatom = trial.g[igroup].random(slump);

                        trial.p[iatom].pos += {1,1,1};
                        change.moved[igroup].push_back(iatom);
                    }

                auto pairpot = Potential::Coulomb() + Potential::HardSphere();
                Nonbonded<Tspace,decltype(pairpot)> pot;
                pot.pairpot = pairpot;
                double du = pot.energy(spc, change);
            }

            void accept() {
                spc.sync(trial, change);
            }

            void reject() {
                trial.sync(spc, change);
            }
        };

        template<class Tspace>
        class MCSimulation {
        private:
            Tspace spc, backup;
            TranslateMove<Tspace> mv;
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