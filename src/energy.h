#pragma once

#include "core.h"
#include "geometry.h"
#include "space.h"
#include "potentials.h"

namespace Faunus {
    namespace Energy {

        class Energybase {
        public:
            //!< Update energy to reflect `newspc`
            virtual void update(Change&) {
            } //!< Update to reflect changes made in other space

            virtual double energy(Change&) {
                return 0;
            } //!< Return energy due to change
        };

        /**
         * @brief Nonbonded energy using a pair-potential
         */
        template<typename Tspace, typename Tpairpot>
        struct Nonbonded : public Energybase {
            Tspace& oldspc;   //!< Ref. to original space
            Tspace& newspc;   //!< Ref. to new or trial space
            Tpairpot pairpot;
            typedef std::vector<int> Tindex;

            Nonbonded(Tspace &oldspc, Tspace &newspc) : oldspc(oldspc), newspc(newspc) {
            }

            template<typename T>
            double i2i(const Tspace &spc, const T &a, const T &b) {
                return pairpot(a, b, spc.geo.vdist(a.pos, b.pos));
            }

            template<typename T>
            double g2g(const Tspace &spc, const T &g1, const T &g2) {
                double u = 0;
                for (auto &i : g1)
                    for (auto &j : g2)
                        u += i2i(spc, i, j);
                return u;
            }

            double index2index(const Tspace &spc, const Tindex &index1, const Tindex &index2) {
                double u = 0;
                for (auto &i : index1)
                    for (auto &j : index2)
                        u += i2i(spc, spc.p[i], spc.p[j]);
                return u;
            }

            double energy(Change &change) override {
                if (change.empty())
                    return 0;

                double uold=0, unew=0;

                using namespace ranges;
                auto moved = change.touchedGroupIndex(); // index of moved groups
                auto fixed = view::ints( 0, int(oldspc.groups.size()) )
                    | view::remove_if(
                            [&moved](int i){return std::binary_search(moved.begin(), moved.end(), i);}
                            ); // index of static groups

                for (auto &d : change.groups) {

                    // moved<->moved
                    /*
                       for (auto &n : change.groups)
                       if (n.first > m.first) {
                       auto &g1 = spc.groups[m.first];
                       auto &g2 = spc.groups[n.first];
                       u += g2g(spc, g1, g2);
                       }*/

                    // moved<->static
                    for (auto i : fixed) {
                        uold += g2g(oldspc, oldspc.groups[d.index], oldspc.groups[i]);
                        unew += g2g(newspc, newspc.groups[d.index], newspc.groups[i]);
                    }
                }
                return unew-uold;
            }
        }; // Nonbonded energy before and after change (returns pair w. uold and unew)

        template<typename Tspace>
        class Hamiltonian : public Energybase {
        private:
            std::vector<std::shared_ptr<Energy::Energybase>> vec;
        public:
            void update(Change &change) override {
                for (auto i : vec)
                    i->update(change);
            }
            double energy(Change &change) override {
                double du=0;
                for (auto i : vec)
                    du += i->energy(change);
                return du;
            } //!< Energy due to changes
        };

    }//namespace
}//namespace
