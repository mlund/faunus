#pragma once

#include "core.h"
#include "geometry.h"
#include "space.h"
#include "potentials.h"

namespace Faunus {
    namespace Energy {

        class Energybase {
            public:
                std::string name;
                virtual double energy(Change&)=0; //!< Return energy due to change
        };

        /**
         * @brief Nonbonded energy using a pair-potential
         */
        template<typename Tspace, typename Tpairpot>
            struct Nonbonded : public Energybase {
                Tspace& spc;   //!< Space to operate on
                Tpairpot pairpot; //!< Pair potential

                Nonbonded(Tspace &spc, const json &j) : spc(spc) {
                    name="nonbonded";
                    pairpot = j;
                }

                template<typename T>
                    inline double i2i(const T &a, const T &b) {
                        return pairpot(a, b, spc.geo.vdist(a.pos, b.pos));
                    }

                template<typename T>
                    double g2g(const T &g1, const T &g2) {
                        double u = 0;
                        for (auto &i : g1)
                            for (auto &j : g2)
                                u += i2i(i, j);
                        return u;
                    }

                template<typename T>
                    double index2index(const T &index1, const T &index2) {
                        double u = 0;
                        for (auto i : index1)
                            for (auto j : index2)
                                u += i2i(spc.p[i], spc.p[j]);
                        return u;
                    }

                double energy(Change &change) override {
                    double u=0;
                    if (!change.empty()) {
                        using namespace ranges;
                        auto moved = change.touchedGroupIndex(); // index of moved groups
                        auto fixed = view::ints( 0, int(spc.groups.size()) )
                            | view::remove_if(
                                    [&moved](int i){return std::binary_search(moved.begin(), moved.end(), i);}
                                    ); // index of static groups

                        // moved<->moved
                        for ( auto i = moved.begin(); i != moved.end(); ++i )
                            for ( auto j=i; ++j != moved.end(); )
                                u += g2g( spc.groups[*i], spc.groups[*j] );

                        // moved<->static
                        for ( auto i : moved)
                            for ( auto j : fixed) {
                                u += g2g(spc.groups[i], spc.groups[j]);
                            }

                        // more todo!
                    }
                    return u;
                }
            }; // Nonbonded energy before and after change (returns pair w. uold and unew)

        template<class T1, class T2>
            void to_json(json &j, const Nonbonded<T1,T2> &nb) {
                j[nb.name] = json(nb.pairpot);
            }

        template<typename Tspace>
            class Hamiltonian : public Energybase, public BasePointerVector<Energybase> {
                public:
                    double energy(Change &change) override {
                        double du=0;
                        for (auto i : this->vec)
                            du += i->energy(change);
                        return du;
                    } //!< Energy due to changes
            };

        template<typename Tspace>
            void to_json(json &j, const Hamiltonian<Tspace> &h) {
                auto& _j = j["energy"];
                for (auto i : h.vec)
                    i->to_json(_j);
            }

    }//namespace
}//namespace
