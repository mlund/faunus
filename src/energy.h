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
                virtual double energy(Change&)=0; //!< energy due to change
                virtual void to_json(json &j) const=0 ; //!< json output
        };

        void to_json(json &j, const Energybase &base) {
            assert(!base.name.empty());
            base.to_json( j[base.name] );
        } //!< Converts any energy class to json object

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
                    using namespace ranges;
                    double u=0;

                    if (!change.empty()) {

                        // did everything change?
                        if (change.all) {
                            // all groups<->all groups
                            //for (size_t i=0; i<spc.groups.size()-1; i++)
                            //    for (size_t j=i+1; j<spc.groups.size(); j++) {
                            //        u+= g2g(spc.groups.at(i), spc.groups.at(j));
                            //    }
                            for ( auto i = spc.groups.begin(); i != spc.groups.end(); ++i ) {
                                for ( auto j=i; ++j != spc.groups.end(); ) {
                                    u += g2g( *i, *j );
                                }
                            }
                            return u;
                        }

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

                void to_json(json &j) const override { j = pairpot; }

            }; //!< Nonbonded, pair-wise additive energy term

        template<typename Tspace>
            class Hamiltonian : public Energybase, public BasePointerVector<Energybase> {
                public:
                    typedef typename Tspace::Tparticle Tparticle;

                    Hamiltonian(Tspace &spc, const json &j) {
                        using namespace Potential;
                        typedef CombinedPairPotential<Coulomb,LennardJones<Tparticle>> Tpairpot;

                        Energybase::name="Hamiltonian";
                        for (auto &m : j.at("energy")) {// loop over move list
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                if (it.key()=="nonbonded_lj")
                                    push_back<Energy::Nonbonded<Tspace,LennardJones<Tparticle>>>(spc, it.value());
                                // additional energies go here...
                            }
                        }
                    }

                    double energy(Change &change) override {
                        double du=0;
                        for (auto i : this->vec)
                            du += i->energy(change);
                        return du;
                    } //!< Energy due to changes

                    void to_json(json &j) const override {
                        for (auto i : this->vec)
                            j.push_back(*i);
                    }
            }; //!< Aggregates and sum energy terms

    }//namespace
}//namespace
