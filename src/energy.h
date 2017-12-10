#pragma once

#include "core.h"
#include "geometry.h"
#include "space.h"
#include "potentials.h"
#include "multipole.h"

namespace Faunus {
    namespace Energy {

        class Energybase {
            public:
                std::string name;
                virtual double energy(Change&)=0; //!< energy due to change
                inline virtual void to_json(json &j) const {}; //!< json output
        };

        void to_json(json &j, const Energybase &base) {
            assert(!base.name.empty());
            base.to_json( j[base.name] );
        } //!< Converts any energy class to json object

        /**
         * @brief Nonbonded energy using a pair-potential
         */
        template<typename Tspace, typename Tpairpot>
            class Nonbonded : public Energybase {
                private:
                    double g2gcnt=0, g2gskip=0;
                protected:
                    double Rc2_g2g=1e20;
                    void to_json(json &j) const override {
                        j = pairpot;
                        json t = json::object();
                        t["g2g"] = { {"cutoff", std::sqrt(Rc2_g2g)} };
                        //t["cutoff_g2g"] = std::sqrt(Rc2_g2g);
                        j.push_back(t);
                    }

                    template<typename T>
                        inline bool cut(const T &g1, const T &g2) {
                            g2gcnt++;
                            if ( spc.geo.sqdist(g1.cm, g2.cm)<Rc2_g2g )
                                return false;
                            g2gskip++;
                            return true;
                        } //!< true if group<->group interaction can be skipped

                    template<typename T>
                        inline double i2i(const T &a, const T &b) {
                            return pairpot(a, b, spc.geo.vdist(a.pos, b.pos));
                        }

                    template<typename T>
                        double g_internal(const T &g) {
                            double u=0;
                            for ( auto i = g.begin(); i != g.end(); ++i )
                                for ( auto j=i; ++j != g.end(); )
                                    u += i2i(*i, *j);
                            return u;
                         }
 
                    template<typename T>
                        double g2g(const T &g1, const T &g2) {
                            double u = 0;
                            if (!cut(g1,g2))
                                for (auto &i : g1)
                                    for (auto &j : g2)
                                        u += i2i(i,j);
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

                public:
                    Tspace& spc;   //!< Space to operate on
                    Tpairpot pairpot; //!< Pair potential

                    Nonbonded(Tspace &spc, const json &j) : spc(spc) {
                        name="nonbonded";
                        pairpot = j;
                        Rc2_g2g = std::pow( j.value("cutoff_g2g", 1e20), 2);
                    }

                    double energy(Change &change) override {
                        using namespace ranges;
                        double u=0;

                        if (!change.empty()) {

                            // did everything change?
                            if (change.all) {
                                for ( auto i = spc.groups.begin(); i != spc.groups.end(); ++i ) {
                                    for ( auto j=i; ++j != spc.groups.end(); )
                                        u += g2g( *i, *j );

                                    u += g_internal(*i);
                                }
                                // more todo here...
                                return u;
                            }

                            // if exactly ONE molecule is changed
                            if (change.groups.size()==1) {
                                auto& d = change.groups[0];

                                // exactly ONE atom is changed
                                if (d.atoms.size()==1) {
                                    auto i = spc.groups[d.index].begin() + d.atoms[0];
                                    for (auto j=spc.p.begin(); j!=spc.p.end(); ++j)
                                        if (i!=j)
                                            u += i2i(*i, *j);
                                    return u;
                                }

                                // everything changed
                                if (d.all) {
                                    for (int i=0; i<int(spc.groups.size()); i++)
                                        if (i!=d.index)
                                            u+=g2g(spc.groups[i], spc.groups[d.index]);
                                    return u;
                                }
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

            }; //!< Nonbonded, pair-wise additive energy term

        template<typename Tspace>
            class Hamiltonian : public Energybase, public BasePointerVector<Energybase> {
                protected:
                    typedef typename Tspace::Tparticle Tparticle;
                    void to_json(json &j) const override {
                        for (auto i : this->vec)
                            j.push_back(*i);
                    }
                public:
                    Hamiltonian(Tspace &spc, const json &j) {
                        using namespace Potential;

                        typedef CombinedPairPotential<CoulombGalore,LennardJones<Tparticle>> CoulombLJ;
                        typedef CombinedPairPotential<CoulombGalore,HardSphere<Tparticle>> CoulombHS;

                        Energybase::name="Hamiltonian";
                        for (auto &m : j.at("energy")) {// loop over move list
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                try {
                                    if (it.key()=="nonbonded_coulomblj")
                                        push_back<Energy::Nonbonded<Tspace,CoulombLJ>>(spc, it.value());
                                    if (it.key()=="nonbonded_coulombhs")
                                        push_back<Energy::Nonbonded<Tspace,CoulombHS>>(spc, it.value());
                                    // additional energies go here...
                                } catch (std::exception &e) {
                                    throw std::runtime_error("Error adding energy '" + it.key() + "': " + e.what());
                                }
                            }
                        }
                    }

                    double energy(Change &change) override {
                        double du=0;
                        for (auto i : this->vec)
                            du += i->energy(change);
                        return du;
                    } //!< Energy due to changes

            }; //!< Aggregates and sum energy terms

    }//namespace
}//namespace
