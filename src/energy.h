#pragma once

#include "core.h"
#include "geometry.h"
#include "space.h"
#include "potentials.h"
#include "multipole.h"
#include <set>

namespace Faunus {
    namespace Energy {

        class Energybase {
            public:
                std::string name;
                std::string cite;
                virtual double energy(Change&)=0; //!< energy due to change
                inline virtual void to_json(json &j) const {}; //!< json output
        };

        void to_json(json &j, const Energybase &base) {
            assert(!base.name.empty());
            if (!base.cite.empty())
                j[base.name]["reference"] = base.cite;
            base.to_json( j[base.name] );
        } //!< Converts any energy class to json object

        template<typename Tspace>
            class Isobaric : public Energybase {
                private:
                    Tspace& spc;
                    double P; // P/kT
                public:
                    Isobaric(const json &j, Tspace &spc) : spc(spc) {
                        name = "isobaric";
                        cite = "Frenkel & Smith 2nd Ed (Eq. 5.4.13)";
                        P = j.value("P/mM", 0.0) * 1.0_mM;
                        if (P<1e-10) {
                            P = j.value("P/Pa", 0.0) * 1.0_Pa;
                            if (P<1e-10)
                                P = j.at("P/atm").get<double>() * 1.0_atm;
                        }
                    }
                    double energy(Change &change) override {
                        if (change.dV || change.all) {
                            double V = spc.geo.getVolume();
                            size_t N=0;
                            for (auto &g : spc.groups)
                                if (!g.empty()) {
                                    if (g.atomic)
                                        N += g.size();
                                    else
                                        N++;
                                }
                            return P*V-(N+1)*std::log(V);
                        } else return 0; 
                    }
                    void to_json(json &j) const override {
                        j["P/atm"] = P / 1.0_atm;
                        j["P/mM"] = P / 1.0_mM;
                        j["P/Pa"] = P / 1.0_Pa;
                        _roundjson(j,5);
                    }
            };

        template<typename Tspace>
            class ExternalPotential : public Energybase {
                protected:
                    typedef typename Tspace::Tpvec Tpvec;
                    typedef typename Tspace::Tparticle Tparticle;
                    Tspace& spc;
                    std::set<int> molids; // molecules to act upon
                    std::function<double(const Tparticle&)> func=nullptr; // energy of single particle

                    template<class Tparticle>
                        double _energy(const Group<Tparticle> &g) const {
                            double u=0;
                            if (molids.find(g.id)!=molids.end())
                                for (auto &p : g) {
                                    u += func(p);
                                    if (std::isnan(u))
                                        break;
                                }
                            return u;
                        } //!< External potential on a single particle
                public:
                    ExternalPotential(const json &j, Tspace &spc) : spc(spc) {
                        name="external";
                        std::vector<std::string> _names = j.at("molecules"); // molecule names
                        auto _ids = names2ids(molecules<Tpvec>, _names);     // names --> molids
                        molids = std::set<int>(_ids.begin(), _ids.end());    // vector --> set
                        if (molids.empty() || molids.size()!=_names.size() )
                            throw std::runtime_error(name + ": molecule list is empty");

                    }

                    double energy(Change &change) override {
                        assert(func!=nullptr);
                        double u=0;
                        if (change.dV || change.all) {
                            for (auto &g : spc.groups) { // check all groups
                                u += _energy(g);
                                if (std::isnan(u))
                                    break;
                            }
                        } else
                            for (auto &d : change.groups) {
                                auto &g = spc.groups.at(d.index); // check specified groups
                                if (d.all)  // check all atoms in group
                                    u += _energy(g);
                                else       // check only specified atoms in group
                                    for (auto i : d.atoms)
                                        u += func( *(g.begin()+i) );
                                if (std::isnan(u))
                                    break;
                            }
                        return u;
                    }

                    void to_json(json &j) const override {
                        _roundjson(j,5);
                    }
            }; //!< Base class for external potentials, acting on particles

        template<typename Tspace, typename base=ExternalPotential<Tspace>>
            class Confine : public base {
                public:
                    enum Variant {sphere, cylinder, cuboid, none};
                    Variant type=none;

                private:
                    Point origo={0,0,0}, dir={1,1,1};
                    double radius, k;
                    std::map<std::string, Variant> m = {
                        {"sphere", sphere}, {"cylinder", cylinder}
                    };

                public:
                    Confine(const json &j, Tspace &spc) : base(j,spc) {
                        base::name = "confine";
                        k = value_inf(j, "k") * 1.0_kJmol; // get floating point; allow inf/-inf
                        type = m.at( j.at("type") );

                        if (type==sphere || type==cylinder) {
                            radius = j.at("radius");
                            origo = j.value("origo", origo);
                            if (type==cylinder)
                                dir = {1,1,0};
                            base::func = [radius=radius, origo=origo, k=k, dir=dir](const typename base::Tparticle &p) {
                                double d2 = (origo-p.pos).cwiseProduct(dir).squaredNorm() - radius*radius;
                                if (d2>0)
                                    return 0.5*k*d2;
                                return 0.0;
                            };
                        }
                    }
            }; //!< Confine particles to a sub-region of the simulation container

        template<typename Tspace>
            class Bonded : public Energybase {
                private:
                    typedef typename Tspace::Tpvec Tpvec;
                    Tspace& spc;
                    typedef std::vector<Potential::BondData> BondVector;
                    BondVector inter;  // inter-molecular bonds
                    std::map<int,BondVector> intra; // intra-molecular bonds

                    void update() {
                        for (size_t i=0; i<spc.groups.size(); i++) {
                            auto &g = spc.groups[i];
                            intra[i] = molecules<Tpvec>.at(g.id).bonds;
                            for (auto &b : intra[i])
                                b.shift( std::distance(spc.p.begin(), g.begin()) );
                        }
                    }

                    double sum( const BondVector &v ) const {
                        double u=0;
                        for (auto &b : v)
                            u += b.energy(spc.p, spc.geo.distanceFunc);
                        return u;
                    }

                public:
                    Bonded(const json &j, Tspace &spc) : spc(spc) {
                        name = "bonded";
                        update();
                        if (j.is_object())
                            if (j.count("bondlist")==1)
                                inter = j["bondlist"].get<BondVector>();
                    }
                    void to_json(json &j) const override {
                        if (!inter.empty())
                            j["bondlist"] = inter;
                        if (!intra.empty()) {
                            json& _j = j["bondlist-intramolecular"];
                            _j = json::array();
                            for (auto &i : intra)
                                for (auto &b : i.second)
                                    _j.push_back(b);
                        }
                    }

                    double energy(Change &c) override {
                        double u=0;
                        if ( !c.empty() ) {
                            u = sum(inter);
                            if ( c.all || c.dV )
                                for (auto& i : intra)
                                    u += sum(i.second);
                            else
                                for (auto &d : c.groups)
                                    u += sum( intra[d.index] );
                        }
                        return u;
                    }; // brute force -- refine this!
            };

        /**
         * @brief Nonbonded energy using a pair-potential
         */
        template<typename Tspace, typename Tpairpot>
            class Nonbonded : public Energybase {
                private:
                    double g2gcnt=0, g2gskip=0;
                protected:
                    double Rc2_g2g=pc::infty;
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
                        double g2all(const T &g1) {
                            double u = 0;
                            for ( auto i = spc.groups.begin(); i != spc.groups.end(); ++i ) {
                                for ( auto j=i; ++j != spc.groups.end(); )
                                    u += g2g( *i, *j );
                                return u;
                            }
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

                    Nonbonded(const json &j, Tspace &spc) : spc(spc) {
                        name="nonbonded";
                        pairpot = j;
                        Rc2_g2g = std::pow( j.value("cutoff_g2g", pc::infty), 2);
                    }

                    double energy(Change &change) override {
                        using namespace ranges;
                        double u=0;

                        if (!change.empty()) {

                            if (change.dV) {
#pragma omp parallel for reduction (+:u) schedule (dynamic) 
                                for ( auto i = spc.groups.begin(); i < spc.groups.end(); ++i ) {
                                    for ( auto j=i; ++j != spc.groups.end(); )
                                        u += g2g( *i, *j );
                                    if (i->atomic)
                                        u += g_internal(*i);
                                }
                                return u;
                            }

                            // did everything change?
                            if (change.all) {
#pragma omp parallel for reduction (+:u) schedule (dynamic) 
                                for ( auto i = spc.groups.begin(); i < spc.groups.end(); ++i ) {
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
                                // WARNING! This does not respect inactive particles!
                                // TODO: Loop over groups instead
                                if (d.atoms.size()==1) {
                                    auto i = spc.groups[d.index].begin() + d.atoms[0];
                                    for (auto j=spc.p.begin(); j!=spc.p.end(); ++j)
                                        if (i!=j)
                                            u += i2i(*i, *j);
                                    return u;
                                }

                                // everything in group changed
                                if (d.all) {
#pragma omp parallel for reduction (+:u) 
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
                        typedef CombinedPairPotential<CoulombGalore,WeeksChandlerAndersen<Tparticle>> CoulombWCA;

                        Energybase::name="hamiltonian";
                        for (auto &m : j.at("energy")) {// loop over move list
                            size_t oldsize = vec.size();
                            for (auto it=m.begin(); it!=m.end(); ++it) {
                                try {
                                    if (it.key()=="nonbonded_coulomblj")
                                        push_back<Energy::Nonbonded<Tspace,CoulombLJ>>(it.value(), spc);

                                    if (it.key()=="nonbonded_coulombhs")
                                        push_back<Energy::Nonbonded<Tspace,CoulombHS>>(it.value(), spc);

                                    if (it.key()=="nonbonded_coulombwca")
                                        push_back<Energy::Nonbonded<Tspace,CoulombWCA>>(it.value(), spc);

                                    if (it.key()=="isobaric")
                                        push_back<Energy::Isobaric<Tspace>>(it.value(), spc);

                                    if (it.key()=="confine")
                                        push_back<Energy::Confine<Tspace>>(it.value(), spc);

                                    if (it.key()=="bonded")
                                        push_back<Energy::Bonded<Tspace>>(it.value(), spc);

                                    // additional energies go here...

                                    if (vec.size()==oldsize)
                                        std::cerr << "warning: ignoring unknown energy '" << it.key() << "'" << endl;

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
