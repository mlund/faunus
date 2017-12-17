#pragma once
#include "core.h"
#include "molecule.h"
#include "geometry.h"
#include "group.h"

namespace Faunus {

    /**
     * @brief Specify change to a new state
     *
     * - If `moved` or `removed` are defined for a group, but are
     *   empty, it is assumed that *all* particles in the group are affected.
     */
    struct Change {
        bool dV = false;    //!< Set to true if there's a volume change
        double all = false; //!< Set to true if *everything* has changed

        struct data {
            int index; //!< Touched group index
            bool all=false; //!< Set to `true` if all particles in group have been updated
            std::vector<int> atoms; //!< Touched atom index w. respect to `Group::begin()`
        }; //!< Properties of changed groups

        std::vector<data> groups; //!< Touched groups by index in group vector

        auto touchedGroupIndex() {
            return ranges::view::transform(groups, [](data &i) -> int {return i.index;});
        } //!< List of moved groups (index)

        void clear()
        {
            dV=false;;
            all=false;
            groups.clear();
            assert(empty());
        } //!< Clear all change data

        bool empty() const
        {
            if (all==false)
                if (dV==false)
                    if (groups.empty())
                        return true;
            return false;
        } //!< Check if change object is empty

    };

    template<class Tgeometry, class Tparticletype>
        struct Space {

            typedef Tparticletype Tparticle;
            typedef Space<Tgeometry,Tparticle> Tspace;
            typedef std::vector<Tparticle> Tpvec;
            typedef Group<Tparticle> Tgroup;
            typedef std::vector<Tgroup> Tgvec;
            typedef Change Tchange;

            typedef std::function<void(Tspace&, const Tchange&)> ChangeTrigger;
            typedef std::function<void(Tspace&, const Tspace&, const Tchange&)> SyncTrigger;

            std::vector<ChangeTrigger> changeTriggers; //!< Call when a Change object is applied
            std::vector<SyncTrigger> onSyncTriggers;   //!< Call when two Space objects are synched

            Tpvec p;       //!< Particle vector
            Tgvec groups;  //!< Group vector
            Tgeometry geo; //!< Container geometry

            //std::vector<MoleculeData<Tpvec>>& molecules; //!< Reference to global molecule list

            Space() {}

            Space(const json &j) {
                new_from_json(j, *this);

                for (auto &i : groups)
                    if (!i.empty())
                        if (!i.atomic)
                            if (geo.sqdist( i.cm,
                                        Geometry::massCenter(i.begin(), i.end(), geo.boundaryFunc, -i.cm) ) > 1e-9 )
                                throw std::runtime_error("space construction error: mass center mismatch");
            }

            void clear() {
                p.clear();
                groups.clear();
            }

            void push_back(int molid, const Tpvec &in) {
                if (!in.empty()) {
                    auto oldbegin = p.begin();
                    p.insert( p.end(), in.begin(), in.end() );
                    if (p.begin() != oldbegin) {// update group iterators if `p` is relocated
                        for (auto &g : groups) {
                            g.relocate(oldbegin, p.begin());
                        }
                    }
                    Tgroup g( p.end()-in.size(), p.end() );
                    g.id = molid;
                    g.atomic = molecules<Tpvec>.at(molid).atomic;
                    g.cm = Geometry::massCenter(in.begin(), in.end(), geo.boundaryFunc, -in.begin()->pos);

                    if (g.atomic==false) {
                        Point cm = Geometry::massCenter(g.begin(), g.end(), geo.boundaryFunc, -g.cm);
                        if (geo.sqdist(g.cm, cm)>1e-9)
                            throw std::runtime_error("space: mass center error upon insertion. Molecule too large?\n");
                    }

                    // inactivate group before insertion?
                    auto molvec = findInactiveMolecules(molid);
                    if (size(molvec) < molecules<Tpvec>.at(molid).Ninactive)
                        g.resize(0);

                    groups.push_back(g);
                    assert( in.size() == groups.back().capacity() );
                }
            } //!< Safely add particles and corresponding group to back

            auto findMolecules(int molid) {
                return groups | ranges::view::filter( [molid](auto &i){ return i.id==molid; } );
            } //!< Range with all groups of type `molid` (complexity: order N)

            auto findInactiveMolecules(int molid) {
                return groups | ranges::view::filter( [molid](auto &i){ return (i.id==molid) && (i.size()!=i.capacity()); } );
            } //!< Range with all (partially) inactive groups of type `molid` (complexity: order N)

            auto findAtoms(int atomid) const {
                return p | ranges::view::filter( [atomid](auto &i){ return i.id==atomid; } );
            } //!< Range with all atoms of type `atomid` (complexity: order N)

            void sync(Tspace &other, const Tchange &change) {

                assert(&other != this);
                assert( p.begin() != other.p.begin());

                if (change.dV || change.all) {
                    geo = other.geo;
                }

                // deep copy *everything*
                if (change.all) {
                    p = other.p; // copy all positions
                    assert( p.begin() != other.p.begin() && "deep copy problem");
                    groups = other.groups;

                    if (!groups.empty())
                        if (groups.front().begin() == other.p.begin())
                            for (auto &i : groups)
                                i.relocate( other.p.begin(), p.begin() );
                }
                else {
                    for (auto &m : change.groups) {

                        auto &g = groups.at(m.index);  // old group
                        auto &gother = other.groups.at(m.index);// new group

                        g.shallowcopy(gother); // copy group data but *not* particles

                        if (m.all) // copy all particles 
                            std::copy( gother.begin(), gother.end(), g.begin() );
                        else // copy only a subset
                            for (auto i : m.atoms)
                                *(g.begin()+i) = *(gother.begin()+i);
                    }
                }
                assert( p.size() == other.p.size() );
                assert( p.begin() != other.p.begin());
            } //!< Copy differing data from other (o) Space using Change object

            json info() {
                json j;
                j["number of particles"] = p.size();
                j["number of groups"] = groups.size();
                j["geometry"] = geo;
                auto& _j = j["groups"];
                for (auto &i : groups) {
                    auto name = molecules<decltype(p)>.at(i.id).name;
                    json tmp, d=i;
                    d.erase("cm");
                    d.erase("id");
                    d.erase("atomic");
                    auto ndx = i.to_index(p.begin());
                    d["index"] = std::to_string(ndx.first)+"-"+std::to_string(ndx.second);
                    tmp[name] = d;
                    _j.push_back( tmp );
                }
                return j;
            }

        };

    template<class Tgeometry, class Tparticle>
        void to_json(json &j, Space<Tgeometry,Tparticle> &spc) {
            j["temperature"] = pc::temperature;
            j["atomlist"] = atoms<Tparticle>;
            j["moleculelist"] = molecules<decltype(spc.p)>;
            j["geometry"] = spc.geo;
            j["groups"] = spc.groups;
            j["particles"] = spc.p;
        } //!< Serialize Space to json object

    template<class Tgeometry, class Tparticletype>
        void new_from_json(const json &j, Space<Tgeometry,Tparticletype> &spc) {
            typedef typename Space<Tgeometry,Tparticletype>::Tpvec Tpvec;
            spc.clear();

            pc::temperature = j.at("temperature").get<double>() * 1.0_K;
            spc.geo = j.at("geometry");
            atoms<Tparticletype> = j.at("atomlist").get<decltype(atoms<Tparticletype>)>();
            molecules<Tpvec> = j.at("moleculelist").get<decltype(molecules<Tpvec>)>();

            if ( j.count("groups")==0 )
            {
                insertMolecules( spc );
            }
            else
            {
                spc.p = j.at("particles").get<Tpvec>();

                if (!spc.p.empty()) {
                    auto begin = spc.p.begin();
                    Group<Tparticletype> g(begin,begin);
                    for (auto &i : j.at("groups")) {
                        g.begin() = begin;
                        from_json(i, g);
                        spc.groups.push_back(g);
                        begin = g.trueend();
                    }
                    if (begin != spc.p.end())
                        throw std::runtime_error("json->space load error");
                }
            }
        } //!< Deserialize json object to Space

    template<typename Tspace>
        void insertMolecules(Tspace &spc) {
            spc.clear();
            assert(spc.geo.getVolume()>0);
            for ( auto& mol : molecules<typename Tspace::Tpvec> ) {
                if (mol.atomic)
                    spc.push_back(mol.id(), mol.getRandomConformation(spc.geo, spc.p));
                else {
                    int n = mol.Ninit;
                    while ( n-- > 0 ) {
                        spc.push_back(mol.id(), mol.getRandomConformation(spc.geo, spc.p));
                    }
                }
            }
        } //!< Insert `Ninit` molecules into space as defined in `molecules`

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Space")
    {
        typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
        typedef Space<Geometry::Cuboid, Tparticle> Tspace;
        Tspace spc1;

        // check molecule insertion
        atoms<typename Tspace::Tpvec>.resize(2);
        CHECK( atoms<typename Tspace::Tpvec>.at(0).mw == 1);
        Tparticle a;
        a.id=0;
        Tspace::Tpvec p(2, a);
        CHECK( p[0].mw == 1);
        p[0].pos.x()=2;
        p[1].pos.x()=3;
        spc1.push_back(1, p);
        CHECK( spc1.p.size()==2 );
        CHECK( spc1.groups.size()==1 );
        CHECK( spc1.groups.front().id==1);
        CHECK( spc1.groups.front().cm.x()==doctest::Approx(2.5));

        // sync groups
        Change c;
        c.all=true;
        c.dV=true;
        c.groups.resize(1);
        c.groups[0].index=0;
        c.groups[0].all=true;
        Tspace spc2;
        spc2.sync(spc1, c);
        CHECK( spc2.p.size()==2 );
        CHECK( spc2.groups.size()==1 );
        CHECK( spc2.groups.front().id==1);
        CHECK( spc2.groups.front().begin() != spc1.groups.front().begin());
        CHECK( spc2.p.front().pos.x()==doctest::Approx(2));

        // nothing should be synched (all==false)
        spc2.p.back().pos.z()=-0.1;
        c.all=false;
        c.groups[0].all=false;
        spc1.sync(spc2, c);
        CHECK( spc1.p.back().pos.z() != -0.1 );

        // everything should be synched (all==true)
        c.groups[0].all=true;
        spc1.sync(spc2, c);
        CHECK( spc1.p.back().pos.z() == doctest::Approx(-0.1) );
    }
#endif

}//namespace
