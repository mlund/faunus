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
        double dV = 0;     //!< Volume change (in different directions)

        struct data {
            int index; //!< Touched group index
            bool all=false; //!< Set to `true` if all particles in group have been updated
            std::vector<int> atoms; //!< Touched atom index w. respect to `Group::begin()`
            std::vector<std::pair<int,int>> activated, deactivated; //!< Range of (de)activated particles
        }; //!< Properties of changed groups

        std::vector<data> groups; //!< Touched groups by index in group vector

        auto touchedGroupIndex() {
            return ranges::view::transform(groups, [](data &i) -> int {return i.index;});
        } //!< List of moved groups (index)

        void clear()
        {
            dV=0;
            groups.clear();
            assert(empty());
        } //!< Clear all change data

        bool empty() const
        {
            if ( groups.empty())
                if ( dV==0 )
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
                    g.cm = Geometry::massCenter(in.begin(), in.end(), geo.boundaryFunc);
                    groups.push_back(g);
                    assert( in.size() == groups.back().size() );
                }
            } //!< Safely add particles and corresponding group to back

            auto findMolecules(int molid) {
                return groups | ranges::view::filter( [molid](auto &i){ return i.id==molid; } );
            } //!< Range with all groups of type `molid` (complexity: order N)

            auto findAtoms(int atomid) const {
                return p | ranges::view::filter( [atomid](auto &i){ return i.id==atomid; } );
            } //!< Range with all atoms of type `atomid` (complexity: order N)

            void sync(Tspace &o, const Tchange &change) {

                if (std::fabs(change.dV)>1e-9)
                    geo = o.geo;

                // if mismatch do a deep copy of everything
                if ( (p.size()!=o.p.size()) || (groups.size()!=o.groups.size())) {
                    p = o.p;
                    groups = o.groups;
                    for (auto &i : groups)
                        i.relocate( o.p.begin(), p.begin() );
                    return;
                }

                assert( p.begin() != o.p.begin());

                for (auto &m : change.groups) {

                    auto &go = groups.at(m.index);  // old group
                    auto &gn = o.groups.at(m.index);// new group

                    //go = gn; // deep copy group
                    //assert( gn.size() == go.size() );
                    assert( gn.capacity() == go.capacity() );

                    if (m.all) // all atoms have moved
                        std::copy( gn.begin(), gn.end(), go.begin() );
                    else // only some atoms have moved
                        for (auto i : m.atoms)
                            *(go.begin()+i) = *(gn.begin()+i);
                }
            } //!< Copy differing data from other (o) Space using Change object

            void applyChange(const Tchange &change) {
                for (auto& f : changeTriggers)
                    f(*this, change);
            }
        };

        template<class Tgeometry, class Tparticletype>
        void from_json(const json &j, Space<Tgeometry,Tparticletype> &p) {
            assert(1==2 && "to be implemented");
        }

        template<typename Tspace>
        void insertMolecules(Tspace &spc) {
            spc.clear();
            for ( auto& mol : molecules<typename Tspace::Tpvec> ) {
                int n = mol.Ninit;
                while ( n-- > 0 ) {
                    spc.push_back(mol.id(), mol.getRandomConformation(spc.geo, spc.p));
                }
            }
        } //!< Insert `Ninit` molecules into space as defined in `molecules`

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Space")
    {
        typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
        typedef Space<Geometry::Cuboid, Tparticle> Tspace;
        Tspace spc1, spc2;

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
        c.groups.resize(1);
        c.groups[0].index=0;
        c.groups[0].all=true;
        spc2.sync(spc1, c);
        CHECK( spc2.p.size()==2 );
        CHECK( spc2.groups.size()==1 );
        CHECK( spc2.groups.front().begin() != spc1.groups.front().begin());
        CHECK( spc2.p.front().pos.x()==doctest::Approx(2));

        // nothing should be synched (all==false)
        spc2.p.back().pos.z()=-0.1;
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
