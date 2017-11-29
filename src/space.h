#pragma once
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

            Space() {}

            Space(const json &j) {
                atoms<Tparticle> = j.at("atomlist").get<decltype(atoms<Tparticle>)>();
                molecules<Tpvec> = j.at("moleculelist").get<decltype(molecules<Tpvec>)>();
                for ( auto& mol : molecules<Tpvec> ) {
                    int n = mol.Ninit;
                    while ( n-- > 0 ) {
                        push_back(mol.id(), mol.getRandomConformation(geo, p));
                    }
                }
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
                    g.id=molid;
                    groups.push_back(g);
                    assert( in.size() == groups.back().size() );
                }
            } //!< Add particles and corresponding group to back

            auto findMolecules(int molid) {
                return groups | ranges::view::filter( [molid](auto &i){ return i.id==molid; } );
            } //!< Range with all groups of type `molid` (complexity: order N)

            auto findAtoms(int atomid) const {
                return p | ranges::view::filter( [atomid](auto &i){ return i.id==atomid; } );
            } //!< Range with all atoms of type `atomid` (complexity: order N)

            void sync(const Tspace &o, const Tchange &change) {

                for (auto &m : change.groups) {

                    auto &go = groups.at(m.index);  // old group
                    auto &gn = o.groups.at(m.index);// new group

                    go = gn; // sync group (*not* the actual elements)

                    assert( gn.size() == go.size() );
                    //assert( go.size()
                    //        < std::max_element(m.atoms.begin(), m.atoms.end()));

                    if (m.all) // all atoms have moved
                        std::copy(gn.begin(), gn.end(), go.begin() );
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
#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Space")
    {
        typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
        typedef Space<Geometry::Cuboid, Tparticle> Tspace;
        Tspace spc;
        spc.p.resize(10);
    }
#endif

}//namespace
