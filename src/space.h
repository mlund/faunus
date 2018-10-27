#pragma once
#include "core.h"
#include "geometry.h"
#include "group.h"
#include "molecule.h"

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
        double du=0;        //!< Additional energy change not captured by Hamiltonian
        bool dNpart=false;      //!< Is the size of groups partially changed

        struct data {
            bool dNpart=false;      //!< Is the size of groups partially changed
            int index; //!< Touched group index
            bool internal=false; //!< True is the internal energy/config has changed
            bool all=false; //!< Set to `true` if all particles in group have been updated
            std::vector<int> atoms; //!< Touched atom index w. respect to `Group::begin()`

            inline bool operator<( const data & a ) const{
                return index < a.index;
            }
        }; //!< Properties of changed groups

        std::vector<data> groups; //!< Touched groups by index in group vector

        auto touchedGroupIndex() {
            // skitfunktion
            return ranges::view::transform(groups, [](data &i) -> int {return i.index;});
        } //!< List of moved groups (index)

        inline void clear()
        {
            du=0;
            dV=false;
            all=false;
            dNpart=false;
            groups.clear();
            assert(empty());
        } //!< Clear all change data

        inline bool empty() const
        {
            if (du==0)
                if (dV==false)
                    if (all==false)
                        if (groups.empty())
                            if (dNpart==false)
                                return true;
            return false;
        } //!< Check if change object is empty

    };

    struct SpaceBase {
        virtual Geometry::GeometryBase& getGeometry()=0;
        virtual void clear()=0;
        virtual void scaleVolume(double Vnew, Geometry::VolumeMethod method=Geometry::ISOTROPIC)=0;
    };

    template<class Tgeometry, class Tparticletype>
        struct Space {

            typedef Tparticletype Tparticle;
            typedef Space<Tgeometry,Tparticle> Tspace;
            typedef std::vector<Tparticle> Tpvec;
            typedef Group<Tparticle> Tgroup;
            typedef std::vector<Tgroup> Tgvec;
            typedef Change Tchange;

            typedef std::function<void(Tspace&, double, double)> ScaleVolumeTrigger;
            typedef std::function<void(Tspace&, const Tchange&)> ChangeTrigger;
            typedef std::function<void(Tspace&, const Tspace&, const Tchange&)> SyncTrigger;

            std::vector<ScaleVolumeTrigger> scaleVolumeTriggers; //!< Call when volume is scaled
            std::vector<ChangeTrigger> changeTriggers; //!< Call when a Change object is applied
            std::vector<SyncTrigger> onSyncTriggers;   //!< Call when two Space objects are synched

            Tpvec p;       //!< Particle vector
            Tgvec groups;  //!< Group vector
            Tgeometry geo; //!< Container geometry

            auto positions() const {
               return ranges::view::transform(p, [](auto &i) -> const Point& {return i.pos;});
            } //!< Iterable range with positions 

            enum Selection {ALL, ACTIVE, INACTIVE};

            void clear() {
                p.clear();
                groups.clear();
            } //!< Clears particle and molecule list

            /*
             * The following is considered:
             *
             * - `groups` vector is expanded with a new group at the end
             * - if the particle vector is relocated, all existing group
             *   iterators are updated to reflect the new memory positions
             */
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

                    if (g.atomic==false) {
                        g.cm = Geometry::massCenter(in.begin(), in.end(), geo.boundaryFunc, -in.begin()->pos);
                        Point cm = Geometry::massCenter(g.begin(), g.end(), geo.boundaryFunc, -g.cm);
                        if (geo.sqdist(g.cm, cm)>1e-9)
                            throw std::runtime_error("space: mass center error upon insertion. Molecule too large?\n");
                    }

                    groups.push_back(g);
                    assert( in.size() == groups.back().capacity() );
                }
            } //!< Safely add particles and corresponding group to back

            auto findMolecules(int molid, Selection sel=ACTIVE) {
                std::function<bool(Tgroup&)> f;
                switch (sel) {
                    case (ALL):
                        f = [molid](Tgroup &i){ return i.id==molid; };
                        break;
                    case (INACTIVE):
                        f = [molid](Tgroup &i){ return (i.id==molid) && (i.size()!=i.capacity()); };
                        break;
                    case (ACTIVE):
                        f = [molid](Tgroup &i){ return (i.id==molid) && (i.size()==i.capacity()); };
                        break;
                }
                return groups | ranges::view::filter(f);
            } //!< Range with all groups of type `molid` (complexity: order N)

            typename decltype(groups)::iterator randomMolecule(int molid, Random &rand, Selection sel=ACTIVE) {
                auto m = findMolecules(molid, sel);
                if (size(m)>0)
                    return groups.begin() + (&*rand.sample( m.begin(), m.end() ) - &*groups.begin());
                return groups.end();
            } //!< Random group; groups.end() if not found

            auto findAtoms(int atomid) const {
                return p | ranges::view::filter( [atomid](auto &i){ return i.id==atomid; } );
            } //!< Range with all atoms of type `atomid` (complexity: order N)

            auto findGroupContaining(const Tparticle &i) {
                return std::find_if( groups.begin(), groups.end(), [&i](auto &g){ return g.contains(i); });
            } //!< Finds the groups containing the given atom

            void sync(Tspace &other, const Tchange &change) {

                assert(&other != this);
                assert( p.begin() != other.p.begin());

                if (change.dV or change.all)
                    geo = other.geo;

                // deep copy *everything*
                if (change.all) {
                    p = other.p; // copy all positions
                    assert( p.begin() != other.p.begin() && "deep copy problem");
                    groups = other.groups;

                    if (not groups.empty())
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

            /*
             * Scales:
             * - positions of free atoms
             * - positions of molecular masscenters
             * - simulation container
             */
            void scaleVolume(double Vnew, Geometry::VolumeMethod method=Geometry::ISOTROPIC) {
                for (auto &g: groups) // remove periodic boundaries
                    if (!g.atomic)
                        g.unwrap(geo.distanceFunc);

                Point scale = geo.setVolume(Vnew, method);

                for (auto& g : groups) {
                    if (!g.empty()) {
                        if (g.atomic) // scale all atoms
                            for (auto& i : g)
                                i.pos = i.pos.cwiseProduct(scale);
                        else { // scale mass center and translate
                            Point delta = g.cm.cwiseProduct(scale) - g.cm;
                            g.cm = g.cm.cwiseProduct(scale);
                            for (auto &i : g) {
                                i.pos += delta;
                                geo.boundary(i.pos);
                            }
                            assert( geo.sqdist( g.cm,
                                        Geometry::massCenter(
                                            g.begin(), g.end(),
                                            geo.boundaryFunc, -g.cm)) < 1e-10 );
                        }
                    }
                }
                double Vold = geo.getVolume();
                // if isochoric, the volume is constant
                if (method==Geometry::ISOCHORIC)
                    Vold = std::pow(Vold,1./3.);

                for (auto f : scaleVolumeTriggers)
                    f(*this, Vold, Vnew);
            } //!< scale space to new volume


            json info() {
                json j = {
                    {"number of particles", p.size()},
                    {"number of groups", groups.size()},
                    {"geometry", geo}
                };
                auto& _j = j["groups"];
                for (auto &i : groups) {
                    auto& name = molecules<decltype(p)>.at(i.id).name;
                    json tmp, d=i;
                    d.erase("cm");
                    d.erase("id");
                    d.erase("atomic");
                    auto ndx = i.to_index(p.begin());
                    if (not i.empty())
                        d["index"] = { ndx.first, ndx.second };
                        //d["index"] = std::to_string(ndx.first)+"-"+std::to_string(ndx.second);
                    tmp[name] = d;
                    _j.push_back( tmp );
                }
                auto& _j2 = j["reactionlist"];
                for (auto &i : reactions<decltype(p)>) {
                    json tmp, d = i;
                    _j2.push_back( i );
                }
                return j;
            }

        }; //!< Space for particles, groups, geometry

    template<class Tgeometry, class Tparticle>
        void to_json(json &j, Space<Tgeometry,Tparticle> &spc) {
            typedef typename Space<Tgeometry,Tparticle>::Tpvec Tpvec;
            j["geometry"] = spc.geo;
            j["groups"] = spc.groups;
            j["particles"] = spc.p;
            j["reactionlist"] = reactions<Tpvec>;
        } //!< Serialize Space to json object

    template<class Tgeometry, class Tparticletype>
        void from_json(const json &j, Space<Tgeometry,Tparticletype> &spc) {
            typedef typename Space<Tgeometry,Tparticletype>::Tpvec Tpvec;
            using namespace std::string_literals;

            try {
                if (atoms.empty())
                    atoms = j.at("atomlist").get<decltype(atoms)>();
                if (molecules<Tpvec>.empty())
                    molecules<Tpvec> = j.at("moleculelist").get<decltype(molecules<Tpvec>)>();
                if (reactions<Tpvec>.empty())
                    if (j.count("reactionlist")>0)
                        reactions<Tpvec> = j.at("reactionlist").get<decltype(reactions<Tpvec>)>();

                spc.clear();
                spc.geo = j.at("geometry");

                if ( j.count("groups")==0 ) {
                    insertMolecules( j.at("insertmolecules"), spc );
                } else {
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
                            throw std::runtime_error("load error");
                    }
                }
                // check correctness of molecular mass centers
                for (auto &i : spc.groups)
                    if (!i.empty())
                        if (!i.atomic)
                            if (spc.geo.sqdist( i.cm,
                                        Geometry::massCenter(i.begin(), i.end(), spc.geo.boundaryFunc, -i.cm) ) > 1e-9 )
                                throw std::runtime_error("mass center mismatch");
            } catch(std::exception& e) {
                throw std::runtime_error("Error while constructing Space from JSON"s + e.what());
            }

        } //!< Deserialize json object to Space

    /**
     * @brief Insert molecules into space
     *
     * Expects an array of molecules, for example:
     *
     * ~~~ yaml
     *     - salt:  { N: 10 }
     *     - water: { N: 256 }
     *     - water: { N: 1, inactive: true }
     * ~~~
     */
    template<typename Tspace>
        void insertMolecules(const json &j, Tspace &spc) {
            typedef typename Tspace::Tpvec Tpvec;
            spc.clear();
            assert(spc.geo.getVolume()>0);
            auto &molvec = molecules<Tpvec>;
            if (j.is_array()) {
                for (auto &m : j) { // loop over array of molecules
                    if (m.is_object() && m.size()==1)
                        for (auto it=m.begin(); it!=m.end(); ++it) {
                            auto mol = findName(molvec, it.key()); // is the molecule defined?
                            if (mol!=molvec.end()) {

                                int N = it.value().at("N").get<int>();  // number of molecules to insert
                                int cnt = N;
                                bool inactive = it.value().value("inactive", false); // active or not?

                                if (mol->atomic) {
                                    typename Tspace::Tpvec p;
                                    p.reserve( N * mol->atoms.size() );
                                    while ( cnt-- > 0 ) {
                                        auto _t = mol->getRandomConformation(spc.geo, spc.p);
                                        p.insert(p.end(), _t.begin(), _t.end());
                                    }
                                    assert(!p.empty());
                                    spc.push_back(mol->id(), p);
                                    if (inactive)
                                        spc.groups.back().resize(0);
                                } else {
                                    while ( cnt-- > 0 ) { // insert molecules
                                        spc.push_back(mol->id(), mol->getRandomConformation(spc.geo, spc.p));
                                        if (inactive)
                                            spc.groups.back().resize(0);
                                    }
                                    // load specific positions for the N added molecules
                                    bool success=false;
                                    std::string file = it.value().value("positions", "");
                                    if (!file.empty()) {
                                        Tpvec p;
                                        if (loadStructure<Tpvec>()(file, p, false)) {
                                            if (p.size() == N*mol->atoms.size()) {
                                                Point offset = it.value().value("translate", Point(0,0,0));
                                                size_t j = spc.p.size() - p.size();
                                                for (auto &i : p) {
                                                    i.pos = i.pos + offset;
                                                    if (spc.geo.collision(i.pos)==false)
                                                        spc.p.at(j++).pos = i.pos;
                                                    else std::cerr << "position outside box" << endl;
                                                }
                                                if (j==p.size()) {
                                                    success=true;
                                                    for (auto g=spc.groups.end()-N; g!=spc.groups.end(); ++g)
                                                        g->cm = Geometry::massCenter(g->begin(), g->end(), spc.geo.boundaryFunc, -g->begin()->pos);
                                                }
                                            } else std::cerr << file + ": wrong number of atoms" << endl;
                                        } else std::cerr << "error opening file '" + file + "'" << endl;
                                        if (success==false)
                                            throw std::runtime_error("error loading positions from '" + file + "'");
                                    }
                                }
                            } else
                                throw std::runtime_error("cannot insert undefined molecule '" + it.key() + "'");
                        }
                }
            } else throw std::runtime_error("'insertmolecules' json entry must be of array type");
        } //!< Insert `N` molecules into space as defined in `insert`

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Space")
    {
        typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
        typedef Space<Geometry::Cuboid, Tparticle> Tspace;
        Tspace spc1;

        // check molecule insertion
        atoms.resize(2);
        CHECK( atoms.at(0).mw == 1);
        Tparticle a;
        a.id=0;
        Tspace::Tpvec p(2, a);
        CHECK( p[0].traits().mw == 1);
        p[0].pos.x()=2;
        p[1].pos.x()=3;
        spc1.push_back(1, p);
        CHECK( spc1.p.size()==2 );
        CHECK( spc1.groups.size()==1 );
        CHECK( spc1.groups.front().id==1);
        CHECK( spc1.groups.front().cm.x()==doctest::Approx(2.5));

        // check `positions()`
        CHECK( &spc1.positions()[0] == &spc1.p[0].pos );

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
