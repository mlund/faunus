#pragma once

#include "core.h"
#include "io.h"
#include "geometry.h"
#include "potentials.h"

namespace Faunus {

    /**
     * @brief Random position and orientation - typical for rigid bodies
     *
     * Molecule inserters take care of generating molecules
     * for insertion into space and can be used in Grand Canonical moves,
     * Widom analysis, and for generating initial configurations.
     * Inserters will not actually insert anything, but rather
     * return a particle vector with proposed coordinates.
     *
     * All inserters are function objects, expecting
     * a geometry, particle vector, and molecule data.
     */
    template<typename TMoleculeData>
        struct RandomInserter
        {
            typedef typename TMoleculeData::TParticleVector Tpvec;
            std::string name;
            Point dir={1,1,1};      //!< Scalars for random mass center position. Default (1,1,1)
            Point offset={0,0,0};   //!< Added to random position. Default (0,0,0)
            bool checkOverlap=true; //!< Set to true to enable container overlap check
            bool rotate=true;       //!< Set to true to randomly rotate molecule when inserted. Default: true
            bool keeppos=false;     //!< Set to true to keep original positions (default: false)
            int maxtrials=2e3;      //!< Maximum number of overlap checks if `checkOverlap==true`

            RandomInserter() : name("random") {}

            Tpvec operator()( Geometry::GeometryBase &geo, const Tpvec &p, TMoleculeData &mol )
            {
                if (std::fabs(geo.getVolume())<1e-20)
                    throw std::runtime_error("geometry has zero volume");
                bool _overlap = true;
                Tpvec v;
                int cnt = 0;
                do
                {
                    if ( cnt++ > maxtrials )
                        throw std::runtime_error("Max. # of overlap checks reached upon insertion.");

                    v = mol.getRandomConformation();

                    if ( mol.atomic )
                    { // insert atomic species
                        for ( auto &i : v )
                        { // for each atom type id
                            QuaternionRotate rot;
                            if ( rotate )
                            {
                                rot.set(2*pc::pi*random(), ranunit(random));
                                i.rotate(rot.first, rot.second);
                            }
                            geo.randompos(i.pos, random);
                            i.pos = i.pos.cwiseProduct(dir) + offset;
                            geo.boundary(i.pos);
                        }
                    }
                    else
                    { // insert molecule
                        if ( keeppos )
                        {                   // keep original positions (no rotation/trans)
                            for ( auto &i : v )              // ...but let's make sure it fits
                                if ( geo.collision(i.pos, 0))
                                    throw std::runtime_error("Error: Inserted molecule does not fit in container");
                        }
                        else
                        {
                            Point cm;                      // new mass center position
                            geo.randompos(cm, random);      // random point in container
                            cm = cm.cwiseProduct(dir);     // apply user defined directions (default: 1,1,1)
                            Geometry::cm2origo(v.begin(), v.end());// translate to origin
                            QuaternionRotate rot;
                            rot.set(random()*2*pc::pi, ranunit(random)); // random rot around random vector
                            for ( auto &i : v )
                            {            // apply rotation to all points
                                if ( rotate ) {
                                    i.rotate(rot.first, rot.second);    // internal atom rotation
                                    i.pos = rot(i.pos) + cm + offset;   // ...and translate
                                }
                                else
                                    i.pos += cm + offset;
                                geo.boundary(i.pos);             // ...and obey boundaries
                            }
                        }
                    }

                    assert(!v.empty());

                    _overlap = false;
                    if ( checkOverlap )              // check for container overlap
                        for ( auto &i : v )
                            if ( geo.collision(i.pos))
                            {
                                _overlap = true;
                                break;
                            }
                }
                while ( _overlap == true );
                return v;
            }
        };

    /**
     * @brief General properties for molecules
     */
    template<class Tpvec>
        class MoleculeData {
            private:
                int _id=-1;
                int _confid=-1;
            public:
                typedef Tpvec TParticleVector;
                typedef typename Tpvec::value_type Tparticle;

                /** @brief Signature for inserted function */
                typedef std::function<Tpvec( Geometry::GeometryBase &, const Tpvec &, MoleculeData<Tpvec> & )> TinserterFunc;
                TinserterFunc inserterFunctor=nullptr;      //!< Function for insertion into space

                int& id() { return _id; } //!< Type id
                const int& id() const { return _id; } //!< Type id

                std::string name;          //!< Molecule name
                std::string structure;     //!< Structure file (pqr|aam|xyz)
                bool atomic=false;         //!< True if atomic group (salt etc.)
                bool rotate=true;          //!< True if molecule should be rotated upon insertion
                bool keeppos=false;        //!< Keep original positions of `structure`
                double activity=0;         //!< Chemical activity (mol/l)
                Point insdir = {1,1,1};    //!< Insertion directions
                Point insoffset = {0,0,0}; //!< Insertion offset

                std::vector<Potential::BondData> bonds;
                std::vector<int> atoms;    //!< Sequence of atoms in molecule (atom id's)
                std::vector<Tpvec> conformations;           //!< Conformations of molecule
                std::discrete_distribution<> confDist;      //!< Weight of conformations

                MoleculeData() {
                    setInserter( RandomInserter<MoleculeData<Tpvec>>() );
                }

                void addConformation( const Tpvec &vec, double weight = 1 )
                {
                    if ( !conformations.empty())
                    {     // resize weights
                        auto w = confDist.probabilities();// (defaults to 1)
                        w.push_back(weight);
                        confDist = std::discrete_distribution<>(w.begin(), w.end());
                    }
                    conformations.push_back(vec);
                    assert(confDist.probabilities().size() == conformations.size());
                } //!< Store a single conformation

                /** @brief Specify function to be used when inserting into space.
                 *
                 * By default a random position and orientation is generator and overlap
                 * with container is avoided.
                 */
                void setInserter( const TinserterFunc &ifunc ) { inserterFunctor = ifunc; };

                /**
                 * @brief Get a random conformation
                 *
                 * This will return the raw coordinates of a random conformation
                 * as loaded from a directory file. The propability of a certain
                 * conformation is dictated by the weight which, by default,
                 * is set to unity. Specify a custom distribution using the
                 * `weight` keyword.
                 */
                Tpvec getRandomConformation()
                {
                    if ( conformations.empty())
                        throw std::runtime_error("No configurations for molecule '" + name +
                                "'. Perhaps you forgot to specity the 'atomic' keyword?");

                    assert(size_t(confDist.max()) == conformations.size() - 1);
                    //assert(atoms.size() == conformations.front().size());

                    _confid = confDist(random.engine); // store the index of the conformation
                    return conformations.at( _confid );
                }

                /**
                 * @brief Get random conformation that fits in container
                 * @param geo Geometry
                 * @param otherparticles Typically `spc.p` is insertion depends on other particle
                 *
                 * By default the molecule is placed at a random position and orientation with
                 * no container overlap using the `RandomInserter` class. This behavior can
                 * be changed by specifying another inserter using `setInserter()`.
                 */
                Tpvec getRandomConformation(Geometry::GeometryBase &geo, Tpvec otherparticles = Tpvec())
                {
                    assert(inserterFunctor!=nullptr);
                    return inserterFunctor(geo, otherparticles, *this);
                }

                /**
                 * @brief Store a single conformation
                 * @param vec Vector of particles
                 * @param weight Relative weight of conformation (default: 1)
                 */
                void pushConformation( const Tpvec &vec, double weight = 1 )
                {
                    if (!vec.empty()) {
                        if ( !conformations.empty())
                        {     // resize weights
                            auto w = confDist.probabilities();// (defaults to 1)
                            w.push_back(weight);
                            confDist = std::discrete_distribution<>(w.begin(), w.end());
                        }
                        conformations.push_back(vec);
                        assert(confDist.probabilities().size() == conformations.size());
                    } else
                        std::cerr << "MoleculeData: attempt to insert empty configuration\n";
                }

                /** @brief Nunber of conformations stored for molecule */
                size_t numConformations() const
                {
                    return conformations.size();
                }

                void loadConformation(const std::string &file)
                {
                    Tpvec v;
                    if (loadStructure<Tpvec>()(file, v, false))
                    {
                        if ( keeppos == false )
                            Geometry::cm2origo( v.begin(), v.end() ); // move to origo
                        pushConformation(v);        // add conformation
                        for ( auto &p : v )           // add atoms to atomlist
                            atoms.push_back(p.id);
                    }
                    if ( v.empty() )
                        throw std::runtime_error("Structure " + structure + " not loaded. Filetype must be .aam/.pqr/.xyz");
                }
        }; // end of class

    template<class Tparticle, class Talloc>
        void to_json(json& j, const MoleculeData<std::vector<Tparticle,Talloc>> &a) {
            j[a.name] = {
                {"activity", a.activity/1.0_molar}, {"atomic", a.atomic},
                {"id", a.id()}, {"insdir", a.insdir}, {"insoffset", a.insoffset},
                {"keeppos", a.keeppos}, {"structure", a.structure}, {"bondlist", a.bonds}
            };
            j[a.name]["atoms"] = json::array();
            for (auto id : a.atoms)
                j[a.name]["atoms"].push_back( atoms<Tparticle>.at(id).name );
        }

    template<class Tparticle, class Talloc>
        void from_json(const json& j, MoleculeData<std::vector<Tparticle,Talloc>> &a) {
            typedef typename std::vector<Tparticle,Talloc> Tpvec;

            try {
                if (j.is_object()==false || j.size()!=1)
                    throw std::runtime_error("invalid json");
                for (auto it=j.begin(); it!=j.end(); ++it) {
                    a.name = it.key();
                    auto& val = it.value();
                    a.insoffset = val.value("insoffset", a.insoffset);
                    a.activity = val.value("activity", a.activity) * 1.0_molar;
                    a.keeppos = val.value("keeppos", a.keeppos);
                    a.atomic = val.value("atomic", a.atomic);
                    a.insdir = val.value("insdir", a.insdir);
                    a.bonds = val.value("bondlist", a.bonds);
                    a.id() = val.value("id", a.id());

                    if (a.atomic) {
                        // read `atoms` list of atom names and convert to atom id's
                        for (auto &i : val.at("atoms").get<std::vector<std::string>>()) {
                            auto it = findName( atoms<Tparticle>, i );
                            if (it == atoms<Tparticle>.end() )
                                throw std::runtime_error("unknown atoms in 'atoms'\n");
                            a.atoms.push_back(it->id());
                        }
                        assert(!a.atoms.empty());

                        // generate config
                        Tpvec v;
                        v.reserve( a.atoms.size() );
                        for ( auto id : a.atoms )
                            v.push_back( atoms<Tparticle>.at(id).p );
                        if (!v.empty())
                            a.pushConformation( v );
                    } else 
                        if (val.count("structure")>0) {
                            json _struct = val["structure"];
                            if (_struct.is_string()) // structure from file
                                a.loadConformation( val.value("structure", a.structure) );
                            else
                                if (_struct.is_array()) { // structure is defined inside json
                                    Tpvec v; // temporary particle vector
                                    v.reserve( _struct.size() );
                                    for (auto &m : _struct)
                                        for (auto i=m.begin(); i!=m.end(); ++i)
                                            if (i->is_object() && i->size()==1) {
                                                auto it = findName( atoms<Tparticle>, i.key() );
                                                if (it == atoms<Tparticle>.end() )
                                                    throw std::runtime_error("unknown atoms in 'structure'");
                                                v.push_back( it->p );     // set properties from atomlist
                                                v.back().pos = i.value(); // set position
                                            }
                                    if (v.empty())
                                        throw std::runtime_error("invalid 'structure' format");
                                    a.pushConformation( v );
                                }
                        }

                    // pass information to inserter
                    auto ins = RandomInserter<MoleculeData<std::vector<Tparticle,Talloc>>>();
                    ins.dir = a.insdir;
                    ins.offset = a.insoffset;
                    ins.keeppos = a.keeppos;
                    a.setInserter(ins);
                }
            } catch(std::exception& e) {
                throw std::runtime_error("JSON->molecule: " + a.name + ": " + e.what());
            }
        }

    template<class Tparticle, class Talloc>
        void from_json(const json& j, std::vector<MoleculeData<std::vector<Tparticle,Talloc>>> &v) {
            v.reserve( v.size() + j.size());
            for (auto &i : j) {
                v.push_back(i);
                v.back().id() = v.size()-1; // id always match vector index
            }
        }

    template<typename Tpvec>
        static std::vector<MoleculeData<Tpvec>> molecules = {}; //!< Global instance of molecule list

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] MoleculeData") {
        using doctest::Approx;

        json j = R"(
            { "moleculelist": [
                { "B": {"activity":0.2, "atomic":true, "insdir": [0.5,0,0], "insoffset": [-1.1, 0.5, 10], "atoms":["A"] } },
                { "A": { "atomic":false } }
            ]})"_json;

        typedef Particle<Radius, Charge, Dipole, Cigar> T;
        typedef std::vector<T> Tpvec;
        typedef MoleculeData<Tpvec> Tmoldata;

        molecules<Tpvec> = j["moleculelist"].get<decltype(molecules<Tpvec>)>(); // fill global instance
        auto &v = molecules<Tpvec>; // reference to global molecule vector

        CHECK(v.size()==2);
        CHECK(v.back().id()==1);
        CHECK(v.back().name=="A"); // alphabetic order in std::map
        CHECK(v.back().atomic==false);

        MoleculeData<Tpvec> m = json(v.front()); // moldata --> json --> moldata

        CHECK(m.name=="B");
        CHECK(m.id()==0);
        CHECK(m.activity==Approx(0.2_molar));
        CHECK(m.atomic==true);
        CHECK(m.insdir==Point(0.5,0,0));
        CHECK(m.insoffset==Point(-1.1,0.5,10));
    }
#endif

}//namespace
