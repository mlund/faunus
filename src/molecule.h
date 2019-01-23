#pragma once

#include <set>
#include "core.h"
#include "io.h"
#include "geometry.h"
#include "potentials.h"

namespace Faunus {

    namespace Potential { struct BondData; }

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
            Point dir={1,1,1};      //!< Scalars for random mass center position. Default (1,1,1)
            Point offset={0,0,0};   //!< Added to random position. Default (0,0,0)
            bool rotate=true;       //!< Set to true to randomly rotate molecule when inserted. Default: true
            bool keeppos=false;     //!< Set to true to keep original positions (default: false)
            int maxtrials=2e4;      //!< Maximum number of container overlap checks
            int confindex=-1;       //!< Index of last used conformation

            Tpvec operator()( Geometry::GeometryBase &geo, const Tpvec&, TMoleculeData &mol )
            {
                int cnt = 0;
                QuaternionRotate rot;
                bool containerOverlap;

                if (std::fabs(geo.getVolume())<1e-20)
                    throw std::runtime_error("geometry has zero volume");

                Tpvec v = mol.conformations.get(); // get random, weighted conformation
                confindex = mol.conformations.index; // lastest index

                do {
                    if ( cnt++ > maxtrials )
                        throw std::runtime_error("Max. # of overlap checks reached upon insertion.");

                    if ( mol.atomic ) { // insert atomic species
                        for ( auto &i : v ) { // for each atom type id
                            if ( rotate ) {
                                rot.set(2*pc::pi*random(), ranunit(random));
                                i.rotate(rot.first, rot.second);
                            }
                            geo.randompos(i.pos, random);
                            i.pos = i.pos.cwiseProduct(dir) + offset;
                            geo.boundary(i.pos);
                        }
                    }
                    else { // insert molecule
                        if ( keeppos ) {                     // keep original positions (no rotation/trans)
                            for ( auto &i : v )              // ...but let's make sure it fits
                                if (geo.collision(i.pos))
                                    throw std::runtime_error("Error: Inserted molecule does not fit in container");
                        } else {
                            Point cm;                       // new mass center position
                            geo.randompos(cm, random);      // random point in container
                            cm = cm.cwiseProduct(dir);      // apply user defined directions (default: 1,1,1)
                            Geometry::cm2origo(v.begin(), v.end());// translate to origin
                            rot.set(random()*2*pc::pi, ranunit(random)); // random rot around random vector
                            if (rotate) {
                                Geometry::rotate(v.begin(), v.end(), rot.first);
                                assert( Geometry::massCenter(v.begin(), v.end()).norm()<1e-6); // cm shouldn't move
                            }
                            for (auto &i : v) {
                                i.pos += cm + offset;
                                geo.boundary(i.pos);
                            }
                        }
                    }

                    if (v.empty())
                        throw std::runtime_error("Nothing to load/insert for molecule '"s + mol.name + "'");

                    // check if molecules / atoms fit inside simulation container
                    containerOverlap = false;
                    for ( auto &i : v )
                        if ( geo.collision(i.pos)) {
                            containerOverlap = true;
                            break;
                        }
                } while (containerOverlap);
                return v;
            }
        };

    /**
     * Possible structure for molecular conformations
     */
    struct Conformation {
        std::vector<Point> positions;
        std::vector<double> charges;

        bool empty() const {
            if (positions.empty())
                if (charges.empty())
                    return true;
            return false;
        }

        template<typename Tpvec>
            Tpvec& toParticleVector(Tpvec &p) const {
                assert(not p.empty() and not empty());
                // copy positions
                if (positions.size()==p.size())
                    for (size_t i=0; i<p.size(); i++)
                        p[i].pos = positions[i];
                // copy charges
                if (charges.size()==p.size())
                    for (size_t i=0; i<p.size(); i++)
                        p[i].charge = charges[i];
                return p;
            } // copy conformation into particle vector
    };
#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Conformation") {
        typedef std::vector<Particle<Charge>> Tpvec;
        Tpvec p(1);
        Conformation c;
        CHECK(c.empty());

        c.positions.push_back( {1,2,3} );
        c.charges.push_back( 0.5 );
        CHECK(not c.empty());

        c.toParticleVector(p);

        CHECK( p[0].pos == Point(1,2,3) );
        CHECK( p[0].charge == 0.5 );
    }
#endif
 
    /**
     * @brief General properties for molecules
     */
    template<class Tpvec>
        class MoleculeData {
            private:
                int _id=-1;
            public:
                typedef Tpvec TParticleVector;
                typedef typename Tpvec::value_type Tparticle;

                typedef std::function<Tpvec( Geometry::GeometryBase&,
                        const Tpvec&, MoleculeData<Tpvec>& )> TinserterFunc;

                TinserterFunc inserterFunctor=nullptr;      //!< Function for insertion into space

                int& id() { return _id; }  //!< Type id
                const int& id() const { return _id; } //!< Type id

                std::string name;          //!< Molecule name
                std::string structure;     //!< Structure file (pqr|aam|xyz)
                bool atomic=false;         //!< True if atomic group (salt etc.)
                bool rotate=true;          //!< True if molecule should be rotated upon insertion
                bool keeppos=false;        //!< Keep original positions of `structure`
                bool rigid=false;          //!< True if particle should be considered as rigid
                double activity=0;         //!< Chemical activity (mol/l)
                Point insdir = {1,1,1};    //!< Insertion directions
                Point insoffset = {0,0,0}; //!< Insertion offset

                std::vector<std::shared_ptr<Potential::BondData>> bonds;
                std::vector<int> atoms;    //!< Sequence of atoms in molecule (atom id's)
                WeightedDistribution<Tpvec> conformations;//!< Conformations of molecule

                MoleculeData() {
                    setInserter( RandomInserter<MoleculeData<Tpvec>>() );
                }

                /** @brief Specify function to be used when inserting into space.
                 *
                 * By default a random position and orientation is generator and overlap
                 * with container is avoided.
                 */
                void setInserter( const TinserterFunc &ifunc ) { inserterFunctor = ifunc; };

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

                void loadConformation(const std::string &file)
                {
                    Tpvec v;
                    if (loadStructure<Tpvec>()(file, v, false))
                    {
                        if ( keeppos == false )
                            Geometry::cm2origo( v.begin(), v.end() ); // move to origo
                        conformations.push_back(v);
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
                {"keeppos", a.keeppos}, {"bondlist", a.bonds},
                {"rigid", a.rigid}
            };
            if (not a.structure.empty())
                j[a.name]["structure"] = a.structure;

            j[a.name]["atoms"] = json::array();
            for (auto id : a.atoms)
                j[a.name]["atoms"].push_back( atoms.at(id).name );
        }

    template<class Tparticle, class Talloc>
        void from_json(const json& j, MoleculeData<std::vector<Tparticle,Talloc>> &a) {
            typedef typename std::vector<Tparticle,Talloc> Tpvec;

            try {
                if (j.is_object()==false || j.size()!=1)
                    throw std::runtime_error("invalid json");
                for (auto it : j.items()) {
                    a.name = it.key();
                    xjson val = it.value(); // keys are deleted after access
                    a.insoffset = val.value("insoffset", a.insoffset);
                    a.activity = val.value("activity", a.activity) * 1.0_molar;
                    a.keeppos = val.value("keeppos", a.keeppos);
                    a.atomic = val.value("atomic", a.atomic);
                    a.insdir = val.value("insdir", a.insdir);
                    a.bonds  = val.value("bondlist", a.bonds);
                    a.rigid = val.value("rigid", a.rigid);
                    a.id() = val.value("id", a.id());

                    if (a.atomic) {
                        // read `atoms` list of atom names and convert to atom id's
                        for (auto &i : val.at("atoms").get<std::vector<std::string>>()) {
                            auto it = findName( atoms, i );
                            if (it == atoms.end() )
                                throw std::runtime_error("unknown atoms in 'atoms'\n");
                            a.atoms.push_back(it->id());
                        }
                        assert(!a.atoms.empty());
                        assert(a.bonds.empty() && "bonds undefined for atomic groups");

                        // generate config
                        Tpvec v;
                        v.reserve( a.atoms.size() );
                        for ( auto id : a.atoms ) {
                            Tparticle _p;
                            _p = atoms.at(id);
                            v.push_back( _p );
                        }
                        if (!v.empty())
                            a.conformations.push_back(v);
                    }  // done handling atomic groups
                    else { // molecular groups
                        if (val.count("structure")>0) {
                            json _struct = val["structure"s];

                            // `structure` is a file name
                            if (_struct.is_string()) // structure from file
                                a.loadConformation( _struct.get<std::string>() );

                            else if (_struct.is_object()) {
                                // `structure` is a fasta sequence
                                if (_struct.count("fasta")) {
                                    Potential::HarmonicBond bond; // harmonic bond
                                    bond.from_json(_struct);      // read 'k' and 'req' from json
                                    std::string fasta = _struct.at("fasta").get<std::string>();
                                    Tpvec v = fastaToParticles<Tpvec>(fasta, bond.req);
                                    if (not v.empty()) {
                                        a.conformations.push_back(v);
                                        for (auto &p : v)
                                            a.atoms.push_back(p.id);
                                        // connect all atoms with harmonic bonds
                                        for (int i=0; i<(int)v.size()-1; i++) {
                                            bond.index = {i,i+1};
                                            a.bonds.push_back(bond.clone());
                                        }
                                    }
                                } // end of fasta handling
                            }

                            // `structure` is a list of atom positions
                            else if (_struct.is_array()) { // structure is defined inside json
                                Tpvec v;
                                a.atoms.clear();
                                v.reserve( _struct.size() );
                                for (auto &m : _struct)
                                    if (m.is_object())
                                        if (m.size()==1)
                                            for (auto& i : m.items()) {
                                                auto it = findName( atoms, i.key() );
                                                if (it == atoms.end())
                                                    throw std::runtime_error("unknown atoms in 'structure'");
                                                v.push_back( *it );     // set properties from atomlist
                                                v.back().pos = i.value(); // set position
                                                a.atoms.push_back(it->id());
                                            }
                                if (v.empty())
                                    throw std::runtime_error("invalid 'structure' format");
                                a.conformations.push_back(v);
                            } // end of position parser
                        } // end of `structure`

                        // read tracjectory w. conformations from disk
                        std::string traj = val.value("traj", std::string() );
                        if ( not traj.empty() ) {
                            a.conformations.clear();
                            FormatPQR::load(traj, a.conformations.vec);
                            if (not a.conformations.empty()) {
                                // create atom list
                                a.atoms.clear();
                                a.atoms.reserve(a.conformations.vec.front().size());
                                for ( auto &p : a.conformations.vec.front())           // add atoms to atomlist
                                    a.atoms.push_back(p.id);

                                // center mass center for each frame to origo assuming whole molecules
                                if (val.value("trajcenter", false)) {
                                    cout << "Centering conformations in trajectory file " + traj + ". ";
                                    for ( auto &p : a.conformations.vec ) // loop over conformations
                                        Geometry::cm2origo( p.begin(), p.end() );
                                    cout << "Done.\n";
                                }

                                // set default uniform weight
                                std::vector<float> w(a.conformations.size(), 1);
                                a.conformations.setWeight(w);

                                // look for weight file
                                std::string weightfile = val.value("trajweight", std::string());
                                if (not weightfile.empty()) {
                                    std::ifstream f( weightfile.c_str() );
                                    if (f) {
                                        w.clear();
                                        w.reserve(a.conformations.size());
                                        double _val;
                                        while ( f >> _val )
                                            w.push_back(_val);
                                        if ( w.size() == a.conformations.size())
                                            a.conformations.setWeight(w);
                                        else
                                            throw std::runtime_error("Number of weights does not match conformations.");
                                    } else
                                        throw std::runtime_error("Weight file " + weightfile + " not found.");
                                }
                            }
                            else
                                throw std::runtime_error("Trajectory " + traj + " not loaded or empty.");
                        } // done handling conformations

                    } // done handling molecular groups

                    // pass information to inserter
                    auto ins = RandomInserter<MoleculeData<std::vector<Tparticle,Talloc>>>();
                    ins.dir = a.insdir;
                    ins.offset = a.insoffset;
                    ins.keeppos = a.keeppos;
                    a.setInserter(ins);

                    // assert that all bonds are *internal*
                    for (auto &bond : a.bonds)
                        for (int i : bond->index)
                            if (i>=a.atoms.size() || i<0)
                                throw std::runtime_error("bonded atom index " + std::to_string(i) + " out of range");
                    // at this stage all given keys should have been accessed. If any are
                    // left, an exception will be thrown.
                    if (not val.empty())
                        throw std::runtime_error("unused key(s):\n"s + val.dump() + usageTip["moleculelist"]);
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

    /*
     * @brief General properties of reactions
     */
    template<class Tpvec>
        class ReactionData {
        private:
            int _id=-1;
        public:
            std::vector<std::string> _reac, _prod;

            typedef Tpvec TParticleVector;
            typedef typename Tpvec::value_type Tparticle;
            typedef std::map<int,int> Tmap;

            Tmap _reacid_m;     // Molecular change, groups. Atomic as Groupwise
            Tmap _reacid_a;     // Atomic change, equivalent of swap/titration
            Tmap _prodid_m;
            Tmap _prodid_a;
            //Tmap _Reac, _Prod;

            bool canonic=false;             //!< Finite reservoir
            int N_reservoir;                //!< Number of molecules in finite reservoir
            double lnK=0;                   //!< Natural logarithm of molar eq. const.
            double pK=0;                    //!< -log10 of molar eq. const.
            std::string name;               //!< Name of reaction
            std::string formula;            //!< Chemical formula
            double weight;                  //!< Statistical weight to be given to reaction in speciation

            bool empty(bool forward) const {
                if (forward)
                    if (canonic)
                        if (N_reservoir <= 0)
                            return true;
                return false;
            }

            std::vector<int> participatingMolecules() const {
                std::vector<int> v;
                v.reserve(_reacid_m.size()+_prodid_m.size());
                for (auto i : _reacid_m)
                    v.push_back(i.first);
                for (auto i : _prodid_m)
                    v.push_back(i.first);
                return v;
             } //!< Returns molids of participating molecules

            bool containsMolecule(int molid) const {
                if (_reacid_m.count(molid)==0)
                    if (_prodid_m.count(molid)==0)
                        return false;
                return true;
            } //!< True of molecule id is part of process

            const Tmap& Molecules2Add(bool forward) const {
                return (forward) ? _prodid_m : _reacid_m;
            } //!< Map for addition depending on direction

            const Tmap& Atoms2Add(bool forward) const {
                return (forward) ? _prodid_a : _reacid_a;
            } //!< Map for addition depending on direction

            auto findAtomOrMolecule(const std::string &name) const {
                auto it_a = findName(atoms, name);
                auto it_m = findName(molecules<Tpvec>, name);
                if (it_m == molecules<Tpvec>.end())
                    if (it_a == atoms.end())
                        throw std::runtime_error("unknown species '" + name + "'");
                return std::make_pair(it_a, it_m);
            } //!< Returns pair of iterators to atomlist and moleculelist. One of them points to end().

        }; //!< End of class

        inline auto parseProcess(const std::string &process) {
            typedef std::vector<std::string> Tvec;
            Tvec v;
            std::string tmp;
            std::istringstream iss(process);
            while (iss>>tmp)
                v.push_back(tmp);
            v.erase(std::remove(v.begin(), v.end(), "+"), v.end());
            auto it = std::find(v.begin(), v.end(), "=");
            if (it==v.end())
                throw std::runtime_error("products and reactants must be separated by '='");
            return std::make_pair( Tvec(v.begin(),it), Tvec(it+1, v.end()) );
        } //!< Parse process string to pair of vectors containing reactant/product species

        template<class Tparticle, class Talloc>
            void from_json(const json& j, ReactionData<std::vector<Tparticle,Talloc>> &a) {

                typedef std::vector<Tparticle,Talloc> Tpvec; // alias for particle vector

                if (j.is_object()==false || j.size()!=1)
                    throw std::runtime_error("Invalid JSON data for ReactionData");

                for (auto &m: molecules<Tpvec>)
                    for (auto &a: atoms)
                        if (m.name == a.name)
                            throw std::runtime_error("Molecules and atoms nust have different names");

                for (auto it=j.begin(); it!=j.end(); ++it) {
                    a.name = it.key();
                    auto& val = it.value();
                    a.canonic = val.value("canonic", false);
                    if (val.count("lnK")==1)
                        a.lnK = val.at("lnK").get<double>();
                    else if (val.count("pK")==1)
                        a.lnK = -std::log(10) * val.at("pK").get<double>();
                    a.pK = - a.lnK / std::log(10);
                    a.N_reservoir = val.value("N_reservoir", a.N_reservoir);

                    // get pair of vectors containing reactant and product species
                    auto process = parseProcess(a.name);
                    a._reac = process.first;
                    a._prod = process.second;

                    for (auto &name : a._reac) { // loop over reactants
                        auto pair = a.findAtomOrMolecule( name );  // {iterator to atom, iterator to mol.}
                        if ( pair.first != atoms.end() ) {
                            if ( pair.first->implicit ) {
                                a.lnK -= std::log(pair.first->activity/1.0_molar);
                            } else {
                                a._reacid_a[ pair.first->id() ]++;
                            }
                        }
                        if ( pair.second != molecules<Tpvec>.end() ) {
                            a._reacid_m[ pair.second->id() ]++;
                            if ( pair.second->activity > 0 ) {
                                a.lnK -= std::log(pair.second->activity/1.0_molar);
                            }
                        }
                    }

                    for (auto &name : a._prod) { // loop over products
                        auto pair = a.findAtomOrMolecule( name );
                        if ( pair.first != atoms.end() ) {
                            if ( pair.first->implicit ) {
                                a.lnK += std::log(pair.first->activity/1.0_molar);
                            } else {
                                a._prodid_a[ pair.first->id() ]++;
                            }
                        }
                        if ( pair.second != molecules<Tpvec>.end() ) {
                            a._prodid_m[ pair.second->id() ]++;
                            if ( pair.second->activity > 0 ) {
                                a.lnK += std::log(pair.second->activity/1.0_molar);
                            }
                        }
                    }
                }
            }

        template<class Tparticle, class Talloc>
            void to_json(json& j, const ReactionData<std::vector<Tparticle,Talloc>> &a) {
                j[a.name] = {
                    {"original pK", a.pK }, {"activity-modified pK", -a.lnK/std::log(10) },
                    {"canonic", a.canonic }, {"N_reservoir", a.N_reservoir },
                    {"products", a._prod },
                    //{"exchange products", a._prodid_m  },
                    {"reactants", a._reac }
                    //{"exchange reactants", a._reacid_m  }
                };
            } //!< Serialize to JSON object

        template<typename Tpvec>
            static std::vector<ReactionData<Tpvec>> reactions = {}; //!< Global instance of reaction list

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] ReactionData") {
        using doctest::Approx;

        json j = R"(
            {
                "atomlist" :
                    [ {"a": { "r":1.1 } } ],
                "moleculelist": [
                    { "B": { "activity":0.2, "atomic":true, "atoms":["a"] } },
                    { "A": { "atomic":false } }
                ],
                "reactionlist": [
                    {"A = B": {"lnK": -10.051, "canonic": true, "N": 100 } }
                ]
            } )"_json;

        typedef std::vector<Particle<Radius, Charge, Dipole, Cigar>> Tpvec;

        atoms = j["atomlist"].get<decltype(atoms)>();
        molecules<Tpvec> = j["moleculelist"].get<decltype(molecules<Tpvec>)>(); // fill global instance

        auto &r = reactions<Tpvec>; // reference to global reaction list
        r = j["reactionlist"].get<decltype(reactions<Tpvec>)>();

        CHECK( r.size()==1 );
        CHECK( r.front().name=="A = B" );
        CHECK( r.front().lnK==-10.051);
    }
#endif

}//namespace
