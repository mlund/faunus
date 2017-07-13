#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <string>
#include <json.hpp>
#include <cmath>
#include <random>
#include <Eigen/Geometry>
#include <range.hpp>

// Eigen<->JSON (de)serialization
namespace Eigen {
    template<typename T>
        void to_json(nlohmann::json& j, const T &p) {
            auto d = p.data();
            for (size_t i = 0; i<p.size(); ++i)
                j.push_back( d[i] ); 
        }
    void from_json(const nlohmann::json& j, Vector3d &p) {
        if ( j.size()==3 ) {
            int i=0;
            for (auto d : j.get<std::vector<double>>())
                p[i++] = d;
            return;
        }
        throw std::runtime_error("JSON->Eigen conversion error");
    }
};

/** @brief Faunus main namespace */
namespace Faunus {

    typedef Eigen::Vector3d Point; //!< 3d vector
    typedef nlohmann::json Tjson;  //!< Json object
    using std::cout;
    using std::endl;

    /** @brief Physical constants */
    namespace PhysicalConstants {
        typedef double T; //!< Float size
        constexpr T infty = std::numeric_limits<T>::infinity(), //!< Numerical infinity
                  pi = 3.141592653589793, //!< Pi
                  e0 = 8.85419e-12,  //!< Permittivity of vacuum [C^2/(J*m)]
                  e = 1.602177e-19,  //!< Absolute electronic unit charge [C] 
                  kB = 1.380658e-23, //!< Boltzmann's constant [J/K]
                  Nav = 6.022137e23, //!< Avogadro's number [1/mol]
                  c = 299792458.0,   //!< Speed of light [m/s]
                  R = kB * Nav;      //!< Molar gas constant [J/(K*mol)] 

        static T temperature=298.15; //!< Temperature (Kelvin)
        static inline T kT() { return temperature*kB; } //!< Thermal energy (Joule)
        static inline T lB( T epsilon_r ) {
            return e*e/(4*pi*e0*epsilon_r*1e-10*kT());
        } //!< Bjerrum length (angstrom)

    }
    namespace pc = PhysicalConstants;

    /**
     * @brief Chemistry units
     *
     * String literals to convert to the following internal units:
     *
     * Property           | Unit
     * :----------------- | :--------------------------
     * Energy             | Thermal energy (kT)
     * Temperature        | Kelvin (K)
     * Length             | Angstrom (A) 
     * Charge             | Electron unit charge (e)
     * Dipole moment      | Electron angstrom (eA)
     * Concentration      | Particles / angstrom^3
     * Pressure           | Particles / angstrom^3
     * Angle              | Radians
     */
    namespace ChemistryUnits
    {
        typedef long double T; //!< Floating point size
        constexpr T operator "" _K( T temp ) { return temp; } //!< Kelvin to Kelvin
        constexpr T operator "" _C( T temp ) { return 273.15 + temp; } //!< Celcius to Kelvin
        constexpr T operator "" _Debye( T mu ) { return mu * 0.208194334424626; }  //!< Debye to eA
        constexpr T operator "" _eA( T mu ) { return mu; } //!< eA to eA
        constexpr T operator "" _Cm( T mu ) { return mu * 1.0_Debye / 3.335640951981520e-30; } //!< Cm to eA
        constexpr T operator "" _angstrom( T l ) { return l; } //!< Angstrom to Angstrom
        constexpr T operator "" _m( T l ) { return l * 1e10; } //!< Meter to angstrom
        constexpr T operator "" _bohr( T l ) { return l * 0.52917721092; } //!< Bohr to angstrom
        constexpr T operator "" _nm( T l ) { return l * 10; } //!< nanometers to angstrom
        constexpr T operator "" _liter( T v ) { return v * 1e27; } //!< liter to angstrom cubed
        constexpr T operator "" _m3( T v ) { return v * 1e30; } //!< -> cubic meter to angstrom cubed
        constexpr T operator "" _mol( T n ) { return n * pc::Nav; } //!< moles to particles
        constexpr T operator "" _molar( T c ) { return c * 1.0_mol / 1.0_liter; } //!< molar to particle / angstrom^3
        constexpr T operator "" _mM( T c ) { return c * 1.0e-3_mol / 1.0_liter; } //!< millimolar to particle / angstrom^3
        constexpr T operator "" _rad( T a ) { return a; } //!< Radians to radians
        constexpr T operator "" _deg( T a ) { return a * 3.14159265358979323846 / 180; } //!< Degrees to radians
        inline T operator "" _Pa( T p ) { return p / pc::kT() / 1.0_m3; } //!< Pascal to particle / angstrom^3
        inline T operator "" _atm( T p ) { return p * 101325.0_Pa; } //!< Atmosphere to particle / angstrom^3
        inline T operator "" _bar( T p ) { return p * 100000.0_Pa; } //!< Bar to particle / angstrom^3
        constexpr T operator "" _kT( T u ) { return u; } //!< kT to kT (thermal energy)
        inline T operator "" _J( T u ) { return u / pc::kT(); } //!< Joule to kT
        inline T operator "" _kJmol( T u ) { return u / pc::kT() / pc::Nav * 1e3; } //!< kJ/mol to kT/particle
        inline T operator "" _kcalmol( T u ) { return u * 4.1868_kJmol; } //!< kcal/mol to kT/particle
        inline T operator "" _hartree( T u ) { return u * 4.35974434e-18_J; } //!< Hartree to kT
    }
    using namespace ChemistryUnits;

    /**
     * @brief Class for handling random number generation
     *
     * By default a deterministic seed is used, but a hardware
     * or saved state can be used via json (de)serialization:
     *
     * ```{.cpp}
     *     Random r1;                                           // deterministic seed
     *     Random r2 = Tjson(r1);                               // copy random engine state from r1
     *     Random r3 = R"( {"randomseed" : "hardware"} )"_json; // non-deterministic seed
     * ```
     */
    struct Random {
        std::mt19937 engine; //!< Random number engine used for all operations
        std::uniform_real_distribution<double> dist01; //!< Uniform real distribution [0,1)

        Random() : dist01(0,1) {}

        double operator()() { return dist01(engine); } //!< Double in uniform range [0,1)

        int range( int min, int max )
        {
            std::uniform_int_distribution<int> d(min, max);
            return d(engine);
        } //!< Integer in uniform range [min:max]

        template<class Titer>
            Titer element( const Titer &beg, const Titer &end )
            {
                auto i = beg;
                std::advance(i, range(0, std::distance(beg, end) - 1));
                return i;
            } //!< Iterator to random element in container (std::vector, std::map, etc)
    };

    void to_json(Tjson &j, const Random &r) {
        std::ostringstream o;
        o << r.engine;
        j["randomseed"] = o.str();
    } //!< Random to Tjson conversion

    void from_json(const Tjson &j, Random &r) {
        if (j.is_object()) {
            auto seed = j.value("randomseed", std::string());
            try {
                if (seed=="hardware")
                    r.engine = decltype(r.engine)(std::random_device()());
                else if (!seed.empty()) {
                    std::stringstream s(seed);
                    s.exceptions( std::ios::badbit | std::ios::failbit );
                    s >> r.engine;
                }
            }
            catch (std::exception &e) {
                std::cerr << "error initializing random from json: " << e.what();
                throw;
            }
        }
    } //!< Tjson to Random conversion

    /**
     * @brief Generate random unit vector
     *
     * Based on the von Neumann method described in
     * *Allen and Tildesley*, page 349, which ensures
     * a uniform distribution on a unit sphere. More information
     * about *Sphere Point Picking* can be found at
     * [MathWorld](http://mathworld.wolfram.com/SpherePointPicking.html).
     *
     * @param rand Function object that takes no arguments and returns a random
     *             float uniformly distributed in the range `[0:1[`.
     */
    Point ranunit( Random &rand )
    {
        Point p;
        double r2;
        do
        {
            for ( size_t i = 0; i < 3; ++i )
                p(i) = 2 * rand() - 1;
            r2 = p.squaredNorm();
        }
        while ( r2 > 1 );
        return p / std::sqrt(r2);
    }

    /** @brief Base class for particle properties */
    struct ParticlePropertyBase {
        virtual void to_json(Tjson &j) const=0; //!< Convert to JSON object
        virtual void from_json(const Tjson &j)=0; //!< Convert from JSON object
        inline void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d&) {}
    };

    template <typename... Ts>
        auto to_json(Tjson&) -> typename std::enable_if<sizeof...(Ts) == 0>::type {} //!< Particle to JSON
    template<typename T, typename... Ts>
        void to_json(Tjson& j, const ParticlePropertyBase &a, const Ts&... rest) {
            a.to_json(j);
            to_json<Ts...>(j, rest...);
        } //!< Particle to JSON

    // JSON --> Particle
    template <typename... Ts>
        auto from_json(const Tjson&) -> typename std::enable_if<sizeof...(Ts) == 0>::type {}
    template<typename T, typename... Ts>
        void from_json(const Tjson& j, ParticlePropertyBase &a, Ts&... rest) {
            a.from_json(j);
            from_json<Ts...>(j, rest...);
        }

    /** @brief Dipole properties
     *
     * Json (de)serialization:
     *
     * ```{.cpp}
     *     Dipole d = R"( "mu":[0,0,1], "mulen":10 )"_json
     * ```
     */
    struct Dipole : public ParticlePropertyBase {
        Point mu={1,0,0}; //!< dipole moment unit vector
        double mulen=0;   //!< dipole moment scalar

        void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d&) {
            cout << "dipole rotate\n";
            mu = q * mu;
        } //!< Rotate dipole moment

        void to_json(Tjson &j) const override {
            j["mulen"] = mulen ;
            j["mu"] = mu;
        }

        void from_json(const Tjson& j) override {
            mulen = j.at("mulen").get<double>();
            mu = j.value("mu", Point(1,0,0) );
        }
    };

    struct Radius : public ParticlePropertyBase {
        double radius; //!< Particle radius
        void to_json(Tjson &j) const override { j["r"] = radius; }
        void from_json(const Tjson& j) override { radius = j.value("r", 0); }
    };

    /** @brief Sphero-cylinder properties */
    struct Cigar : public ParticlePropertyBase {
        double sclen=0;       //!< Sphero-cylinder length
        Point scdir = {1,0,0};//!< Sphero-cylinder direction unit vector

        void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d&) {
            cout << "cigar rotate\n";
            scdir = q * scdir;
        } //!< Rotate sphero-cylinder

        void to_json(Tjson &j) const override {
            j["sclen"] = sclen;
            j["scdir"] = scdir;
        }
        void from_json(const Tjson& j) override {
            sclen = j.at("sclen").get<double>();
            scdir = j.value("scdir", Point(1,0,0) );
        }
    };

    /** @brief Particle
     *
     * In addition to basic particle properties (id, charge, weight), arbitrary
     * features can be added via template arguments:
     *
     * ```{.cpp}
     *     Particle<Dipole> p
     *     p = R"( {"mulen":2.8} )"_json; // Json --> Particle
     *     std::cout << Tjson(p); // Particle --> Json --> stream
     * ```
     *
     * Currently the following properties are implemented:
     *
     * Keyword   | Class       |  Description
     * --------- | ----------- | --------------------
     *  `id`     | `Particle`  |  Type id (int)
     *  `mw`     | `Particle`  |  molecular weight (g/mol)
     *  `q`      | `Particle`  |  valency (e)
     *  `r`      | `Radius`  |  radius (angstrom)
     *  `mu`     | `Dipole`    |  dipole moment unit vector (array)
     *  `mulen`  | `Dipole`    |  dipole moment scalar (eA)
     *  `scdir`  | `Cigar`     |  Sphero-cylinder direction unit vector (array)
     *  `sclen`  | `Cigar`     |  Sphero-cylinder length (A)
     */
    template<typename... Properties>
        class Particle : public Properties... {
            private:
                // rotate internal coordinates
                template <typename... Ts>
                    auto _rotate(const Eigen::Quaterniond&, const Eigen::Matrix3d&) -> typename std::enable_if<sizeof...(Ts) == 0>::type {}
                template<typename T, typename... Ts>
                    void _rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m, T &a, Ts&... rest) {
                        a.rotate(q, m);
                        _rotate<Ts...>(q, m, rest...);
                    }

            public:
                Point pos;       //!< Particle position vector
                int id=-1;       //!< Particle type
                double _radius=0, //!< Particle radius
                       charge=0, //!< Particle charge
                       mw=0;     //!< Particle molecular weight

                Particle() : Properties()... {}

                void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) {
                    _rotate<Properties...>(q, m, dynamic_cast<Properties&>(*this)...);
                } //!< Rotate all internal coordinates if needed
        };

    template<typename... Properties>
        void to_json(Tjson& j, const Particle<Properties...> &a) {
            j["q"] = a.charge; j["mw"]=a.mw;
            //j["r"]=a.radius;
            j["id"]=a.id;
            j["pos"] = a.pos; 
            to_json<Properties...>(j, Properties(a)... );
        }

    template<typename... Properties>
        void from_json(const Tjson& j, Particle<Properties...> &a) {
            a.charge = j.value<double>("q", 0.0); 
            //a.radius = j.value<double>("r", 0.0); 
            a.mw = j.value<double>("mw", 0.0); 
            a.id = j.value<int>("id",-1);
            a.pos = j.value("pos", Point(0,0,0));
            from_json<Properties...>(j, dynamic_cast<Properties&>(a)...);
        }

    /** @brief Simulation geometries and related operations */
    namespace Geometry {

        struct GeometryBase {
            virtual void setVolume(double, const std::vector<double>&)=0; //!< Set volume
            virtual double getVolume(int=3) const=0; //!< Get volume
            virtual void randompos( Point&, std::function<double()>& ) const=0; //!< Generate random position
            virtual Point vdist( const Point &a, const Point &b ) const=0; //!< (Minimum) distance between two points
            virtual void boundary( Point &a ) const=0; //!< Apply boundary conditions

            std::function<void(Point&)> boundaryFunc; //!< Functor for boundary()

            GeometryBase() {
                boundaryFunc = std::bind( &GeometryBase::boundary, this, std::placeholders::_1);
            }

        }; //!< Base class for all geometries

        /** @brief Cuboidal box */
        class Box : public GeometryBase {
            protected:
                Point len, //!< side length
                      len_half, //!< half side length
                      len_inv; //!< inverse side length

            public:

                void setLength( const Point &l ) {
                    len = l;
                    len_half = l*0.5;
                    len_inv = l.cwiseInverse();
                } //!< Set cuboid side length

                void setVolume(double V, const std::vector<double> &s={1,1,1}) override {
                    double l = std::cbrt(V);
                    setLength( {l,l,l} );
                }

                double getVolume(int dim=3) const override {
                    return len.x() * len.y() * len.z();
                }

                void randompos( Point &m, std::function<double()> &rand ) const override {
                    m.x() = (rand()-0.5) * this->len.x();
                    m.y() = (rand()-0.5) * this->len.y();
                    m.z() = (rand()-0.5) * this->len.z();
                }
        };

        /** @brief Periodic boundary conditions */
        template<bool X=true, bool Y=true, bool Z=true>
            struct PBC : public Box {

                Point vdist( const Point &a, const Point &b ) const override
                {
                    Point r(a - b);
                    if (X) {
                        if ( r.x() > len_half.x())
                            r.x() -= len.x();
                        else if ( r.x() < -len_half.x())
                            r.x() += len.x();
                    }
                    if (Y) {
                        if ( r.y() > len_half.y())
                            r.y() -= len.y();
                        else if ( r.y() < -len_half.y())
                            r.y() += len.y();
                    }
                    if (Z) {
                        if ( r.z() > len_half.z())
                            r.z() -= len.z();
                        else if ( r.z() < -len_half.z())
                            r.z() += len.z();
                    }
                    return r;
                } //!< (Minimum) distance between two points

                template<typename T=double>
                    inline int anint( T x ) const
                    {
                        return int(x > 0.0 ? x + 0.5 : x - 0.5);
                    } //!< Round to int

                void boundary( Point &a ) const override
                {
                    if (X)
                        if ( std::fabs(a.x()) > len_half.x())
                            a.x() -= len.x() * anint(a.x() * len_inv.x());
                    if (Y)
                        if ( std::fabs(a.y()) > len_half.y())
                            a.y() -= len.y() * anint(a.y() * len_inv.y());
                    if (Z)
                        if ( std::fabs(a.z()) > len_half.z())
                            a.z() -= len.z() * anint(a.z() * len_inv.z());
                } //!< Apply boundary conditions

            };

        /** @brief Cylindrical cell */
        class Cylinder : public PBC<false,false,true> {
            private:
                double r, r2, diameter, len;
                typedef PBC<false,false,true> base;
            public:
                //std::function<void(Point&)> boundaryFunc;

                //Cylinder() {
                //    using namespace std::placeholders;
                //    boundaryFunc = std::bind( &Cylinder::boundary, this, _1);
                //}

                void setRadius(double radius, double length) {
                    len = length;
                    r = radius;
                    r2 = r*r;
                    diameter = 2*r;
                    Box::setLength( { diameter, diameter, len } ); 
                } //!< Set radius

                void setVolume(double V, const std::vector<double> &s={}) override {
                    r = std::sqrt( V / (pc::pi * len) );
                    setRadius( r2, len);
                }

                double getVolume(int dim=3) const override {
                    if (dim==1)
                        return len;
                    if (dim==2)
                        return pc::pi * r2;
                    return r2 * pc::pi * len;
                }

                void randompos( Point &m, std::function<double()>& rand ) const override
                {
                    double l = r2 + 1;
                    m.z() = (rand()-0.5) * len;
                    while ( l > r2 )
                    {
                        m.x() = (rand()-0.5) * diameter;
                        m.y() = (rand()-0.5) * diameter;
                        l = m.x() * m.x() + m.y() * m.y();
                    }
                }

        };

        /** @brief Spherical cell */
        class Sphere : public PBC<false,false,false> {
            private:
                double r;
            public:
        };

        /** @brief Translate a particle vector by a vector */
        template<class T>
            void translate( std::vector<T> &p, const Point &d, const std::function<void(Point&)>& boundary=[](Point&){} )
            {
                for ( auto &i : p )
                {
                    i.pos += d;
                    boundary(i.pos);
                }
            }

        enum class weight { MASS, CHARGE, GEOMETRIC };

        template<class Titer>
            Point anyCenter( const std::pair<Titer,Titer> &iter, weight k) {
                typedef typename Titer::value_type Tparticle;
                std::function<double(Tparticle&)> f;
            } //!< Mass, charge, or geometric center of a collection of particles

        using Cuboid = PBC<true, true, true>; //!< Cuboid w. PBC in all directions
        using Cuboidslit = PBC<true, true, false>; //!< Cuboidal slit w. PBC in XY directions

    } //geometry namespace

    namespace Analysis {

        /** @brief Example analysis */
        template<class T, class Enable = void>
            struct analyse {
                void sample(T &p) {
                    std::cout << "not a dipole!" << std::endl;
                } //!< Sample
            }; // primary template

        /** @brief Example analysis */
        template<class T>
            struct analyse<T, typename std::enable_if<std::is_base_of<Dipole, T>::value>::type> {
                void sample(T &p) {
                    std::cout << "dipole!" << std::endl;
                } //!< Sample
            }; // specialized template

    }//namespace

    struct Group : public Range::range {
        Point cm;     //!< Mass center
        int id;       //!< Group ID
        bool atomic;  //!< True if container species are atomic
        Group() {}
        Group(int beg, int end) : Range::range(beg, end+1, 1) {}

        int front() const {
            assert( !empty() );
            return *begin();
        }
        int back() const {
            assert( !empty() );
            return *(end()-1);
        }
    };

    /**
     * @brief General properties for atoms
     */
    struct AtomData {
        int id;            //!< Internal id
        double activity=0, //!< Chemical potential (mol/l)
               charge=0,   //!< Valency (e)
               mw=1,       //!< Molecular weight (g/mol)
               radius=0;   //!< Radius (A)
        std::string name;  //!< Name
    };

    void from_json(const Tjson& j, AtomData& a) {
        if ( j.is_array() )
            if ( j.size()==2 ) {
                a.name     = j[0].get<std::string>();
                a.activity = j[1].value("activity", 0.0) * 1.0_molar;
                a.charge   = j[1].value("q", 0.0);
                a.mw       = j[1].value("mw", 1.0);
                a.radius   = j[1].value("r", 0.0) * 1.0_angstrom;
                return;
            }
        throw std::runtime_error("Invalid JSON data for AtomData");
    }

    /**
     * @brief General properties for molecules
     */
    template<class Tpvec>
        struct MoleculeData {
            int id = -1,               //!< Internal molecule id
                Ninit = 0;             //!< Number of initial molecules
            std::string name;          //!< Molecule name
            bool atomic=false,         //!< True if atomic group (salt etc.)
                 rotate=true;          //!< True if molecule should be rotated upon insertion
            double activity=0;         //!< Chemical activity (mol/l)
            Point insdir = {1,1,1},    //!< Insertion directions
                  insoffset = {0,0,0}; //!< Insertion offset
            std::vector<int> atoms;    //!< Sequence of atoms in molecule (atom id's)
            std::vector<Tpvec> conformations;           //!< Conformations of molecule
            std::discrete_distribution<> confDist;      //!< Weight of conformations

            /**
             * @brief Store a single conformation
             * @param vec Vector of particles
             * @param weight Relative weight of conformation (default: 1)
             */
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
            }

        };

    template<class T>
        void from_json(const Tjson& j, MoleculeData<T>& m) {
            if ( j.is_array() )
                if ( j.size()==2 ) {
                    m.name   = j[0].get<std::string>();
                    m.Ninit  = j[1].value("Ninit", 0);
                    m.atomic = j[1].value("atomic", false);
                    m.activity = j[1].value("activity", 0.0) * 1.0_molar;
                }
        }

    template<class Tparticle>
        class Space {
            typedef std::vector<Tparticle> Tpvec;
            static std::vector<AtomData> atoms;
            static std::vector<MoleculeData<Tpvec>> molecules;
            Tpvec p;
            std::vector<Group> g;

            auto findMolecule(int molid) const {
                return std::find_if( g.begin(), g.end(),
                        [&molid](const Group &g) { return g.id==molid; } );
            }

            void insert(int molid, int N=1) {
            }

            template<class Tindex>
                void erase(int molid, Tindex ndx );

            auto findAtomData(const std::string &name) const {
                return std::find_if( atoms.begin(), atoms.end(),
                        [&name](const AtomData &a) { return a.name==name; } );
            } //!< Iterator to AtomData with `name` -- `end()` if not found

        };

    /**
     * @brief Class for specifying changes to Space
     *
     * - If `moved` or `removed` are defined for a group, but are
     *   empty, it is assumed that *all* particles in the group are affected.
     * - The total volume change is defined as `|dV| = dV.norm()`.
     */
    template<class Tpvec>
        struct Change
        {
            typedef int groupindex; //!< Group index in group vector
            typedef int atomindex;  //!< Atom index in particle vector
            typedef int molid;      //!< Molecule id (unique)
            Point dV = {0,0,0};     //!< Volume change (in different directions)
            std::map<groupindex, std::vector<atomindex>> moved;  //!< Moved particles, grouped by group index
            std::map<groupindex, std::vector<atomindex>> removed;//!< Removed particles, grouped by group index
            std::map<molid, Tpvec> inserted; //!< Particles to be inserted, grouped by *molecule id*

            void clear()
            {
                dV.setZero();
                moved.clear();
                removed.clear();
                inserted.clear();
                assert(empty());
            } //!< Clear all change data

            bool empty() const
            {
                if ( moved.empty())
                    if ( removed.empty())
                        if ( inserted.empty())
                            if ( dV.squaredNorm()<1e-9 )
                                return true;
                return false;
            } //!< Check if change object is empty
        };


}//namespace

using namespace Faunus;
using namespace std;

int main() {
    using DipoleParticle = Particle<Radius, Dipole, Cigar>;
    using PointParticle = Particle<>;
    typedef DipoleParticle Tparticle;

    Tparticle p, dst;
    p.charge = 1.0;
    p.sclen=10.2;
    p.mulen = -0.5;
    p.radius = 999;

    Tjson j = p;

    std::cout << j << "\n";

    dst = j;

    std::cout << Tjson(dst) << "\n";

    Geometry::Cuboid cuboid;
    Geometry::Cylinder cyl;

    p.pos = {1,0,1};
    dst.pos = {0,0,0};

    cout << "d = " << cuboid.vdist( dst.pos, p.pos ) << endl;

    Eigen::Quaterniond q;
    Eigen::Matrix3d m;

    p.rotate(q, m);
    vector<decltype(p)> vec;
    vec.push_back(p);
    Geometry::translate(vec, {0,0,0});

    Random r = R"( {"randomseed" : "hardware"} )"_json;

    cout << r() << endl;

    Group rr(0,11);

    for (auto i : rr)
        cout << i << endl;

    cout << rr.front() << " - " << rr.back() << endl;

    //analyse<Tparticle> a;
    //a.sample(p);
}

