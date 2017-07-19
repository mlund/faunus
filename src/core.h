#pragma once

#include <vector>
#include <map>
#include <tuple>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <string>
#include <cmath>
#include <random>
#include <Eigen/Geometry>
#include <json.hpp>
#include <range/v3/all.hpp>

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

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("Check unit conversion and string literals")
    {
        using doctest::Approx;
        pc::temperature = 298.15_K;
        CHECK( 1.0e-10_m == 1 );
        CHECK( (1/1.0_Debye) == Approx(4.8032) );
        CHECK( 1.0_Debye == Approx( 3.33564e-30_Cm ) );
        CHECK( 1.0_Debye == Approx( 0.20819434_eA ) );
        CHECK( 360.0_deg == Approx( 2*std::acos(-1) ) );
        CHECK( (1.0_mol / 1.0_liter) == Approx( 1.0_molar ) );
        CHECK( 1.0_bar == Approx( 0.987_atm ) );
        CHECK( 1.0_atm == Approx( 101325._Pa) );
        CHECK( 1.0_kT == Approx( 2.47897_kJmol ) );
        CHECK( 1.0_hartree == Approx( 2625.499_kJmol ) );
    }
#endif

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
     *
     * @warning
     * Random number seeding in C++11 has issues, probably not important for molecular
     * simulations. More here:
     * https://probablydance.com/2016/12/29/random_seed_seq-a-small-utility-to-properly-seed-random-number-generators-in-c/
     */
    struct Random {
        std::mt19937 engine; //!< Random number engine used for all operations
        std::uniform_real_distribution<double> dist01; //!< Uniform real distribution [0,1)

        Random() : dist01(0,1) {}

        void seed() { engine = std::mt19937(std::random_device()()); }

        double operator()() { return dist01(engine); } //!< Double in uniform range [0,1)

        int range( int min, int max )
        {
            std::uniform_int_distribution<int> d(min, max);
            return d(engine);
        } //!< Integer in uniform range [min:max]

        template<class Titer>
            Titer sample(const Titer &beg, const Titer &end)
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

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("Random")
    {
        Random slump;
        int min=10, max=0, N=1e6;
        double x=0;
        for (int i=0; i<N; i++) {
            int j = slump.range(0,9);
            if (j<min) min=j;
            if (j>max) max=j;
            x+=j;
        }
        CHECK( min==0 );
        CHECK( max==9 );
        CHECK( std::fabs(x/N) == doctest::Approx(4.5).epsilon(0.01) );

        Random r1 = R"( {"randomseed" : "hardware"} )"_json; // non-deterministic seed
        Random r2; // default is a deterministic seed
        CHECK( r1() != r2() );
        Random r3 = Tjson(r1); // r1 --> json --> r3
        CHECK( r1() == r3() );

        // check if random_device works
        Random a, b;
        CHECK( a() == b() );
        a.seed();
        b.seed();
        CHECK( a() != b() );
    }
#endif

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

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] ranunit") {
    }
#endif

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
        void from_json(const Tjson& j) override { radius = j.value("r", 0.0); }
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
            j["id"]=a.id;
            j["pos"] = a.pos; 
            to_json<Properties...>(j, Properties(a)... );
        }

    template<typename... Properties>
        void from_json(const Tjson& j, Particle<Properties...> &a) {
            a.charge = j.value<double>("q", 0.0); 
            a.mw = j.value<double>("mw", 0.0); 
            a.id = j.value<int>("id",-1);
            a.pos = j.value("pos", Point(0,0,0));
            from_json<Properties...>(j, dynamic_cast<Properties&>(a)...);
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("Particle") {
        Particle<Dipole,Cigar,Radius> p1, p2;
        p1.id = 100;
        p1.mw = 0.2;
        p1.pos = {1,2,3};
        p1.charge = -0.8;
        p1.radius = 7.1;
        p1.mu = {0,0,1};
        p1.mulen = 2.8;
        p1.scdir = {0.1, 0.3, 1.9};
        p1.sclen = 0.5;

        p2 = Tjson(p1); // p1 --> json --> p2
        CHECK( Tjson(p1) == Tjson(p2) ); // identical json objects?
        CHECK( p2.id == 100 );
    }
#endif

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

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("Cylinder") {
            Cylinder c;
            c.setRadius( 1.0, 1/pc::pi );
            CHECK( c.getVolume() == doctest::Approx( 1.0 ) );
        }
#endif

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


    template<typename T>
        void inline swap_to_back(T first, T last, T end) {
            while (end-- > last)
                std::iter_swap(first++,end);
        } //!< Move range [first:last] to [end] by swapping elements

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] swap_to_back") {
        typedef std::vector<int> _T;
        _T v = {1,2,3,4};

        swap_to_back( v.begin(), v.end(), v.end() );
        CHECK( v==_T({1,2,3,4}) );

        std::sort(v.begin(), v.end());
        swap_to_back( v.begin()+1, v.begin()+3, v.end() );
        CHECK( v==_T({1,4,3,2}) );
    }
#endif

    template<class T>
        struct IterRange : std::pair<T,T> { 
            using std::pair<T,T>::pair;
            T& begin() { return this->first; }
            T& end() { return this->second; }
            const T& begin() const { return this->first; }
            const T& end() const { return this->second; }
            size_t size() const { return std::distance(this->first, this->second); }
            void resize(size_t n) { end() += n; }
            bool empty() const { return this->first==this->second; }
            void clear() { this->second = this->first; }
            std::pair<int,int> to_index(T reference) {
                return {std::distance(reference, begin()), std::distance(reference, end())};
            } //!< Returns index pair relative to reference
        }; //!< Turns a pair of iterators into a range

    /**
     * @brief Turns a pair of iterators into an elastic range
     *
     * The elastic range is a range where elements can be deactivated
     * and later activated without inserting/erasing:
     *
     * - Just deactivated elements are moved to `end()` and can be retrieved from there.
     * - Just activated elements are placed at `end()-n`.
     * - The true size the range is given by `capacity()`
     */
    template<class T>
        class ElasticRange : public IterRange<T> {
            private:
                T _trueend;
            public:
                typedef typename T::value_type Tvalue;
                using IterRange<T>::begin;
                using IterRange<T>::end;
                using IterRange<T>::size;

                ElasticRange(T begin, T end) : IterRange<T>({begin,end}), _trueend(end) {}

                size_t capacity() const { return std::distance(begin(), _trueend); }

                auto inactive() const { return IterRange<T>({ end(), _trueend}); }

                void deactivate(T first, T last) {
                    size_t n = std::distance(first,last);
                    assert(n>=0);
                    assert(first>=begin() && last<=end() );
                    std::rotate( begin(), last, end() );
                    end() -= n;
                    assert(size() + inactive().size() == capacity());
                } //!< Deactivate particles by moving to end, reducing the effective size

                void activate(T first, T last) {
                    size_t n = std::distance(first,last);
                    std::rotate( end(), first, _trueend );
                    end() += n;
                    assert(size() + inactive().size() == capacity());
                } //!< Activate previously deactivated elements
        };

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("ElasticRange") {
        std::vector<int> v = {10,20,30,40,50,60};
        ElasticRange< typename std::vector<int>::iterator> r(v.begin(), v.end());
        CHECK( r.size() == 6 );
        CHECK( r.empty() == false );
        CHECK( r.size() == r.capacity() );
        *r.begin() += 1;
        CHECK( v[0]==11 );

        r.deactivate( r.begin(), r.end() ); 
        CHECK( r.size() == 0 );
        CHECK( r.empty() == true );
        CHECK( r.capacity() == 6 );
        CHECK( r.begin() == r.end() );

        r.activate( r.inactive().begin(), r.inactive().end() );
        CHECK( r.size() == 6);
        CHECK( std::is_sorted(r.begin(), r.end() ) == true ); // back to original

        r.deactivate( r.begin()+1, r.begin()+3 ); 
        CHECK( r.size() == 4 );
        CHECK( std::find(r.begin(), r.end(), 20)==r.end() );
        CHECK( std::find(r.begin(), r.end(), 30)==r.end() );
        CHECK( *r.end()==20); // deactivated elements can be retrieved from `end()`
        CHECK( *(r.end()+1)==30);

        auto ipair = r.to_index( v.begin() );
        CHECK( ipair.first == 0);
        CHECK( ipair.second == 4);

        r.activate( r.end(), r.end()+2 );
        CHECK( *(r.end()-2)==20); // activated elements can be retrieved from `end()-n`
        CHECK( *(r.end()-1)==30);
        CHECK( r.size() == 6);
    }
#endif

    template<class Titer>
        struct Group : public ElasticRange<Titer> {
            typedef ElasticRange<Titer> base;
            typedef typename base::Tvalue Tvalue;
            int id=-1;           //!< Type id
            bool atomic=false;   //!< Is it an atomic group?
            Point cm={0,0,0};    //!< Mass center
            Group(Titer begin, Titer end) : base(begin,end) {}

            auto filter( std::function<bool(Tvalue&)> f ) const {
                return *this | ranges::view::filter(f); 
            }

            Group& operator=(const Group &o) {
                if (this->capacity() != o.capacity())
                    throw std::runtime_error("capacity mismatch");
                this->resize(o.size());
                id = o.id;
                atomic = o.atomic;
                cm = o.cm;
                return *this;
            }
        };

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Group") {
        struct particle { int id; };
        std::vector<particle> p = { {1}, {2}, {1}, {4}, {1} };
        Group<typename std::vector<particle>::iterator> g(p.begin(), p.end());

        auto slice = g.filter( [](particle &i){return i.id==1;} );

        CHECK( std::distance(slice.begin(), slice.end()) == 3 );
    }
#endif

    /**
     * @brief Class for specifying changes to Space
     *
     * - If `moved` or `removed` are defined for a group, but are
     *   empty, it is assumed that *all* particles in the group are affected.
     */
    template<class Tpvec>
        struct Change
        {

            double dV = 0;     //!< Volume change (in different directions)

            struct data {
                typedef std::pair<int,int> Tindexpair;
                bool all=false;
                std::vector<int> index;
                std::vector<Tindexpair> activated, deactivated;
            };

            std::map<int, data> groups;

            void clear()
            {
                dV = 0;
                groups.clear();
                assert(empty());
            } //!< Clear all change data

            bool empty() const
            {
                if ( groups.empty())
                    if ( std::fabs(dV)>0 )
                        return true;
                return false;
            } //!< Check if change object is empty
        };

    template<class Tgeometry, class Tparticle>
        struct Space {

            typedef Space<Tgeometry,Tparticle> Tspace;
            typedef std::vector<Tparticle> Tpvec;

            typedef Group<typename Tpvec::iterator> Tgroup;
            typedef std::vector<Tgroup> Tgvec;

            typedef Change<Tpvec> Tchange;

            typedef std::function<void(Tspace&, const Tchange&)> ChangeTrigger;
            typedef std::function<void(Tspace&, const Tspace&, const Tchange&)> SyncTrigger;

            std::vector<ChangeTrigger> changeTriggers; //!< Call when a Change object is applied
            std::vector<SyncTrigger> onSyncTriggers;   //!< Call when two Space objects are synched

            Tpvec p;       //!< Particle vector
            Tgvec groups;  //!< Group vector
            Tgeometry geo; //!< Container geometry

            template<typename T>
                std::vector<int> findIndex(const std::vector<T> &v, int id) const {
                    std::vector<int> ndx;
                    ndx.reserve( v.size() );
                    for (size_t i=0; i<v.size(); i++)
                        if (v[i].id==id)
                            ndx.push_back(i);
                    return ndx;
                } //!< ordered vector of index pointing to `v`-elements of type `id`

            template<typename T>
                int random(const std::vector<T> &v, int id, Random &r) const {
                    auto i = findIndex(v, id);
                    if (i.empty())
                        return -1;
                    return *r.sample(i.begin(), i.end());
                } //!< find index to random element of `id`. Returns -1 if not found.

            void sync(const Tspace &o, const Tchange &change) {

                for (auto &m : change.groups) {

                    auto &gold = groups[m.first];
                    auto &gnew = o.groups[m.first];

                    gold = gnew;

                    // whole group has moved
                    if (m.second.all) {
                        auto it = gold.begin();
                        for (auto &i : gnew)
                            *it++ = i;
                        assert(it==gold.end());
                    }
                }

                    /*
                    for (auto &m : change.groups) {
                        int n=0; // total number of activated particles
                        auto &g = groups[m.first];
                        assert( g.size() < o.groups[m.first].size() );
                        for (auto &rng : m.second) {
                            auto begin = p.begin()+rng.first;
                            auto end = p.begin()+rng.second;
                            n += end-begin;
                            g.activate(begin, end);
                        }
                        assert( g.size() == o.groups[m.first].size() );

                        // activated elements end up just below end()
                        int end = std::distance(p, g.end());
                        for (int i=end-n; i<end; i++)
                            p[i] = o.p[i];
                    }
                    if (std::fabs(change.dV)>0)
                        geo = o.geo;

                    for (auto& m : change.moved) { // loop over moved groups
                        groups[ m.first ] = o.groups[ m.first ];
                        for (auto i : m.second)
                            p[i] = o.p[i];

                    }*/

            } //!< Copy differing data from other (o) Space using Change object

            void applyChange(const Tchange &change) {
                for (auto& f : changeTriggers)
                    f(*this, change);
            }

            template<class Tindex>
                void erase(int molid, Tindex ndx );

            //auto findAtomData(const std::string &name) const {
            //    return std::find_if( atoms.begin(), atoms.end(),
            //            [&name](const AtomData &a) { return a.name==name; } );
            //} //!< Iterator to AtomData with `name` -- `end()` if not found
        };

}//end of faunus namespace

