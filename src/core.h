#pragma once

#include <limits>
#include <vector>
#include <map>
#include <tuple>
#include <iostream>
#include <sstream>
#include <fstream>
#include <type_traits>
#include <string>
#include <cmath>
#include <random>
#include <memory>
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
    template<class T>
        void from_json(const nlohmann::json& j, Eigen::Matrix<T,3,1> &p) {
            if ( j.size()==3 ) {
                int i=0;
                for (auto d : j.get<std::vector<T>>())
                    p[i++] = d;
                return;
            }
            throw std::runtime_error("JSON->Eigen conversion error");
        }
};

/** @brief Faunus main namespace */
namespace Faunus {

    typedef Eigen::Vector3d Point; //!< 3d vector
    typedef nlohmann::json json;  //!< Json object

    using std::cout;
    using std::endl;
    using std::fabs;
    using std::exp;
    using std::sqrt;
    using std::log;

    template<class T1, class T2>
        int distance(T1 first, T2 last) {
            return std::distance( &(*first), &(*last) );
        } //!< Distance between two arbitrary contiguous iterators

    template<class T>
        int size(T &rng) {
            return ranges::distance(rng.begin(), rng.end());
        } //!< Size of arbitrary range

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] distance")
    {
        std::vector<long long int> v = {10,20,30,40,30};
        auto rng = v | ranges::view::filter( [](int i){return i==30;} );
        CHECK( Faunus::distance(v.begin(), rng.begin()) == 2 );
        auto it = rng.begin();
        CHECK( Faunus::distance(v.begin(), ++it) == 4 );
    }
#endif

    inline json merge( const json &a, const json &b ) {
        json result = a.flatten();
        json tmp = b.flatten();
        for ( auto it = tmp.begin(); it != tmp.end(); ++it )
            result[it.key()] = it.value();
        return result.unflatten();
    } //!< Merge two json objects

    inline json openjson( const std::string &file ) {
        json js;
        std::ifstream f (file );
        if ( f ) {
            try {
                f >> js;
            }
            catch(std::exception& e) {
                throw std::runtime_error("Syntax error in JSON file " + file + ": " + e.what());
            }
        }
        else
            throw std::runtime_error("Cannot find or read JSON file " + file);
        return js;
    } //!< Read json file into json object (w. syntax check)

    double _round(double x, int n=3) {
        std::stringstream o;
        o << std::setprecision(n) << x;
        return std::stod(o.str());
    } //!< Round to n number of significant digits

    void _roundjson(json &j, int n=3) {
        if (j.is_object())
            for (auto &i : j)
                if (i.is_number_float())
                    i = _round(i,n);
    } // round float objects to n number of significant digits

    double value_inf(const json &j, const std::string &key) {
        auto it = j.find(key);
        if (it==j.end())
            throw std::runtime_error("unknown json key '" + key + "'");
        else
            if (it->is_string()) {
                if (*it=="inf")
                    return std::numeric_limits<double>::infinity();
                if (*it=="-inf")
                    return -std::numeric_limits<double>::infinity();
                throw std::runtime_error("value must be number or 'inf'");
            }
        return double(*it);
    } //!< Extract floating point from json and allow for 'inf' and '-inf'

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
        static T temperature=300; //!< Temperature (Kelvin)
        static inline T kT() { return temperature*kB; } //!< Thermal energy (Joule)
        static inline T lB( T epsilon_r ) {
            return e*e/(4*pi*e0*epsilon_r*1e-10*kT());
        } //!< Bjerrum length (angstrom)
        static inline T lB2epsr( T bjerrumlength ) {
            return lB(bjerrumlength);
        } //!< lB --> relative dielectric constant
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
        constexpr T operator "" _angstrom3( T l ) { return l; } //!< Angstrom^3 to Angstrom^3
        constexpr T operator "" _m( T l ) { return l * 1e10; } //!< Meter to angstrom
        constexpr T operator "" _bohr( T l ) { return l * 0.52917721092; } //!< Bohr to angstrom
        constexpr T operator "" _nm( T l ) { return l * 10; } //!< nanometers to angstrom
        constexpr T operator "" _liter( T v ) { return v * 1e27; } //!< liter to angstrom cubed
        constexpr T operator "" _m3( T v ) { return v * 1e30; } //!< -> cubic meter to angstrom cubed
        constexpr T operator "" _mol( T n ) { return n * pc::Nav; } //!< moles to particles
        constexpr T operator "" _molar( T c ) { return c * 1.0_mol / 1.0_liter; } //!< molar to particle / angstrom^3
        constexpr T operator "" _mM( T c ) { return c * 1.0e-3_mol / 1.0_liter; } //!< millimolar to particle / angstrom^3
        constexpr T operator "" _rad( T a ) { return a; } //!< Radians to radians
        constexpr T operator "" _deg( T a ) { return a * pc::pi / 180; } //!< Degrees to radians
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
    TEST_CASE("[Faunus] Units and string literals")
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

    namespace u8 {
        const std::string angstrom = "\u00c5";   //!< Angstrom symbol
        const std::string beta = "\u03b2";       //!< Greek beta
        const std::string cubed = "\u00b3";      //!< Superscript 3
        const std::string cuberoot = "\u221b";   //!< Cubic root
        const std::string degrees = "\u00b0";    //!< Degrees
        const std::string Delta = "\u0394";      //!< Greek Delta
        const std::string epsilon = "\u03f5";    //!< Greek epsilon
        const std::string epsilon_m = "\u03b5";  //!< Greek epsilon (minuscule)
        const std::string gamma = "\u0263";      //!< Greek gamma
        const std::string Gamma = "\u0393";      //!< Greek capital gamma
        const std::string infinity="\u221E";     //!< Infinity
        const std::string kappa = "\u03ba";      //!< Greek kappa
        const std::string mu = "\u03bc";         //!< Greek mu
        const std::string partial = "\u2202";    //!< Partial derivative
        const std::string percent = "\ufe6a";    //!< Percent sign
        const std::string pm = "\u00b1";         //!< Plus minus sign
        const std::string rho = "\u03C1";        //!< Greek rho
        const std::string rootof = "\u221a";     //!< Square root sign
        const std::string squared = "\u00b2";    //!< Superscript 2
        const std::string sigma = "\u03c3";      //!< Greek sigma
        const std::string superminus = "\u207b"; //!< Superscript minus (-)
        const std::string subr = "\u1D63";       //!< Subscript "r"
        const std::string theta = "\u03b8";      //!< Greek theta

        inline std::string bracket( const std::string &s )
        {
            return "\u27e8" + s + "\u27e9";
        }
    } //!< Unicode

    struct Tensor : public Eigen::Matrix3d {
        typedef Eigen::Matrix3d base;

        Tensor() {
            base::setZero();
        } //!< Constructor, clear data

        Tensor( double xx, double xy, double xz, double yy, double yz, double zz ) {
            (*this) << xx, xy, xz, xy, yy, yz, xz, yz, zz;
        } //!< Construct from input

        template<typename T>
            Tensor( const Eigen::MatrixBase<T> &other ) : base(other) {}

        template<typename T>
            Tensor &operator=( const Eigen::MatrixBase<T> &other )
            {
                base::operator=(other);
                return *this;
            }

        void rotate( const base &m ) {
            (*this) = m * (*this) * m.transpose();
        } //!< Rotate using rotation matrix. Remove?
        void eye() { *this = base::Identity(3, 3); }
    }; //!< Tensor class

    void to_json(nlohmann::json& j, const Tensor &t) {
        j = { t(0,0), t(0,1), t(0,2), t(1,1), t(1,2), t(2,2) };
    } //!< Tensor -> Json

    void from_json(const nlohmann::json& j, Tensor &t) {
        if ( j.size()!=6 || !j.is_array() )
            throw std::runtime_error("Json->Tensor: array w. exactly six coefficients expected.");
        t = Tensor(j[0],j[1],j[2],j[3],j[4],j[5]);
    } //!< Json -> Tensor

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Tensor") {
        using doctest::Approx;
        Tensor Q1 = Tensor(1,2,3,4,5,6);
        Tensor Q2 = json(Q1); // Q1 --> json --> Q2
        CHECK( json(Q1) == json(Q2) ); // Q1 --> json == json <-- Q2 ?
        CHECK( Q2 == Tensor(1,2,3,4,5,6) );

        auto m = Eigen::AngleAxisd( pc::pi/2, Point(0,1,0) ).toRotationMatrix();
        Q1.rotate(m);
        CHECK( Q1(0,0) == Approx(6) );
        CHECK( Q1(0,1) == Approx(5) );
        CHECK( Q1(0,2) == Approx(-3) );
        CHECK( Q1(1,1) == Approx(4) );
        CHECK( Q1(1,2) == Approx(-2) );
        CHECK( Q1(2,2) == Approx(1) );
    }
#endif

    /**
     * Example code:
     *
     * ```{.cpp}
     *     Random r1;                                     // default deterministic seed
     *     Random r2 = json(r1);                          // copy engine state
     *     Random r3 = R"( {"seed" : "hardware"} )"_json; // non-deterministic seed
     *     Random r1.seed();                              // non-deterministic seed
     * ```
     */
    struct Random {
        std::mt19937 engine; //!< Random number engine used for all operations
        std::uniform_real_distribution<double> dist01; //!< Uniform real distribution [0,1)

        inline Random() : dist01(0,1) {}

        inline void seed() { engine = std::mt19937(std::random_device()()); }

        inline double operator()() { return dist01(engine); } //!< Double in uniform range [0,1)

        inline int range( int min, int max )
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
    }; //!< Class for handling random number generation

    void to_json(json &j, const Random &r) {
        std::ostringstream o;
        o << r.engine;
        j["seed"] = o.str();
    } //!< Random to json conversion

    void from_json(const json &j, Random &r) {
        if (j.is_object()) {
            auto seed = j.value("seed", std::string());
            try {
                if (seed=="default" || seed=="fixed")
                    return;
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
    } //!< json to Random conversion

    static Random random; // global instance of Random

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Random")
    {
        Random slump; // local instance

        CHECK( slump() == random() );

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

        Random r1 = R"( {"seed" : "hardware"} )"_json; // non-deterministic seed
        Random r2; // default is a deterministic seed
        CHECK( r1() != r2() );
        Random r3 = json(r1); // r1 --> json --> r3
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
     * @brief Convert cartesian- to spherical-coordinates
     * @note Input (x,y,z), output \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [-\pi,\pi) \f$, and \f$ \phi\in [0,\pi] \f$.
     */
    inline Point xyz2rtp(const Point &p, const Point &origin={0,0,0}) {
        Point xyz = p - origin;
        double radius = xyz.norm();
        return {
            radius,
                std::atan2( xyz.y(), xyz.x() ),
                std::acos( xyz.z()/radius) };
    }

    /**
     * @brief Convert spherical- to cartesian-coordinates
     * @param origin The origin to be added (optional)
     *
     * @note Input \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [0,2\pi) \f$, and \f$ \phi\in [0,\pi] \f$, and output (x,y,z).
     */
    inline Point rtp2xyz(const Point &rtp, const Point &origin = {0,0,0}) {
        return origin + rtp.x() * Point(
                std::cos(rtp.y()) * std::sin(rtp.z()),
                std::sin(rtp.y()) * std::sin(rtp.z()),
                std::cos(rtp.z()) );
    }

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] spherical coordinates") {
        using doctest::Approx;

        Point sph1 = {2, 0.5, -0.3};
        auto pnt1 = rtp2xyz(sph1); // sph --> cart
        auto sph2 = xyz2rtp(pnt1); // cart --> sph

        CHECK( pnt1.norm() == Approx(2));
        CHECK( sph1.x() == Approx(sph2.x()));
        //CHECK( sph1.y() == Approx(sph2.y()));
        //CHECK( sph1.z() == Approx(sph2.z()));
    }
#endif

    Point ranunit_neuman(Random &rand)
    {
        double r2;
        Point p;
        do {
            p = {rand()-0.5, rand()-0.5, rand()-0.5};
            r2 = p.squaredNorm();
        } while ( r2 > 0.25 );
        return p / std::sqrt(r2);
    } //!< Random unit vector using Neuman's method ("sphere picking")

    Point ranunit_polar(Random &rand) {
        return rtp2xyz( {1, 2*pc::pi*rand(), std::acos(2*rand()-1)} );
    } //!< Random unit vector using polar coordinates ("sphere picking")

    const auto& ranunit = ranunit_polar; //!< Alias for default random unit vector function

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] ranunit_neuman") {
        Random r;
        int n=2e5;
        Point rtp(0,0,0);
        for (int i=0; i<n; i++)
            rtp += xyz2rtp( ranunit_neuman(r) );
        rtp = rtp / n;
        CHECK( rtp.x() == doctest::Approx(1) );
        CHECK( rtp.y() == doctest::Approx(0).epsilon(0.005) ); // theta [-pi:pi] --> <theta>=0
        CHECK( rtp.z() == doctest::Approx(pc::pi/2).epsilon(0.005) );// phi [0:pi] --> <phi>=pi/2
    }
#endif

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] ranunit_polar") {
        Random r;
        int n=2e5;
        Point rtp(0,0,0);
        for (int i=0; i<n; i++)
            rtp += xyz2rtp( ranunit_polar(r) );
        rtp = rtp / n;
        CHECK( rtp.x() == doctest::Approx(1) );
        CHECK( rtp.y() == doctest::Approx(0).epsilon(0.005) ); // theta [-pi:pi] --> <theta>=0
        CHECK( rtp.z() == doctest::Approx(pc::pi/2).epsilon(0.005) );// phi [0:pi] --> <phi>=pi/2
    }
#endif

    /**
     * @brief Quaternion rotation routine using the Eigen library
     */
    struct QuaternionRotate : public std::pair<Eigen::Quaterniond, Eigen::Matrix3d> {

        typedef std::pair<Eigen::Quaterniond, Eigen::Matrix3d> base;
        //using base::pair;
        using base::first;
        using base::second;

        double angle=0; //!< Rotation angle

        QuaternionRotate() {};

        QuaternionRotate(double angle, Point u) { set(angle,u); };

        void set(double angle, Point u) {
            this->angle = angle;
            u.normalize();
            first = Eigen::AngleAxisd(angle, u);
            second << 0, -u.z(), u.y(), u.z(), 0, -u.x(), -u.y(), u.x(), 0;
            second =
                Eigen::Matrix3d::Identity() + second * std::sin(angle)
                + second * second * (1 - std::cos(angle));

            // Quaternion can be converted to rotation matrix:
            // second = first.toRotationMatrix()
        }

        Point operator()( Point a, std::function<void(Point&)> boundary = [](Point &i){}, const Point &shift={0,0,0} ) const
        {
            a = a - shift;
            boundary(a);
            a = first * a + shift;
            boundary(a);
            return a;
            // https://www.cc.gatech.edu/classes/AY2015/cs4496_spring/Eigen.html
        } //!< Rotate point w. optional PBC boundaries

        auto operator()( const Eigen::Matrix3d &a ) const
        {
            return second * a * second.transpose();
        } //!< Rotate matrix/tensor
    };

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] QuaternionRotate")
    {
        using doctest::Approx;
        QuaternionRotate qrot;
        Point a = {1,0,0};
        qrot.set( pc::pi/2, {0,1,0} ); // rotate around y-axis
        CHECK( qrot.angle == Approx(pc::pi/2) );

        SUBCASE("rotate using quaternion") {
            a = qrot(a); // rot. 90 deg.
            CHECK( a.x() == Approx(0) );
            a = qrot(a); // rot. 90 deg.
            CHECK( a.x() == Approx(-1) );
        }

        SUBCASE("rotate using rotation matrix") {
            a = qrot.second * a;
            CHECK( a.x() == Approx(0) );
            a = qrot.second * a;
            CHECK( a.x() == Approx(-1) );
        }
    }
#endif

    /** @brief Base class for particle properties */
    struct ParticlePropertyBase {
        virtual void to_json(json &j) const=0; //!< Convert to JSON object
        virtual void from_json(const json &j)=0; //!< Convert from JSON object
        inline void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d&) {}
    };

    template <typename... Ts>
        auto to_json(json&) -> typename std::enable_if<sizeof...(Ts) == 0>::type {} //!< Particle to JSON
    template<typename T, typename... Ts>
        void to_json(json& j, const ParticlePropertyBase &a, const Ts&... rest) {
            a.to_json(j);
            to_json<Ts...>(j, rest...);
        } //!< Particle to JSON

    // JSON --> Particle
    template <typename... Ts>
        auto from_json(const json&) -> typename std::enable_if<sizeof...(Ts) == 0>::type {}
    template<typename T, typename... Ts>
        void from_json(const json& j, ParticlePropertyBase &a, Ts&... rest) {
            a.from_json(j);
            from_json<Ts...>(j, rest...);
        }

    struct Radius : public ParticlePropertyBase {
        double radius=0; //!< Particle radius
        void to_json(json &j) const override { j["r"] = radius; }
        void from_json(const json& j) override { radius = j.value("r", 0.0); }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    }; //!< Radius property

    struct Charge : public ParticlePropertyBase {
        double charge=0; //!< Particle radius
        void to_json(json &j) const override { j["q"] = charge; }
        void from_json(const json& j) override { charge = j.value("q", 0.0); }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    }; //!< Charge (monopole) property

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
            mu = q * mu;
        } //!< Rotate dipole moment

        void to_json(json &j) const override {
            j["mulen"] = mulen ;
            j["mu"] = mu;
        }

        void from_json(const json& j) override {
            mulen = j.value("mulen", 0.0);
            mu = j.value("mu", Point(1,0,0) );
        }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    struct Quadrupole : public ParticlePropertyBase {
        Tensor Q;      //!< Quadrupole
        void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) { Q.rotate(m); } //!< Rotate dipole moment
        void to_json(json &j) const override { j["Q"] = Q; }
        void from_json(const json& j) override { Q = j.value("Q", Q); }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    }; // Quadrupole property

    struct Cigar : public ParticlePropertyBase {
        double sclen=0;       //!< Sphero-cylinder length
        Point scdir = {1,0,0};//!< Sphero-cylinder direction unit vector

        void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d&) {
            scdir = q * scdir;
        } //!< Rotate sphero-cylinder

        void to_json(json &j) const override {
            j["sclen"] = sclen;
            j["scdir"] = scdir;
        }
        void from_json(const json& j) override {
            sclen = j.value("sclen", 0.0);
            scdir = j.value("scdir", Point(1,0,0) );
        }
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    }; //!< Sphero-cylinder properties

    /** @brief Particle
     *
     * In addition to basic particle properties (id, charge, weight), arbitrary
     * features can be added via template arguments:
     *
     * ```{.cpp}
     *     Particle<Dipole> p
     *     p = R"( {"mulen":2.8} )"_json; // Json --> Particle
     *     std::cout << json(p); // Particle --> Json --> stream
     * ```
     *
     * Currently the following properties are implemented:
     *
     * Keyword   | Class       |  Description
     * --------- | ----------- | --------------------
     *  `id`     | `Particle`  |  Type id (int)
     *  `pos`    | `Particle`  |  Position (Point)
     *  `mw`     | `Particle`  |  Weight (double)
     *  `q`      | `Charge`    |  valency (e)
     *  `r`      | `Radius`    |  radius (angstrom)
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
                        a.rotate(q,m);
                        _rotate<Ts...>(q, m, rest...);
                    }

            public:
                int id=-1;         //!< Particle id/type
                Point pos={0,0,0}; //!< Particle position vector
                double mw=1;       //!< Molecular weight

                Particle() : Properties()... {}

                void rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) {
                    _rotate<Properties...>(q, m, dynamic_cast<Properties&>(*this)...);
                } //!< Rotate all internal coordinates if needed

                EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };

    template<typename... Properties>
        void to_json(json& j, const Particle<Properties...> &a) {
            j["id"]=a.id;
            j["mw"]= a.mw;
            j["pos"] = a.pos;
            to_json<Properties...>(j, Properties(a)... );
        }

    template<typename... Properties>
        void from_json(const json& j, Particle<Properties...> &a) {
            a.id = j.value("id", a.id);
            a.mw = j.value("mw", a.mw);
            a.pos = j.value("pos", a.pos);
            from_json<Properties...>(j, dynamic_cast<Properties&>(a)...);
        }

    using ParticleAllProperties = Particle<Radius,Dipole,Charge,Quadrupole,Cigar>;

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Particle") {
        using doctest::Approx;
        ParticleAllProperties p1, p2;
        p1.id = 100;
        p1.pos = {1,2,3};
        p1.charge = -0.8;
        p1.radius = 7.1;
        p1.mu = {0,0,1};
        p1.mulen = 2.8;
        p1.scdir = {-0.1, 0.3, 1.9};
        p1.sclen = 0.5;
        p1.Q = Tensor(1,2,3,4,5,6);

        p2 = json(p1); // p1 --> json --> p2

        CHECK( json(p1) == json(p2) ); // p1 --> json == json <-- p2 ?

        CHECK( p2.id == 100 );
        CHECK( p2.pos == Point(1,2,3) );
        CHECK( p2.charge == -0.8 );
        CHECK( p2.radius == 7.1 );
        CHECK( p2.mu == Point(0,0,1) );
        CHECK( p2.mulen == 2.8 );
        CHECK( p2.scdir == Point(-0.1, 0.3, 1.9) );
        CHECK( p2.sclen == 0.5 );
        CHECK( p2.Q == Tensor(1,2,3,4,5,6) );

        // check of all properties are rotated
        QuaternionRotate qrot( pc::pi/2, {0,1,0} );
        p1.mu = p1.scdir = {1,0,0};
        p1.rotate( qrot.first, qrot.second );

        CHECK( p1.mu.x() == Approx(0) );
        CHECK( p1.mu.z() == Approx(-1) );
        CHECK( p1.scdir.x() == Approx(0) );
        CHECK( p1.scdir.z() == Approx(-1) );

        CHECK( p1.Q(0,0) == Approx(6) );
        CHECK( p1.Q(0,1) == Approx(5) );
        CHECK( p1.Q(0,2) == Approx(-3) );
        CHECK( p1.Q(1,1) == Approx(4) );
        CHECK( p1.Q(1,2) == Approx(-2) );
        CHECK( p1.Q(2,2) == Approx(1) );
    }
#endif

    /**
     * @brief Eigen::Map facade to data members in STL container
     *
     * No data is copied and modifications of the Eigen object
     * modifies the original container and vice versa.
     *
     * Example:
     *
     *    std::vector<Tparticle> v(10);
     *    auto m1 = asEigenVector(v.begin, v.end(), &Tparticle::pos);    --> 10x3 maxtix view
     *    auto m2 = asEigenMatrix(v.begin, v.end(), &Tparticle::charge); --> 10x1 vector view
     *
     * @warning Be careful that objects are properly aligned and divisible with `sizeof<double>`
     */
    template<typename dbl=double, class iter, class memberptr>
        auto asEigenMatrix(iter begin, iter end, memberptr m) {
            typedef typename std::iterator_traits<iter>::value_type T;
            static_assert( sizeof(T) % sizeof(dbl) == 0, "value_type size must multiples of double");
            const size_t s = sizeof(T) / sizeof(dbl);
            const size_t cols = sizeof(((T *) 0)->*m) / sizeof(dbl);
            typedef Eigen::Matrix<dbl, Eigen::Dynamic, cols> Tmatrix;
            return Eigen::Map<Tmatrix, 0, Eigen::Stride<1,s>>((dbl*)&(*begin.*m), end-begin, cols).array();
        }

    template<typename dbl=double, class iter, class memberptr>
        auto asEigenVector(iter begin, iter end, memberptr m) {
            typedef typename std::iterator_traits<iter>::value_type T;
            static_assert( std::is_same<dbl&, decltype(((T *) 0)->*m)>::value, "member must be a scalar");
            return asEigenMatrix<dbl>(begin, end, m).col(0);
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] asEigenMatrix") {
        using doctest::Approx;
        typedef Particle<Radius, Charge, Dipole, Cigar> T;
        std::vector<T> v(4);
        v[0].pos.x()=5;
        v[1].pos.y()=10;
        v[2].pos.z()=2;
        auto m = asEigenMatrix(v.begin(), v.end(), &T::pos);

        CHECK( m.cols()==3 );
        CHECK( m.rows()==4 );
        CHECK( m.row(0).x() == 5 );
        CHECK( m.row(1).y() == 10 );
        CHECK( m.row(2).z() == 2 );
        CHECK( m.sum() == 17);
        m.row(0).z()+=0.5;
        CHECK( v[0].pos.z() == Approx(0.5) );

        v[2].charge = 2;
        v[3].charge = -12;
        auto m2 = asEigenVector(v.begin()+1, v.end(), &T::charge);
        CHECK( m2.cols()==1 );
        CHECK( m2.rows()==3 );
        CHECK( m2.col(0).sum() == Approx(-10) );
    }
#endif

    template<class Trange>
        auto findName(Trange &rng, const std::string &name) {
            return std::find_if( rng.begin(), rng.end(), [&name](auto &i){ return i.name==name; });
        } //!< Returns iterator to first element with member `name` matching input

    template<class Trange>
        std::vector<int> names2ids(Trange &rng, const std::vector<std::string> &names) {
            std::vector<int> index;
            index.reserve(names.size());
            for (auto &n : names) {
                auto it = findName(rng, n);
                if (it!=rng.end())
                    index.push_back(it->id());
                else
                    throw std::runtime_error("name '" + n + "' not found");
            }
            return index;
        } //!< Convert vector of names into vector in id's from Trange (exception if not found)

    /**
     * @brief General properties for atoms
     */
    template<class T /** Particle type */>
        struct AtomData {
            public:
                T p;               //!< Particle with generic properties
                std::string name;  //!< Name
                double eps=0;      //!< LJ epsilon [kJ/mol] (pair potentials should convert to kT)
                double activity=0; //!< Chemical activity [mol/l]
                double alphax=0;   //!< Excess polarisability (unit-less)
                double dp=0;       //!< Translational displacement parameter [angstrom]
                double dprot=0;    //!< Rotational displacement parameter [degrees]
                double mw=1;       //!< Weight
                double sigma=0;    //!< Diamater for e.g Lennard-Jones etc. [angstrom]
                double tension=0;  //!< Surface tension [kT/Å^2]
                double tfe=0;      //!< Transfer free energy [J/mol/angstrom^2/M]
                int& id() { return p.id; } //!< Type id
                const int& id() const { return p.id; } //!< Type id
        };

    template<class T>
        void to_json(json& j, const AtomData<T> &a) {
            auto& _j = j[a.name];
            _j = a.p;
            _j["activity"] = a.activity / 1.0_molar;
            _j["alphax"] = a.alphax;
            _j["dp"] = a.dp / 1.0_angstrom;
            _j["dprot"] = a.dprot / 1.0_rad;
            _j["eps"] = a.eps / 1.0_kJmol;
            _j["mw"] = a.mw;
            _j["sigma"] = a.sigma / 1.0_angstrom;
            _j["tension"] = a.tension * 1.0_angstrom*1.0_angstrom / 1.0_kJmol;
            _j["tfe"] = a.tfe * 1.0_angstrom*1.0_angstrom * 1.0_molar / 1.0_kJmol;
        }

    template<class T>
        void from_json(const json& j, AtomData<T>& a) {
            if (j.is_object()==false || j.size()!=1)
                throw std::runtime_error("Invalid JSON data for AtomData");
            for (auto it=j.begin(); it!=j.end(); ++it) {
                a.name = it.key();
                auto& val = it.value();
                a.p = val;
                a.activity = val.value("activity", a.activity) * 1.0_molar;
                a.alphax   = val.value("alphax", a.alphax);
                a.dp       = val.value("dp", a.dp) * 1.0_angstrom;
                a.dprot    = val.value("dprot", a.dprot) * 1.0_rad;
                a.eps      = val.value("eps", a.eps) * 1.0_kJmol;
                a.mw       = val.value("mw", a.mw);
                a.sigma    = val.value("sigma", 0.0) * 1.0_angstrom;
                a.sigma    = 2*val.value("r", 0.5*a.sigma) * 1.0_angstrom;
                a.tension  = val.value("tension", a.tension) * 1.0_kJmol / (1.0_angstrom*1.0_angstrom);
                a.tfe      = val.value("tfe", a.tfe) * 1.0_kJmol / (1.0_angstrom*1.0_angstrom*1.0_molar);
            }
        }

    /**
     * @brief Construct vector of atoms from json array
     *
     * Items are added to existing items while if an item
     * already exists, it will be overwritten.
     */
    template<class T>
        void from_json(const json& j, std::vector<AtomData<T>> &v) {
            auto it = j.find("atomlist");
            json _j = ( it==j.end() ) ? j : *it;
            v.reserve( v.size() + _j.size() );

            for ( auto &i : _j ) {
                if ( i.is_string() ) // treat ax external file to load
                    from_json( openjson(i.get<std::string>()), v );
                else if ( i.is_object() ) {
                    AtomData<T> a = i;
                    auto it = findName( v, a.name );
                    if ( it==v.end() ) 
                        v.push_back( a ); // add new atom
                    else
                        *it = a;
                }
            }
            for (size_t i=0; i<v.size(); i++)
                v[i].id() = i; // id must match position in vector
        }

    template<typename Tparticle>
        static std::vector<AtomData<Tparticle>> atoms = {}; //!< Global instance of atom list

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] AtomData") {
        using doctest::Approx;

        json j = R"({ "atomlist" : [
             { "A": { "r":1.1 } },
             { "B": { "activity":0.2, "eps":0.05, "dp":9.8, "dprot":3.14, "mw":1.1, "tfe":0.98, "tension":0.023 } }
             ]})"_json;

        typedef Particle<Radius, Charge, Dipole, Cigar> T;

        atoms<T> = j["atomlist"].get<decltype(atoms<T>)>();
        auto &v = atoms<T>; // alias to global atom list

        CHECK(v.size()==2);
        CHECK(v.front().id()==0);
        CHECK(v.front().name=="A"); // alphabetic order in std::map
        CHECK(v.front().sigma==Approx(2*1.1e-10_m));

        AtomData<T> a = json(v.back()); // AtomData -> JSON -> AtomData

        CHECK(a.name=="B");
        CHECK(a.id()==1);
        CHECK(a.id()==a.p.id);

        CHECK(a.activity==Approx(0.2_molar));
        CHECK(a.eps==Approx(0.05_kJmol));
        CHECK(a.dp==Approx(9.8));
        CHECK(a.dprot==Approx(3.14));
        CHECK(a.mw==Approx(1.1));
        CHECK(a.tfe==Approx(0.98_kJmol/1.0_angstrom/1.0_angstrom/1.0_molar));
        CHECK(a.tension==Approx(0.023_kJmol/1.0_angstrom/1.0_angstrom));

        auto it = findName(v, "B");
        CHECK( it->id() == 1 );
        it = findName(v, "unknown atom");
        CHECK( it==v.end() );
    }
#endif

    struct json_io_base {
        std::string name;
        std::string cite;

        inline void _to_json(json&) const {};
        inline void _from_json(const json&) {};
    };
}//end of faunus namespace
