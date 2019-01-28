#pragma once
#include <iomanip>
#include <limits>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <type_traits>
#include <string>
#include <cmath>
#include <random>
#include <memory>
#include <chrono>
#include <Eigen/Geometry>
#include <nlohmann/json.hpp>

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <range/v3/view.hpp>
#endif
 
#ifdef DOCTEST_LIBRARY_INCLUDED
#include "units.h"
#endif
#include "random.h"

// Eigen<->JSON (de)serialization
namespace Eigen {

    template<typename T>
        void to_json(nlohmann::json& j, const T &p) {
            auto d = p.data();
            for (int i=0; i<(int)p.size(); ++i)
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
}

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

    using namespace std::string_literals;

    template<class T1, class T2>
        int distance(T1 first, T2 last) {
            return std::distance( &(*first), &(*last) );
        } //!< Distance between two arbitrary contiguous iterators

    template<class T>
        int size(T &rng) {
            return std::distance(rng.begin(), rng.end());
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

    json merge( const json &a, const json &b ); //!< Merge two json objects
    json openjson( const std::string &file, bool=true); //!< Read json file into json object (w. syntax check)

    /**
     * @brief Check for unknown keys in JSON object
     * @param j JSON object to check
     * @param okkeys Valid keys
     * @param exception If true a runtime error will be thrown if unknown key is found
     */
    bool assertKeys(const json&, const std::vector<std::string>&, bool=true);

    struct xjson : private json {
        xjson(const json &j);
        bool empty() const;
        size_type count(const std::string &key) const;
        inline auto dump(int w=-1) const { return json::dump(w); }

        void clear();
        json at(const std::string &key);
        json operator[](const std::string &key);
        void erase(const std::string &key);

        template<class T> T value(const std::string &key, const T &fallback) {
            return (count(key)>0) ? at(key).get<T>() : fallback;
        }
    }; //!< Like json, but delete entries after access

    double _round(double x, int n=3); //!< Round to n number of significant digits
    void _roundjson(json &j, int n=3); // round float objects to n number of significant digits
    double value_inf(const json &j, const std::string &key); //!< Extract floating point from json and allow for 'inf' and '-inf'

    /**
     * @brief Class for showing help based on input errors
     *
     * If no valid database files are found using `load()`,
     * or of the key is not found, an empty string is
     * returned by the call operator. The idea is that this
     * functionality is completely optional.
     */
    class TipFromTheManual {
        private:
            json db; // database
            Random slump;
        public:
            bool tip_already_given=false;
            void load(const std::vector<std::string>&);
            std::string operator[](const std::string&);
    };

    extern TipFromTheManual usageTip; // global instance

    struct Tensor : public Eigen::Matrix3d {
        typedef Eigen::Matrix3d base;

        Tensor() {
            base::setZero();
        } //!< Constructor, clear data

        Tensor( double xx, double xy, double xz, double yy, double yz, double zz ); //!< Construct from input
        void rotate( const base &m ); //!< Rotate using rotation matrix. Remove?
        void eye();

        template<typename T>
            Tensor( const Eigen::MatrixBase<T> &other ) : base(other) {}

        template<typename T>
            Tensor &operator=( const Eigen::MatrixBase<T> &other )
            {
                base::operator=(other);
                return *this;
            }
    }; //!< Tensor class

    void to_json(nlohmann::json& j, const Tensor &t); //!< Tensor -> Json
    void from_json(const nlohmann::json& j, Tensor &t); //!< Json -> Tensor

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
     * @brief Quaternion rotation routine using the Eigen library
     */
    struct QuaternionRotate : public std::pair<Eigen::Quaterniond, Eigen::Matrix3d> {

        typedef std::pair<Eigen::Quaterniond, Eigen::Matrix3d> base;
        //using base::pair;
        using base::first;
        using base::second;

        double angle=0; //!< Rotation angle

        QuaternionRotate();
        QuaternionRotate(double angle, Point u);
        void set(double angle, Point u);
        Point operator()( Point a, std::function<void(Point&)> boundary = [](Point&){}, const Point &shift={0,0,0} ) const; //!< Rotate point w. optional PBC boundaries
        auto operator()( const Eigen::Matrix3d &a ) const; //!< Rotate matrix/tensor
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

    /**
     * @brief Generate reference view to data members in container
     *
     *     struct data { int N=0; };
     *     std::vector<data> vec(2);               // original vector
     *     auto Nvec = member_view(vec, &data::N); // ref. to all N's in vec
     *     for (int &i : Nvec)                     // modify N values
     *        i=1;
     *     assert( vec[0].N == vec[1].N == 1 );    // original vector was changed
     */
    template<typename Container, typename Value, typename Member>
        auto member_view(Container &p, Value Member::*m) {
            std::vector<std::reference_wrapper<Value>> v;
            v.reserve( p.size() );
            std::transform(p.begin(), p.end(), std::back_inserter(v), [&](auto &i) -> Value& {return i.*m;});
            return v;
        }

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

#ifdef _DOCTEST_LIBRARY_INCLUDED
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

    /** Add growing suffix filename until non-existing name is found */
    std::string addGrowingSuffix(const std::string&);


}//end of faunus namespace
