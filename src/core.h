#pragma once
#include <vector>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <range/v3/distance.hpp>

// forward declare logger
namespace spdlog {
    class logger;
}

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
    struct Random;

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
            return ranges::distance(rng.begin(), rng.end());
        } //!< Size of arbitrary range

    json merge( const json &a, const json &b ); //!< Merge two json objects
    json openjson( const std::string &file, bool=true); //!< Read json file into json object (w. syntax check)

    /**
     * @brief Check for unknown keys in JSON object
     * @param j JSON object to check
     * @param okkeys Valid keys
     * @param exception If true a runtime error will be thrown if unknown key is found
     */
    bool assertKeys(const json&, const std::vector<std::string>&, bool=true);

    /**
     * @brief Like json, but delete entries after access
     *
     * Only selected functions from the json class is exposed and
     * and accessing a key it will be deleted afterwards, i.e. a key
     * can be used only once. This can be used to check for unknown
     * keys as the object should be zero after being processed by
     * i.e. `from_json` or similar.
     */
    struct SingleUseJSON : private json {
        SingleUseJSON(const json &j);
        bool empty() const;
        size_type count(const std::string &key) const;
        inline auto dump(int w=-1) const { return json::dump(w); }
        inline bool is_object() const { return json::is_object(); }

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
        public:
            bool tip_already_given=false;
            void load(const std::vector<std::string>&);
            std::string operator[](const std::string&);
    };

    extern TipFromTheManual usageTip; // global instance

    extern std::shared_ptr<spdlog::logger> faunus_logger; // global instance
    extern std::shared_ptr<spdlog::logger> mcloop_logger; // global instance

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

    /**
     * @brief Convert cartesian- to cylindrical-coordinates
     * @note Input (x,y,z), output \f$ (r,\theta, h) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [-\pi,\pi) \f$,
     * and \f$ h\in (-\infty,\infty] \f$.
     */
    Point xyz2rth(const Point &, const Point &origin = {0, 0, 0}, const Point &dir = {0, 0, 1},
                  const Point &dir2 = {1, 0, 0});

    /**
     * @brief Convert cartesian- to spherical-coordinates
     * @note Input (x,y,z), output \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [-\pi,\pi) \f$,
     * and \f$ \phi\in [0,\pi] \f$.
     */
    Point xyz2rtp(const Point &, const Point &origin = {0, 0, 0});

    /**
     * @brief Convert spherical- to cartesian-coordinates
     * @param origin The origin to be added (optional)
     * @note Input \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [0,2\pi) \f$, and \f$ \phi\in
     * [0,\pi] \f$, and output (x,y,z).
     */
    Point rtp2xyz(const Point &rtp, const Point &origin = {0, 0, 0});

    /** Add growing suffix filename until non-existing name is found */
    std::string addGrowingSuffix(const std::string&);

    Point ranunit(Random &,
                  const Point &dir = {1, 1, 1}); //!< Random unit vector using Neuman's method ("sphere picking")
    Point ranunit_polar(Random &);               //!< Random unit vector using polar coordinates ("sphere picking")

}//end of faunus namespace
