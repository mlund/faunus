#pragma once

#include <vector>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

// forward declare logger
namespace spdlog {
    class logger;
}

extern template class nlohmann::basic_json<>;

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
    struct SingleUseJSON : public json {
        SingleUseJSON(const json &);
        bool empty() const;
        size_type count(const std::string &) const;
        std::string dump(int=-1) const;
        bool is_object() const;

        void clear();
        json at(const std::string &);
        json operator[](const std::string &);
        void erase(const std::string &);

        template<class T> T value(const std::string &key, const T &fallback) {
            return (count(key)>0) ? at(key).get<T>() : fallback;
        }
    }; //!< Like json, but delete entries after access

    double _round(double, int n=3); //!< Round to n number of significant digits
    void _roundjson(json &, int n=3); // round float objects to n number of significant digits
    double value_inf(const json &, const std::string &); //!< Extract floating point from json and allow for 'inf' and '-inf'

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
            std::shared_ptr<Random> random;

          public:
            std::string buffer; // accumulate output here
            bool quiet = true;  // if operator[] returns empty string
            TipFromTheManual();
            bool asciiart = true;
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
     *     std::vector<data> v(2);               // original vector
     *     auto Nvec = member_view(v.begin(), v.end(), &data::N); // ref. to all N's in vec
     *     for (int &i : Nvec)                     // modify N values
     *        i=1;
     *     assert( v[0].N == v[1].N == 1 );    // original vector was changed
     */
    template<typename iter, typename T, typename Member>
        auto member_view(iter begin, iter end, T Member::*m) {
            std::vector<std::reference_wrapper<T>> v;
            v.reserve( std::distance(begin, end));
            std::transform(begin, end, std::back_inserter(v), [&](auto &i) -> T& {return i.*m;});
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
 * @brief Returns filtered *view* of iterator range
 */
template<typename iter, typename T=typename std::iterator_traits<iter>::value_type>
auto filter(iter begin, iter end, std::function<bool(T&)> unarypredicate) {
    std::vector<std::reference_wrapper<T>> v;
    for (auto& i=begin; i!=end; i++)
        if (unarypredicate(i))
            v.push_back(i);
    return v;
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
