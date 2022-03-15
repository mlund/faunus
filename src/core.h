#pragma once

#include <vector>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <spdlog/fmt/fmt.h>
#include <range/v3/range/concepts.hpp>

// forward declare logger
namespace spdlog {
class logger;
}

extern template class nlohmann::basic_json<>;

/** @brief Faunus main namespace */
namespace Faunus {

using Point = Eigen::Vector3d;          //!< 3D vector used for positions, velocities, forces etc.
using PointVector = std::vector<Point>; //!< Vector of 3D vectors
using json = nlohmann::json;            //!< JSON object
class Random;

/** Concept for a range of points */
template <class T>
concept RequirePoints = ranges::cpp20::range<T> && std::is_same_v<ranges::cpp20::range_value_t<T>, Point>;

/** Concept for an iterator to a `Point` */
template <class T>
concept RequirePointIterator = std::is_convertible_v<ranges::cpp20::iter_value_t<T>, Point>;

using namespace std::string_literals;

json loadJSON(const std::string& filename); //!< Read json filename into json object (w. syntax check)

/**
 * @brief Like json, but delete entries after access
 *
 * Only selected functions from the json class is exposed and
 * and accessing a key it will be deleted afterwards, i.e. a key
 * can be used only once. This can be used to check for unknown
 * keys as the object should be zero after being processed by
 * i.e. `from_json` or similar.
 *
 * @todo This class should be retired and handled by JSON schema instead
 */
struct SingleUseJSON : public json {
    SingleUseJSON(const json&);
    bool empty() const;
    size_type count(const std::string&) const;
    std::string dump(int = -1) const;
    bool is_object() const;

    void clear();
    json at(const std::string&);
    json operator[](const std::string&);
    void erase(const std::string&);

    template <class T> T value(const std::string& key, const T& fallback) {
        return (count(key) > 0) ? at(key).get<T>() : fallback;
    }
};

double roundValue(double value, int number_of_digits = 3); //!< Round to n number of significant digits
void roundJSON(json& j, int number_of_digits = 3);         //!< Round float objects to n number of significant digits
double getValueInfinity(const json& j,
                        const std::string& key); //!< Extract floating point from json and allow for 'inf' and '-inf'

/**
 * @brief Returns a key-value pair from a JSON object which contains a single key.
 *
 * JSON objects having a single key-value pair are a common pattern in JSON configuration used in Faunus. This
 * function provides a convenient way to handle it.
 *
 * @param j JSON object
 * @return tuple [key as a string, value as a JSON]
 * @throw std::runtime_error  when not a JSON object or the object is empty or the object contains more than a
 *   single value
 */
std::tuple<const std::string&, const json&> jsonSingleItem(const json& j);

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
    json database; // database
    std::shared_ptr<Random> random;
    bool tip_already_given = false;

  public:
    std::string buffer; // accumulate output here
    bool quiet = true;  // if true, operator[] returns empty string
    bool asciiart = true;
    TipFromTheManual();
    void load(const std::vector<std::string>&);
    std::string operator[](std::string_view key);
    void pick(const std::string&);
};

extern TipFromTheManual usageTip;                     // global instance
extern std::shared_ptr<spdlog::logger> faunus_logger; // global instance
extern std::shared_ptr<spdlog::logger> mcloop_logger; // global instance

/**
 * @brief Eigen::Map facade to data members in STL container
 *
 * No data is copied and modifications of the Eigen object
 * modifies the original container and vice versa.
 *
 * Example:
 *
 *    std::vector<Tparticle> v(10);
 *    auto m1 = asEigenVector(v.begin, v.end(), &Tparticle::pos);    --> 10x3 maxtrix view
 *    auto m2 = asEigenMatrix(v.begin, v.end(), &Tparticle::charge); --> 10x1 vector view
 *
 * @warning Be careful that objects are properly aligned and divisible with `sizeof<double>`
 */
template <typename dbl = double, class iter, class memberptr> auto asEigenMatrix(iter begin, iter end, memberptr m) {
    using T = typename std::iterator_traits<iter>::value_type;
    static_assert(sizeof(T) % sizeof(dbl) == 0, "value_type size must multiples of double");
    constexpr size_t s = sizeof(T) / sizeof(dbl);
    constexpr size_t cols = sizeof((static_cast<T*>(0))->*m) / sizeof(dbl);
    using Tmatrix = Eigen::Matrix<dbl, Eigen::Dynamic, cols>;
    return Eigen::Map<Tmatrix, 0, Eigen::Stride<1, s>>((dbl*)&(*begin.*m), end - begin, cols).array();
}

template <typename dbl = double, class iter, class memberptr> auto asEigenVector(iter begin, iter end, memberptr m) {
    using T = typename std::iterator_traits<iter>::value_type;
    static_assert(std::is_same<dbl&, decltype((static_cast<T*>(0))->*m)>::value, "member must be a scalar");
    return asEigenMatrix<dbl>(begin, end, m).col(0);
}

/**
 * @brief Convert cartesian- to cylindrical-coordinates
 * @note Input (x,y,z), output \f$ (r,\theta, h) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [-\pi,\pi) \f$,
 * and \f$ h\in (-\infty,\infty] \f$.
 */
Point xyz2rth(const Point&, const Point& origin = {0, 0, 0}, const Point& dir = {0, 0, 1},
              const Point& dir2 = {1, 0, 0});

/**
 * @brief Convert cartesian- to spherical-coordinates
 * @note Input (x,y,z), output \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [-\pi,\pi) \f$,
 * and \f$ \phi\in [0,\pi] \f$.
 */
Point xyz2rtp(const Point&, const Point& origin = {0, 0, 0});

/**
 * @brief Convert spherical- to cartesian-coordinates
 * @param origin The origin to be added (optional)
 * @note Input \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [0,2\pi) \f$, and \f$ \phi\in
 * [0,\pi] \f$, and output (x,y,z).
 */
Point rtp2xyz(const Point& rtp, const Point& origin = {0, 0, 0});

std::string addGrowingSuffix(const std::string&); //!< Add growing suffix filename until non-existing name is found

Point randomUnitVector(
    Random& rand,
    const Point& directions = Point::Ones()); //!< Random unit vector using Neuman's method ("sphere picking")

Point randomUnitVectorPolar(Random& rand); //!< Random unit vector using polar coordinates ("sphere picking")

//! Common ancestor of Faunus specific runtime errors
struct GenericError : public std::runtime_error {
    explicit GenericError(const std::exception& e);
    explicit GenericError(const std::runtime_error& e);
    explicit GenericError(const std::string& msg);
    explicit GenericError(const char* msg);
    template <class... Args>
    explicit GenericError(std::string_view fmt, const Args&... args)
        : std::runtime_error(fmt::vformat(fmt, fmt::make_format_args(args...))) {}
};

//! Exception to be thrown when parsing json configuration
struct ConfigurationError : public GenericError {
    using GenericError::GenericError;
    const json& attachedJson() const;
    ConfigurationError& attachJson(const json& j);

  private:
    json attached_json;
};

//! Exception to be thrown on IO errors
struct IOError : public GenericError {
    using GenericError::GenericError;
};

/**
 * @brief Nicely displays nested exceptions using a logger
 *
 * @param logger  logger to display the exception with
 * @param e  exception to display (possibly nested)
 * @param level  internal counter for recursion
 */
void displayError(spdlog::logger& logger, const std::exception& e, int level = 0);
} // namespace Faunus
