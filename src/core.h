#pragma once

#include <format>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <ranges>
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

using Point = Eigen::Vector3d;          //!< 3D vector used for positions, velocities, forces etc.
using PointVector = std::vector<Point>; //!< Vector of 3D vectors
using json = nlohmann::json;            //!< JSON object
class Random;

namespace Geometry {
//! Function to apply PBC to a position
using BoundaryFunction = std::function<void(Point&)>;
//! Function to calculate the (minimum) distance between two points
using DistanceFunction = std::function<Point(const Point&, const Point&)>;
} // namespace Geometry

/** Concept for a range of points */
template <class T>
concept RequirePoints =
    std::ranges::range<T> && std::is_same_v<std::ranges::range_value_t<T>, Point>;

/** Concept for an iterator to a `Point` */
template <class T>
concept RequirePointIterator = std::is_convertible_v<std::iter_value_t<T>, Point>;

using namespace std::string_literals;

extern std::shared_ptr<spdlog::logger> faunus_logger; // global instance
extern std::shared_ptr<spdlog::logger> mcloop_logger; // global instance

//! Common ancestor of Faunus specific runtime errors
struct GenericError : public std::runtime_error
{
    explicit GenericError(const std::exception& e);
    explicit GenericError(const std::runtime_error& e);
    explicit GenericError(const std::string& msg);
    explicit GenericError(const char* msg);

    template <class... Args>
    explicit GenericError(std::format_string<Args...> fmt, Args&&... args)
        : std::runtime_error(std::format(fmt, std::forward<Args>(args)...))
    {
    }
};

//! Exception to be thrown when parsing json configuration
struct ConfigurationError : public GenericError
{
    using GenericError::GenericError;
    [[nodiscard]] const json& attachedJson() const;
    ConfigurationError& attachJson(const json& j);

  private:
    json attached_json;
};

//! Exception to be thrown on IO errors
struct IOError : public GenericError
{
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

// Include moved items for backwards compatibility
#include "aux/json_support.h"
#include "aux/eigensupport.h"
#include "aux/coordinates.h"
#include "aux/usagetip.h"
#include "electrolyte.h"
