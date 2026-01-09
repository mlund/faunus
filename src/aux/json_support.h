#pragma once

#include <optional>
#include <string>
#include <tuple>
#include <nlohmann/json.hpp>

namespace Faunus {

using json = nlohmann::json;

json loadJSON(
    const std::string& filename); //!< Read json filename into json object (w. syntax check)

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
struct SingleUseJSON : public json
{
    explicit SingleUseJSON(const json&);
    [[nodiscard]] bool empty() const;
    [[nodiscard]] size_type count(const std::string&) const;
    [[nodiscard]] std::string dump(int = -1) const;
    [[nodiscard]] bool is_object() const;

    void clear();
    json at(const std::string&);
    json operator[](const std::string&);
    void erase(const std::string&);

    template <class T> T value(const std::string& key, const T& fallback)
    {
        return (count(key) > 0) ? at(key).get<T>() : fallback;
    }
};

double roundValue(double value,
                  int number_of_digits = 3); //!< Round to n number of significant digits
void roundJSON(json& j,
               int number_of_digits = 3); //!< Round float objects to n number of significant digits
double getValueInfinity(
    const json& j,
    const std::string& key); //!< Extract floating point from json and allow for 'inf' and '-inf'

/**
 * @brief Returns a key-value pair from a JSON object which contains a single key.
 *
 * JSON objects having a single key-value pair are a common pattern in JSON configuration used in
 * Faunus. This function provides a convenient way to handle it.
 *
 * @param j JSON object
 * @return tuple [key as a string, value as a JSON]
 * @throw std::runtime_error  when not a JSON object or the object is empty or the object contains
 * more than a single value
 */
std::tuple<const std::string&, const json&> jsonSingleItem(const json& j);

/// Extract JSON value associated with `key` into `std::optional`
///
/// If the value does not exist, return `std::nullopt`
///
/// @throws If value exists but cannot be extracted as `T`.
template <typename T> std::optional<T> get_optional(const json& j, std::string_view key)
{
    if (const auto it = j.find(key); it != j.end()) {
        return it->get<T>(); // may throw exception
    }
    return std::nullopt;
}

} // namespace Faunus
