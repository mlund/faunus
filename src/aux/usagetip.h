#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <nlohmann/json.hpp>

namespace Faunus {

class Random;

/**
 * @brief Class for showing help based on input errors
 *
 * If no valid database files are found using `load()`,
 * or of the key is not found, an empty string is
 * returned by the call operator. The idea is that this
 * functionality is completely optional.
 */
class TipFromTheManual
{
  private:
    nlohmann::json database; // database
    std::unique_ptr<Random> random;
    bool tip_already_given = false;

  public:
    std::string output_buffer; // accumulate output here
    bool quiet = true;         // if true, operator[] returns empty string
    bool asciiart = true;
    TipFromTheManual();
    void load(const std::vector<std::string>&);
    std::string operator[](std::string_view key);
    void pick(const std::string&);
};

extern TipFromTheManual usageTip; // global instance

} // namespace Faunus
