#include "random.h"
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream>

namespace Faunus {

void from_json(const nlohmann::json &j, Random &random) {
    if (j.is_object()) {
        auto seed = j.value("seed", std::string());
        try {
            if (seed == "default" or seed == "fixed") { // use default seed, i.e. do nothing
                return;
            } else if (seed == "hardware") { // use hardware seed
                random.engine = decltype(random.engine)(std::random_device()());
            } else if (!seed.empty()) { // read engine state
                std::stringstream stream(seed);
                stream.exceptions(std::ios::badbit | std::ios::failbit);
                stream >> random.engine;
            }
        } catch (std::exception &e) {
            std::cerr << "could not initialize random engine - falling back to fixed seed." << std::endl;
        }
    }
}

void to_json(nlohmann::json &j, const Random &random) {
    std::ostringstream stream;
    stream << random.engine; // dump engine state to stream
    j["seed"] = stream.str();
}

void Random::seed() { engine = std::mt19937(std::random_device()()); }

Random::Random() : dist01(0, 1) {}

double Random::operator()() { return dist01(engine); }

Random random; // Global instance
} // namespace Faunus
