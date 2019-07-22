#include "random.h"
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream>

namespace Faunus {

    void from_json(const nlohmann::json &j, Random &r) {
        if (j.is_object()) {
            auto seed = j.value("seed", std::string());
            try {
                if (seed=="default" or seed=="fixed")
                    return;
                else if (seed=="hardware")
                    r.engine = decltype(r.engine)(std::random_device()());
                else if (!seed.empty()) {
                    std::stringstream s(seed);
                    s.exceptions( std::ios::badbit | std::ios::failbit );
                    s >> r.engine;
                }
            }
            catch (std::exception &e) {
                std::cerr << "could not initialize random - falling back to fixed seed." << std::endl;
            }
        }
    }

    void to_json(nlohmann::json &j, const Random &r) {
        std::ostringstream o;
        o << r.engine;
        j["seed"] = o.str();
    }

    void Random::seed() { engine = std::mt19937(std::random_device()()); }

    Random::Random() : dist01(0,1) {}

    double Random::operator()() { return dist01(engine); }

    int Random::range(int min, int max) {
        std::uniform_int_distribution<int> d(min, max);
        return d(engine);
    }

    Random random; // Global instance
}
