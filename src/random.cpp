#include <doctest/doctest.h>
#include "random.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <string>
#include <sstream>

namespace Faunus {

void from_json(const nlohmann::json& j, Random& rng)
{
    if (j.is_object()) {
        auto seed = j.value("seed", std::string());
        try {
            if (seed == "default" or seed == "fixed") { // use default seed, i.e. do nothing
                return;
            }
            else if (seed == "hardware") { // use hardware seed
                rng.engine = decltype(rng.engine)(std::random_device()());
            }
            else if (!seed.empty()) { // read engine state
                std::stringstream stream(seed);
                stream.exceptions(std::ios::badbit | std::ios::failbit);
                stream >> rng.engine;
            }
        }
        catch (std::exception& e) {
            std::cerr << "could not initialize rng engine - falling back to fixed seed."
                      << std::endl;
        }
    }
}

void to_json(nlohmann::json& j, const Random& random)
{
    std::ostringstream stream;
    stream << random.engine; // dump engine state to stream
    j["seed"] = stream.str();
    if constexpr (std::is_same<RandomNumberEngine, std::mt19937>::value) {
        j["engine"] = "Mersenne Twister (std::mt19937)";
    }
    else {
        j["engine"] = "Permuted Congruential Generator (pcg32)";
    }
}

void Random::seed()
{
    engine = RandomNumberEngine(std::random_device()());
}

Random::Random()
    : dist01(0, 1)
{
}

double Random::operator()()
{
    return dist01(engine);
}

Random random; // Global instance
} // namespace Faunus

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("[Faunus] Random")
{
    using namespace Faunus;
    Random slump, slump2; // local instances

    CHECK_EQ(slump(), slump2()); // deterministic initialization by default; the global random
                                 // variable cannot be used for comparison as its state is not reset
                                 // at the beginning of each test case

    int min = 10, max = 0, N = 1e6;
    double x = 0;
    for (int i = 0; i < N; i++) {
        int j = slump.range(0, 9);
        if (j < min)
            min = j;
        if (j > max)
            max = j;
        x += j;
    }
    CHECK_EQ(min, 0);
    CHECK_EQ(max, 9);
    CHECK_EQ(std::fabs(x / N), doctest::Approx(4.5).epsilon(0.01));

    Random r1 = R"( {"seed" : "hardware"} )"_json; // non-deterministic seed
    Random r2;                                     // default is a deterministic seed
    CHECK((r1() != r2()));
    Random r3 = nlohmann::json(r1); // r1 --> json --> r3
    CHECK_EQ(r1(), r3());

    // check if random_device works
    Random a, b;
    CHECK_EQ(a(), b());
    a.seed();
    b.seed();
    CHECK((a() != b()));
}

TEST_CASE("[Faunus] WeightedDistribution")
{
    using namespace Faunus;
    WeightedDistribution<double> v;

    v.push_back(0.5);
    CHECK_EQ(v.getLastIndex(), 0);
    CHECK_EQ(v.size(), 1);

    v.push_back(0.1, 4);
    CHECK_EQ(v.getLastIndex(), 1);
    CHECK_EQ(v.size(), 2);
    CHECK(not v.empty());

    int N = 1e4;
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += v.sample(Faunus::random.engine);
    }
    CHECK_EQ(sum / N, doctest::Approx((0.5 * 1 + 0.1 * 4) / (1 + 4)).epsilon(0.05));

    std::vector<int> weights = {2, 1};
    v.setWeight(weights.begin(), weights.end());
    sum = 0;
    for (int i = 0; i < N; i++) {
        sum += v.sample(Faunus::random.engine);
    }
    CHECK_EQ(sum / N, doctest::Approx((0.5 * 2 + 0.1 * 1) / (2 + 1)).epsilon(0.05));

    weights = {2, 1, 1};
    CHECK_THROWS(v.setWeight(weights.begin(), weights.end()));

    v.clear();
    CHECK(v.empty());
}
#endif
