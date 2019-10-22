#include "bonds.h"
#include "units.h"

namespace Faunus {
namespace Potential {

using doctest::Approx;

TEST_SUITE_BEGIN("Bonds");
TEST_CASE("[Faunus] BondData") {
    std::shared_ptr<BondData> b;

    // exact match required
    CHECK_THROWS(b = R"({ "harmoNIC": {"index":[2,3], "k":0.5, "req":2.1}} )"_json;);

    // test harmonic
    SUBCASE("HarmonicBond") {
        json j = R"({ "harmonic": {"index":[2,3], "k":0.5, "req":2.1}} )"_json;
        b = j;
        CHECK(j == json(b));
        CHECK_THROWS(b = R"({"harmonic": { "index":[2], "k":0.5, "req":2.1}} )"_json);
        CHECK_THROWS(b = R"({"harmonic": { "index":[2,3], "req":2.1}} )"_json);
        CHECK_THROWS(b = R"({"harmonic": { "index":[2,3], "k":2.1}} )"_json);
    }

    // test fene
    SUBCASE("FENEBond") {
        json j = R"({"fene": { "index":[2,3], "k":1, "rmax":2.1 }} )"_json;
        b = j;
        CHECK(j == json(b));
        CHECK_THROWS(b = R"({"fene": { "index":[2,3,4], "k":1, "rmax":2.1}} )"_json);
        CHECK_THROWS(b = R"({"fene": { "index":[2,3], "rmax":2.1}} )"_json);
        CHECK_THROWS(b = R"({"fene": { "index":[2,3], "k":1}} )"_json);
    }

    // test fene+wca
    SUBCASE("FENEWCABond") {
        json j = R"({"fene+wca": { "index":[2,3], "k":1, "rmax":2.1, "eps":2.48, "sigma":2}} )"_json;
        b = j;
        CHECK(j == json(b));
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3,4], "k":1, "rmax":2.1, "eps":2.48, "sigma":2}} )"_json);
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3], "rmax":2.1, "eps":2.48, "sigma":2}} )"_json);
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3], "k":1, "eps":2.48, "sigma":2}} )"_json);
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3], "k":1, "rmax":2.1, "eps":2.48}} )"_json);
        CHECK_THROWS(b = R"({"fene+wca": { "index":[2,3], "k":1, "rmax":2.1, "sigma":2}} )"_json);
    }

    // test harmonic
    SUBCASE("HarmonicTorsion") {
        json j = R"({ "harmonic_torsion": {"index":[0,1,2], "k":0.5, "aeq":60}} )"_json;
        b = j;
        CHECK(j == json(b));
        CHECK_THROWS(b = R"({"harmonic_torsion": { "index":[2], "k":0.5, "aeq":2.1}} )"_json);
        CHECK_THROWS(b = R"({"harmonic_torsion": { "index":[0,1,2], "aeq":2.1}} )"_json);
        CHECK_THROWS(b = R"({"harmonic_torsion": { "index":[0,1,3], "k":2.1}} )"_json);
    }

    // test bond filter
    SUBCASE("filterBonds()") {
        std::vector<std::shared_ptr<BondData>> bonds = {
            R"({"fene":      {"index":[2,3], "k":1, "rmax":2.1, "eps":2.48}} )"_json,
            R"({"harmonic" : {"index":[2,3], "k":0.5, "req":2.1} } )"_json};
        auto filt = filterBonds(bonds, BondData::HARMONIC);
        CHECK(filt.size() == 1);
        CHECK(filt[0]->type() == BondData::HARMONIC);
        CHECK(filt[0] == bonds[1]); // filt should contain references to bonds
    }
}
TEST_SUITE_END();
} // namespace Potential
} // namespace Faunus
