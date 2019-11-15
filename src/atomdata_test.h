#include "atomdata.h"
#include "units.h"

namespace Faunus {

using doctest::Approx;

TEST_SUITE_BEGIN("AtomData");

TEST_CASE("[Faunus] AtomData") {
    using doctest::Approx;

    json j = R"({ "atomlist" : [
             { "A": { "sigma": 2.5, "pactivity":2, "eps_custom": 0.1 } },
             { "B": { "r":1.1, "activity":0.2, "eps":0.05, "dp":9.8, "dprot":3.14, "mw":1.1, "tfe":0.98, "tension":0.023 } }
             ]})"_json;

    pc::temperature = 298.15_K;
    atoms = j["atomlist"].get<decltype(atoms)>();
    auto &v = atoms; // alias to global atom list

    CHECK_EQ(v.size(), 2);
    CHECK_EQ(v.front().id(), 0);
    CHECK_EQ(v.front().name, "A");                             // alphabetic order in std::map
    CHECK(v.front().interaction.get("sigma") == Approx(2.5));      // raw number, no units
    CHECK(v.front().interaction.get("eps_custom") == Approx(0.1)); // raw number, no units

    CHECK_EQ(std::isnan(v.front().interaction.get("eps_unknown")), true);
    // CHECK_THROWS_AS_MESSAGE(v.front().interaction.get("eps_unknown"), std::runtime_error, "unknown atom property");
    CHECK(v.front().sigma == Approx(2.5e-10_m));
    CHECK(v.front().activity == Approx(0.01_molar));
    CHECK(v.back().tfe == Approx(0.98_kJmol / (1.0_angstrom * 1.0_angstrom * 1.0_molar)));

    AtomData a = json(v.back()); // AtomData -> JSON -> AtomData

    CHECK_EQ(a.name, "B");
    CHECK_EQ(a.id(), 1);
    CHECK(a.activity == Approx(0.2_molar));
    CHECK(a.interaction.get("sigma") == Approx(2.2)); // raw number, no units
    CHECK(a.interaction.get("eps") == Approx(0.05));  // raw number, no units
    CHECK(a.dp == Approx(9.8));
    CHECK(a.dprot == Approx(3.14));
    CHECK(a.mw == Approx(1.1));
    CHECK(a.tfe == Approx(0.98_kJmol / 1.0_angstrom / 1.0_angstrom / 1.0_molar));
    CHECK(a.tension == Approx(0.023_kJmol / 1.0_angstrom / 1.0_angstrom));

    auto it = findName(v, "B");
    CHECK_EQ(it->id(), 1);
    it = findName(v, "unknown atom");
    CHECK_EQ(it, v.end());
}
TEST_SUITE_END();
} // namespace Faunus
