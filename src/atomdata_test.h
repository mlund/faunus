#include "atomdata.h"
#include "units.h"

namespace Faunus {

using doctest::Approx;

TEST_SUITE_BEGIN("AtomData");

TEST_CASE("[Faunus] AtomData") {
    using doctest::Approx;

    json j = R"({ "atomlist" : [
             { "A": { "r":1.1, "pactivity":2 } },
             { "B": { "activity":0.2, "eps":0.05, "dp":9.8, "dprot":3.14, "mw":1.1, "tfe":0.98, "tension":0.023 } }
             ]})"_json;

    pc::temperature = 298.15;

    atoms = j["atomlist"].get<decltype(atoms)>();
    auto &v = atoms; // alias to global atom list

    CHECK(v.size() == 2);
    CHECK(v.front().id() == 0);
    CHECK(v.front().name == "A"); // alphabetic order in std::map
    CHECK(v.front().sigma == Approx(2 * 1.1e-10_m));
    CHECK(v.front().activity == Approx(0.01_molar));
    CHECK(v.back().tfe == Approx(0.98_kJmol / (1.0_angstrom * 1.0_angstrom * 1.0_molar)));

    AtomData a = json(v.back()); // AtomData -> JSON -> AtomData

    CHECK(a.name == "B");
    CHECK(a.id() == 1);
    CHECK(a.activity == Approx(0.2_molar));
    CHECK(a.eps == Approx(0.05_kJmol));
    CHECK(a.dp == Approx(9.8));
    CHECK(a.dprot == Approx(3.14));
    CHECK(a.mw == Approx(1.1));
    CHECK(a.tfe == Approx(0.98_kJmol / 1.0_angstrom / 1.0_angstrom / 1.0_molar));
    CHECK(a.tension == Approx(0.023_kJmol / 1.0_angstrom / 1.0_angstrom));

    auto it = findName(v, "B");
    CHECK(it->id() == 1);
    it = findName(v, "unknown atom");
    CHECK(it == v.end());
}
TEST_SUITE_END();
} // namespace Faunus
