#include "forcemove.h"
#include "units.h"
#include "externalpotential.h"

namespace Faunus {
namespace Move {

using doctest::Approx;

TEST_SUITE_BEGIN("ForceMove");

TEST_CASE("[Faunus] Integrator") {
    class DummyEnergy : public Energy::Energybase {
        double energy(Change &) override {return 0.0; }
    };
    Space spc;
    DummyEnergy energy;

    auto ld = LangevinVelocityVerlet(spc, energy, R"({"friction": 8.0, "time_step": 2.0 })"_json);
    json j_ld = ld;
    CHECK_EQ(j_ld["friction_coefficient"], 8.0);
    CHECK_EQ(j_ld["time_step"], 2.0);
}

TEST_CASE("[Faunus] LangevinDynamics") {
    class DummyEnergy : public Energy::Energybase {
        double energy(Change &) override {return 0.0; }
    };
    Space spc;
    DummyEnergy energy;

    SUBCASE("[Faunus] JSON init") {
      json j_in =  R"({"nsteps": 100, "integrator": {"time_step": 0.001, "friction": 2.0}})"_json;
      LangevinDynamics ld(spc, energy);
      ld.from_json(j_in);
      json j_out = ld;

      //CHECK_EQ(j_out, j_in);
    }
}
TEST_SUITE_END();
} // namespace Move
} // namespace Faunus
