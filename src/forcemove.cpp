#include "forcemove.h"
#include "random.h"
#include "energy.h"

namespace Faunus::Move {

TEST_SUITE_BEGIN("ForceMove");

/**
 * @brief Compute a single dimension contribution to the mean square thermal speed of a particle, i.e., compute 〈v²_x〉.
 *
 * This is equivalent to a one third of a total mean square thermal speed of a particle in three dimensions,
 * i.e., 〈v²〉 = 3 × 〈v²_x〉.
 *
 * @param mass  particle's mass in g/mol
 * @return mean square thermal speed in (Å/ps)²
 */
static inline double meanSquareSpeedComponent(T mass) {
    return (pc::kT() / mass * 1.0_kg) * ((1.0_m * 1.0_m) / (1.0_s * 1.0_s));
}

// =============== IntegratorBase  ===============

IntegratorBase::IntegratorBase(Space &spc, Energy::Energybase &energy) : spc(spc), energy(energy) {}

void from_json(const json &j, IntegratorBase &i) { i.from_json(j); }
void to_json(json &j, const IntegratorBase &i) { i.to_json(j); }

// =============== LangevinVelocityVerlet ===============

LangevinVelocityVerlet::LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy)
    : IntegratorBase(spc, energy) {}

LangevinVelocityVerlet::LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy, double time_step,
                                               double friction_coefficient)
    : IntegratorBase(spc, energy), time_step(time_step), friction_coefficient(friction_coefficient) {}

LangevinVelocityVerlet::LangevinVelocityVerlet(Space &spc, Energy::Energybase &energy, const json &j)
    : LangevinVelocityVerlet::LangevinVelocityVerlet(spc, energy) {
    from_json(j);
}

inline Point LangevinVelocityVerlet::positionIncrement(const Point& velocity) {
    return 0.5 * time_step * velocity;
}

inline Point LangevinVelocityVerlet::velocityIncrement(const Point& force, const double mass) {
    // As forces are in kT per ångström units (a hybrid between reduced energy units and absolute units), we use
    // the mean square speed to compute acceleration from the force and as a conversion factor.
    // Dimension analysis: (ps * 1 / Å) * (Å^2 / ps^2) = Å / ps.
    return 0.5 * time_step * force * meanSquareSpeedComponent(mass);
}

/**
 * @param velocity Initial velocity
 * @param mass Particle Mass in g/mol
 * @return Updated velocity
 */
inline Point LangevinVelocityVerlet::velocityFluctuationDissipation(const Point &velocity, const double mass) {
    const double prefactor = std::exp(-friction_coefficient * time_step); // Ornstein-Uhlenbeck process prefactor
    return (prefactor * velocity) +
           random_vector(random.engine) * std::sqrt((1.0 - prefactor * prefactor) * meanSquareSpeedComponent(mass));
}

/**
 * @param velocities Vector of velocities
 * @param forces Vector of forces
 * @note
 * Using rangesv3, raw loops can be avoided and allow for future c++17 execution policies:
 *
 *     auto rng = zip(spc.p, velocities, forces);
 *     for (auto&& [particle, velicity, force] : rng) { ... };
 *     std::for_each(rng.begin(), rng.end(), [](auto &&tuple){
 *         auto&& [particle, velicity, force] = tuple;
 *         ...
 *     }};
 * @todo Splitting scheme still hard-coded to 'BAOAB'
 */
void LangevinVelocityVerlet::step(PointVector &velocities, PointVector &forces) {
    assert(spc.numParticles(Space::ACTIVE) == forces.size());
    assert(spc.numParticles(Space::ACTIVE) == velocities.size());

    auto zipped = ranges::views::zip(spc.activeParticles(), forces, velocities);

    for (auto &&[particle, force, velocity] : zipped) {
        const auto mass = particle.traits().mw;
        velocity += velocityIncrement(force, mass);                // B step
        particle.pos += positionIncrement(velocity);               // A step
        velocity = velocityFluctuationDissipation(velocity, mass); // O step
        particle.pos += positionIncrement(velocity);               // A step
        spc.geo.boundary(particle.pos);
    }
    std::fill(forces.begin(), forces.end(), Point::Zero()); // forces must be updated ...
    energy.force(forces);                                   // ... before each B step

    for (auto &&[particle, force, velocity] : zipped) {
        velocity += velocityIncrement(force, particle.traits().mw); // B step
    }
}

void LangevinVelocityVerlet::from_json(const json &j) {
    time_step = j.at("time_step").get<double>() * 1.0_ps;
    friction_coefficient = j.at("friction").get<double>() / 1.0_ps;
}

void LangevinVelocityVerlet::to_json(json &j) const {
    j = {{"time_step", time_step / 1.0_ps}, {"friction", friction_coefficient * 1.0_ps}};
}

TEST_CASE("[Faunus] Integrator") {
    class DummyEnergy : public Energy::Energybase {
        double energy(Change &) override { return 0.0; }
    };
    Space spc;
    DummyEnergy energy;

    auto ld = LangevinVelocityVerlet(spc, energy, R"({"friction": 8.0, "time_step": 2.0 })"_json);
    json j_ld = ld;
    CHECK_EQ(j_ld["friction"], 8.0);
    CHECK_EQ(j_ld["time_step"], 2.0);
}

// =============== ForceMoveBase ===============

ForceMoveBase::ForceMoveBase(Space &spc, std::shared_ptr<IntegratorBase> integrator, unsigned int nsteps)
    : integrator(integrator), number_of_steps(nsteps), spc(spc) {
    forces.reserve(spc.p.size());
    velocities.reserve(spc.p.size());
    resizeForcesAndVelocities();
    repeat = 1;
}

/**
 * @return Number of activate particles
 *
 * Upon resizing, new elements in `forces` and `velocities` are zeroed.
 */
size_t ForceMoveBase::resizeForcesAndVelocities() {
    const auto num_active_particles = spc.numParticles(Space::ACTIVE);
    forces.resize(num_active_particles, Point::Zero());
    velocities.resize(num_active_particles, Point::Zero());
    return num_active_particles;
}

void ForceMoveBase::_move(Change &change) {
    change.clear();
    change.all = true;
    resizeForcesAndVelocities();
    for (unsigned int step = 0; step < number_of_steps; ++step) {
        integrator->step(velocities, forces);
    }
    for (auto &group : spc.groups) { // update mass centers before returning to Monte Carlo
        group.updateMassCenter(spc.geo.getBoundaryFunc());
    }
}

void ForceMoveBase::_to_json(json &j) const {
    j = {{"nsteps", number_of_steps}};
    j["integrator"] = *integrator;
}

void ForceMoveBase::_from_json(const json &j) {
    number_of_steps = j.at("nsteps").get<unsigned int>();
    integrator->from_json(j["integrator"]);
    generateVelocities();
}

double ForceMoveBase::bias(Change &, double, double) {
    return pc::neg_infty; // always accept the move
}

void ForceMoveBase::generateVelocities() {
    NormalRandomVector random_vector; // generator of random 3d vector from a normal distribution
    const auto particles = spc.activeParticles();
    resizeForcesAndVelocities();
    std::transform(particles.begin(), particles.end(), velocities.begin(), [&](auto &particle) {
        return random_vector(random.engine) * std::sqrt(meanSquareSpeedComponent(particle.traits().mw));
    });
    std::fill(forces.begin(), forces.end(), Point::Zero());
}

const PointVector &ForceMoveBase::getForces() const { return forces; }
const PointVector &ForceMoveBase::getVelocities() const { return velocities; }

// =============== LangevinMove ===============

LangevinDynamics::LangevinDynamics(Space &spc, std::shared_ptr<IntegratorBase> integrator, unsigned int nsteps)
    : ForceMoveBase::ForceMoveBase(spc, integrator, nsteps) {
    name = "langevin_dynamics";
}

LangevinDynamics::LangevinDynamics(Space &spc, Energy::Energybase &energy)
    : LangevinDynamics::LangevinDynamics(spc, std::make_shared<LangevinVelocityVerlet>(spc, energy), 0) {}

LangevinDynamics::LangevinDynamics(Space &spc, Energy::Energybase &energy, const json &j)
    : LangevinDynamics::LangevinDynamics(spc, energy) {
    from_json(j);
}

void LangevinDynamics::_to_json(json &j) const {
    ForceMoveBase::_to_json(j);

}
void LangevinDynamics::_from_json(const json &j) {
    ForceMoveBase::_from_json(j);
}

TEST_CASE("[Faunus] LangevinDynamics") {
    class DummyEnergy : public Energy::Energybase {
        double energy(Change &) override { return 0.0; }
    };
    Space spc;
    DummyEnergy energy;

    SUBCASE("Velocity and force initialization") {
        spc.p.resize(10);                                        // 10 particles in total
        spc.groups.emplace_back(spc.p.begin(), spc.p.end() - 1); // 9 active particles
        LangevinDynamics ld(spc, energy);
        CHECK(ld.getForces().capacity() >= 10);
        CHECK(ld.getVelocities().capacity() >= 10);
        CHECK_EQ(ld.getForces().size(), 9);
        CHECK_EQ(ld.getVelocities().size(), 9);
    }

    SUBCASE("JSON init") {
        json j_in = R"({"nsteps": 100, "integrator": {"time_step": 0.001, "friction": 2.0}})"_json;
        LangevinDynamics ld(spc, energy);
        ld.from_json(j_in);
        json j_out = ld;
        // CHECK_EQ(j_out, j_in);
    }
}

TEST_SUITE_END();

} // end of namespace
