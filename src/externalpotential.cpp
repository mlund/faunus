#include "externalpotential.h"
#include "multipole.h"
#include "aux/eigensupport.h"
#include "functionparser.h"
#include "space.h"
#include "spdlog/spdlog.h"

namespace Faunus {
namespace Energy {

// ------------ Energybase -------------

void Energybase::to_json(json &) const {}

/**
 * @param other_energy Other energy instance to copy data from
 * @param change Describes the difference with the other energy term
 */
void Energybase::sync([[maybe_unused]] Energybase* other_energy, [[maybe_unused]] const Change& change) {}

void Energybase::init() {}
void Energybase::force([[maybe_unused]] PointVector& forces) {}

void to_json(json &j, const Energybase &base) {
    assert(not base.name.empty());
    if (base.timer)
        j[base.name]["relative time"] = base.timer.result();
    if (not base.citation_information.empty())
        j[base.name]["reference"] = base.citation_information;
    base.to_json(j[base.name]);
}

// ------------ ExternalPotential -------------

/**
 * @param group Group to calculate energy of
 * @return Energy of group in kT
 *
 * - Calculates the interaction of a whole group with the applied external potential
 * - The group is ignored if not part of the `molecule_id_list`
 * - If `act_on_mass_center` is true, the external potential is applied on a
 *   fictitious particle placed at the COM and with a net-charge of the group.
 */
double ExternalPotential::groupEnergy(const Group<Particle> &group) const {
    double u = 0;
    if (molecule_ids.find(group.id) != molecule_ids.end()) {
        if (act_on_mass_center && group.isMolecular()) { // apply only to center of mass
            if (group.size() == group.capacity()) {    // only apply if group is active
                Particle mass_center;                  // temp. particle representing molecule
                mass_center.charge = Faunus::monopoleMoment(group.begin(), group.end());
                mass_center.pos = group.mass_center;
                return externalPotentialFunc(mass_center);
            }
        } else {
            for (auto &particle : group) { // loop over active particles
                u += externalPotentialFunc(particle);
                if (std::isnan(u)) {
                    break;
                }
            }
        }
    }
    return u;
}

ExternalPotential::ExternalPotential(const json &j, Space &spc) : space(spc) {
    name = "external";
    act_on_mass_center = j.value("com", false);
    molecule_names = j.at("molecules").get<decltype(molecule_names)>(); // molecule names
    auto _ids = Faunus::names2ids(Faunus::molecules, molecule_names);   // names --> molids
    molecule_ids = std::set<int>(_ids.begin(), _ids.end());             // vector --> set
    if (molecule_ids.empty()) {
        throw std::runtime_error(name + ": molecule list is empty");
    }
}
double ExternalPotential::energy(Change &change) {
    assert(externalPotentialFunc != nullptr);
    double energy = 0.0;
    if (change.volume_change or change.everything or change.matter_change) {
        for (auto &group : space.groups) { // loop over all groups
            energy += groupEnergy(group);
            if (not std::isfinite(energy)) {
                break; // stop summing if not finite
            }
        }
    } else {
        for (auto &group_change : change.groups) {             // loop over all changed groups
            auto& group = space.groups.at(group_change.group_index); // check specified groups
            if (group_change.all or act_on_mass_center) {      // check all atoms in group
                energy += groupEnergy(group);                  // groupEnergy also checks for molecule id
            } else {                                           // only specified atoms in group
                if (molecule_ids.find(group.id) != molecule_ids.end()) {
                    for (int index : group_change.relative_atom_indices) { // loop over changed atoms in group
                        energy += externalPotentialFunc(group[index]);
                    }
                }
            }
            if (not std::isfinite(energy)) {
                break; // stop summing if not finite
            }
        }
    }
    return energy; // in kT
}
void ExternalPotential::to_json(json &j) const {
    j["molecules"] = molecule_names;
    j["com"] = act_on_mass_center;
}

TEST_CASE("[Faunus] ExternalPotential") {
    using doctest::Approx;
    Faunus::atoms = R"([
        { "A": { "sigma": 4.0, "tfe": 1.0 } },
        { "B": { "sigma": 2.4, "tfe": 1.0 } }
    ])"_json.get<decltype(atoms)>();

    Faunus::molecules = R"([
        { "M": { "atoms": ["A", "B"], "atomic": true } }
    ])"_json.get<decltype(molecules)>();

    json j = R"({
        "geometry": {"type": "sphere", "radius": 100 },
        "insertmolecules": [ { "M": { "N": 1 } } ]
    })"_json;

    SUBCASE("ParticleSelfEnergy") {
        Space spc = j;
        ParticleSelfEnergy pot(spc, [](const Particle &) { return 0.5; });
        Change change;
        change.everything = true; // if both particles have changed
        CHECK(pot.energy(change) == Approx(0.5 + 0.5));
    }
}

// ------------ Confine -------------

Confine::Confine(const json& j, Space& spc) : ExternalPotential(j, spc) {
    name = "confine";
    k = getValueInfinity(j, "k") * 1.0_kJmol; // get floating point; allow inf/-inf
    type = m.at(j.at("type"));

    if (type == sphere or type == cylinder) {
        radius = j.at("radius");
        origo = j.value("origo", origo);
        scale = j.value("scale", scale);
        if (type == cylinder)
            dir = {1, 1, 0};
        externalPotentialFunc = [&radius = radius, origo = origo, k = k, dir = dir](const Particle &p) {
            double d2 = (origo - p.pos).cwiseProduct(dir).squaredNorm() - radius * radius;
            if (d2 > 0)
                return 0.5 * k * d2;
            return 0.0;
        };

        // If volume is scaled, also scale the confining radius by adding a trigger
        // to `Space::scaleVolume()`
        if (scale) {
            spc.scaleVolumeTriggers.push_back([&radius = radius]([[maybe_unused]] Space& spc, double Vold,
                                                                 double Vnew) { radius *= std::cbrt(Vnew / Vold); });
        }
    }

    if (type == cuboid) {
        low = j.at("low").get<Point>();
        high = j.at("high").get<Point>();
        externalPotentialFunc = [low = low, high = high, k = k](const Particle &p) {
            double u = 0;
            Point d = low - p.pos;
            for (int i = 0; i < 3; ++i)
                if (d[i] > 0)
                    u += d[i] * d[i];
            d = p.pos - high;
            for (int i = 0; i < 3; ++i)
                if (d[i] > 0)
                    u += d[i] * d[i];
            return 0.5 * k * u;
        };
    }
}
void Confine::to_json(json &j) const {
    if (type == cuboid)
        j = {{"low", low}, {"high", high}};
    if (type == sphere or type == cylinder)
        j = {{"radius", radius}};
    if (type == sphere) {
        j["origo"] = origo;
        j["scale"] = scale;
    }
    for (auto &i : m)
        if (i.second == type)
            j["type"] = i.first;
    j["k"] = k / 1.0_kJmol;
    ExternalPotential::to_json(j);
    roundJSON(j, 5);
}

// ------------ ExternalAkesson -------------

ExternalAkesson::ExternalAkesson(const json& j, Space& spc) : ExternalPotential(j, spc) {
    name = "akesson";
    citation_information = "doi:10/dhb9mj";
    nstep = j.at("nstep").get<unsigned int>();
    dielectric_constant = j.at("epsr").get<double>();
    fixed_potential = j.value("fixed", false);
    phi_update_interval = j.value("nphi", 10);

    half_box_length_z = 0.5 * spc.geometry.getLength().z();
    bjerrum_length = pc::bjerrumLength(dielectric_constant);

    dz = j.value("dz", 0.2); // read z resolution
    charge_profile.setResolution(dz, -half_box_length_z, half_box_length_z);
    rho.setResolution(dz, -half_box_length_z, half_box_length_z);
    phi.setResolution(dz, -half_box_length_z, half_box_length_z);

    filename = j.value("file", "mfcorr.dat"s);
    load_rho();
    externalPotentialFunc = [&phi = phi](const Particle& particle) { return particle.charge * phi(particle.pos.z()); };
}

double ExternalAkesson::energy(Change &change) {
    if (not fixed_potential) {              // phi(z) unconverged, keep sampling
        if (state == MonteCarloState::ACCEPTED) { // only sample on accepted configs
            num_density_updates++;
            if (num_density_updates % nstep == 0) {
                update_rho();
            }
            if (num_density_updates % nstep * phi_update_interval == 0) {
                update_phi();
            }
        }
    }
    return ExternalPotential::energy(change);
}

ExternalAkesson::~ExternalAkesson() {
    // save only if still updating and if energy type is `ACCEPTED_MONTE_CARLO_STATE`,
    // that is, accepted configurations (not trial)
    if (not fixed_potential and state == MonteCarloState::ACCEPTED) {
        save_rho();
    }
}

void ExternalAkesson::to_json(json &j) const {
    j = {{"lB", bjerrum_length}, {"dz", dz},       {"nphi", phi_update_interval}, {"epsr", dielectric_constant},
         {"file", filename},     {"nstep", nstep}, {"Nupdates", num_rho_updates}, {"fixed", fixed_potential}};
    ExternalPotential::to_json(j);
    roundJSON(j, 5);
}

void ExternalAkesson::save_rho() {
    if (auto stream = std::ofstream(filename); stream) {
        stream.precision(16);
        stream << rho;
    } else
        throw std::runtime_error("cannot save file '"s + filename + "'");
}

void ExternalAkesson::load_rho() {
    if (auto stream = std::ifstream(filename); stream) {
        rho << stream;
        update_phi();
    } else {
        faunus_logger->warn("density file {} not loaded", filename);
    }
}

/**
 * This is Eq. 15 of the mol. phys. 1996 paper by Greberg et al.
 * (sign typo in manuscript: phi^infty(z) should be "-2*pi*z" on page 413, middle)
 */
double ExternalAkesson::phi_ext(double z, double a) const {
    double a2 = a * a, z2 = z * z;
    return -2 * pc::pi * z - 8 * a * std::log((std::sqrt(2 * a2 + z2) + a) / std::sqrt(a2 + z2)) +
           2 * z * (0.5 * pc::pi + std::asin((a2 * a2 - z2 * z2 - 2 * a2 * z2) / std::pow(a2 + z2, 2)));
}

void ExternalAkesson::sync(Energybase* energybase, const Change&) {
    if (fixed_potential) {
        return;
    }
    if (auto* other = dynamic_cast<ExternalAkesson*>(energybase)) {
        if (other->state == MonteCarloState::ACCEPTED) { // only trial energy (new) requires sync
            if (num_density_updates != other->num_density_updates) {
                assert(num_density_updates < other->num_density_updates);
                num_density_updates = other->num_density_updates;
                rho = other->rho;
                phi = other->phi;
            }
        }
    } else {
        throw std::runtime_error("akesson sync error");
    }
}

void ExternalAkesson::update_rho() {
    num_rho_updates++;
    Point L = space.geometry.getLength();
    if (L.x() not_eq L.y() or 0.5 * L.z() != half_box_length_z) {
        throw std::runtime_error("Requires box Lx=Ly and Lz=const.");
    }
    charge_profile.clear();
    for (auto &group : space.groups) { // loop over all groups
        for (auto &particle : group) { // ...and their active particles
            charge_profile(particle.pos.z()) += particle.charge;
        }
    }
    double area = L.x() * L.y();
    for (double z = -half_box_length_z; z <= half_box_length_z; z += dz) {
        rho(z) += charge_profile(z) / area;
    }
}

void ExternalAkesson::update_phi() {
    auto L = space.geometry.getLength();
    double a = 0.5 * L.x();
    for (double z = -half_box_length_z; z <= half_box_length_z; z += dz) {
        double s = 0;
        for (double zn = -half_box_length_z; zn <= half_box_length_z; zn += dz) {
            if (!rho(zn).empty()) {
                s += rho(zn).avg() * phi_ext(std::fabs(z - zn), a); // Eq. 14 in Greberg's paper
            }
        }
        phi(z) = bjerrum_length * s;
    }
}

// ------------ createGouyChapman -------------

std::function<double(const Particle &)> createGouyChapmanPotential(const json &j, const Geometry::Chameleon &geo) {
    if (geo.boundaryConditions().direction.z() != Geometry::Boundary::FIXED) {
        throw std::runtime_error("Gouy-Chapman requires non-periodicity in z-direction");
    }
    double rho = 0; // surface charge density (charge per area)
    double bjerrum_length = pc::bjerrumLength(j.at("epsr").get<double>());
    double molarity = j.at("molarity").get<double>();
    double kappa = 1.0 / Faunus::debyeLength(molarity, {1, 1}, bjerrum_length);
    double phi0 = j.value("phi0", 0.0); // Unitless potential = beta*e*phi0
    if (std::fabs(phi0) > 0) {
        rho = std::sqrt(2.0 * molarity / (pc::pi * bjerrum_length)) *
              std::sinh(0.5 * phi0); // Evans&Wennerstrom,Colloidal Domain p. 138-140
    } else {                         // phi0 was not provided
        double area_per_charge = j.value("rhoinv", 0.0);
        if (std::fabs(area_per_charge) > 0) {
            rho = 1.0 / area_per_charge;
        } else {
            rho = j.at("rho").get<double>();
        }
        phi0 = 2.0 * std::asinh(rho * std::sqrt(0.5 * bjerrum_length * pc::pi / molarity)); // [Evans..]
    }
    double gamma0 = std::tanh(phi0 / 4.0); // assuming z=1 [Evans..]

    faunus_logger->trace("generated Gouy-Chapman potential with {} A^2/charge ", 1.0 / rho);

    if (j.value("linearise", false)) {
        return [=, &geo](const Particle &p) {
            double surface_z_pos = -0.5 * geo.getLength().z();
            return p.charge * phi0 * std::exp(-kappa * std::fabs(surface_z_pos - p.pos.z()));
        };
    } else {
        return [=, &geo](const Particle &p) {
            double surface_z_pos = -0.5 * geo.getLength().z();
            double x = gamma0 * std::exp(-kappa * std::fabs(surface_z_pos - p.pos.z()));
            return 2.0 * p.charge * std::log((1.0 + x) / (1.0 - x));
        };
    }
}

TEST_CASE("[Faunus] Gouy-Chapman") {
    using doctest::Approx;
    Geometry::Slit slit(50, 50, 50);
    Geometry::Chameleon geometry(slit, Geometry::Variant::SLIT);
    json j = {{"molarity", 0.1}, {"epsr", 80}, {"linearise", false}, {"rhoinv", 100.0}};
    auto phi = Energy::createGouyChapmanPotential(j, geometry);
    Particle p;
    p.charge = 1.0;
    p.pos = {0, 0, -25};                            // potential at charged surface
    CHECK(phi(p) == doctest::Approx(0.2087776151)); // = phi_0

    p.pos = {0, 0, 0}; // potential at mid-plane
    CHECK(phi(p) == doctest::Approx(0.0160227029));

    j = {{"molarity", 0.1}, {"epsr", 80}, {"linearise", false}, {"phi0", 0.2087776151}};
    phi = Energy::createGouyChapmanPotential(j, geometry);
    CHECK(phi(p) == doctest::Approx(0.0160227029));

    j = {{"molarity", 0.1}, {"epsr", 80}, {"linearise", false}, {"rho", 0.01}};
    phi = Energy::createGouyChapmanPotential(j, geometry);
    CHECK(phi(p) == doctest::Approx(0.0160227029));

    j = {{"molarity", 0.1}, {"epsr", 80}, {"linearise", true}, {"rho", 0.01}};
    phi = Energy::createGouyChapmanPotential(j, geometry);
    CHECK(phi(p) == doctest::Approx(0.0160371645));
}

// ------------ CustomExternal -------------

CustomExternal::CustomExternal(const json &j, Space &spc) : ExternalPotential(j, spc), json_input_backup(j) {
    name = "customexternal";
    auto &constants = json_input_backup["constants"];
    if (std::string function = j.at("function"); function == "gouychapman") {
        externalPotentialFunc = createGouyChapmanPotential(constants, spc.geometry);
    } else if (function == "some-new-potential") { // add new potentials here
        // func = createSomeNewPotential(...);
    } else { // nothing found above; assume `function` is an expression
        if (constants == nullptr) {
            constants = json::object();
        }
        constants["e0"] = pc::e0;
        constants["kB"] = pc::kB;
        constants["kT"] = pc::kT();
        constants["Nav"] = pc::Nav;
        constants["T"] = pc::temperature;
        expr = std::make_unique<ExprFunction<double>>();
        expr->set(
            j,
            {{"q", &particle_data.charge}, {"x", &particle_data.x}, {"y", &particle_data.y}, {"z", &particle_data.z}});
        externalPotentialFunc = [&](const Particle &a) {
            particle_data.x = a.pos.x();
            particle_data.y = a.pos.y();
            particle_data.z = a.pos.z();
            particle_data.charge = a.charge;
            return expr->operator()();
        };
    }
}
void CustomExternal::to_json(json &j) const {
    j = json_input_backup;
    ExternalPotential::to_json(j);
}

// ------------- ParticleSelfEnergy ---------------

/*
 * Upon construction, make sure the ExternalPotential base class loop
 * over all groups and particles (com=false)
 */
ParticleSelfEnergy::ParticleSelfEnergy(Space &spc, std::function<double(const Particle &)> selfEnergy)
    : ExternalPotential({{"molecules", {"*"}}, {"com", false}}, spc) {
    assert(selfEnergy && "selfEnergy is not callable");
    externalPotentialFunc = selfEnergy;
#ifndef NDEBUG
    // test if self energy can be called
    assert(not Faunus::atoms.empty());
    Particle myparticle;
    myparticle.id=0;
    if (this->externalPotentialFunc) {
        double u = this->externalPotentialFunc(myparticle);
        assert(std::isfinite(u));
    }
#endif
    name = "particle-self-energy";
}

} // namespace Energy
} // namespace Faunus
