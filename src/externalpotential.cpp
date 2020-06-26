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

void Energybase::sync(Energybase *, Change &) {}

void Energybase::init() {}

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
        if (act_on_mass_center and not group.atomic) { // apply only to center of mass
            if (group.size() == group.capacity()) {    // only apply if group is active
                Particle mass_center;                  // temp. particle representing molecule
                mass_center.charge = Faunus::monopoleMoment(group.begin(), group.end());
                mass_center.pos = group.cm;
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
    if (change.dV or change.all or change.dN) {
        for (auto &group : space.groups) { // loop over all groups
            energy += groupEnergy(group);
            if (not std::isfinite(energy)) {
                break; // stop summing if not finite
            }
        }
    } else {
        for (auto &group_change : change.groups) {             // loop over all changed groups
            auto &group = space.groups.at(group_change.index); // check specified groups
            if (group_change.all or act_on_mass_center) {      // check all atoms in group
                energy += groupEnergy(group);                  // groupEnergy also checks for molecule id
            } else {                                           // only specified atoms in group
                if (molecule_ids.find(group.id) != molecule_ids.end()) {
                    for (int index : group_change.atoms) { // loop over changed atoms in group
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

// ------------ Confine -------------

Confine::Confine(const json &j, Tspace &spc) : ExternalPotential(j, spc) {
    name = "confine";
    k = value_inf(j, "k") * 1.0_kJmol; // get floating point; allow inf/-inf
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
        if (scale)
            spc.scaleVolumeTriggers.push_back(
                [&radius = radius](Tspace &, double Vold, double Vnew) { radius *= std::cbrt(Vnew / Vold); });
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
    _roundjson(j, 5);
}

// ------------ ExternalAkesson -------------

ExternalAkesson::ExternalAkesson(const json &j, Tspace &spc) : ExternalPotential(j, spc) {
    name = "akesson";
    citation_information = "doi:10/dhb9mj";

    SingleUseJSON _j = j; // json variant where items are deleted after access
    _j.erase("com");
    _j.erase("molecules");

    nstep = _j.at("nstep").get<unsigned int>();
    epsr = _j.at("epsr").get<double>();
    fixed = _j.value("fixed", false);
    nphi = _j.value("nphi", 10);

    halfz = 0.5 * spc.geo.getLength().z();
    lB = pc::bjerrumLength(epsr);

    dz = _j.value("dz", 0.2); // read z resolution
    Q.setResolution(dz, -halfz, halfz);
    rho.setResolution(dz, -halfz, halfz);
    phi.setResolution(dz, -halfz, halfz);

    filename = _j.value("file", "mfcorr.dat"s);
    load();

    externalPotentialFunc = [&phi = phi](const typename Tspace::Tparticle &p) { return p.charge * phi(p.pos.z()); };

    if (not _j.empty()) // throw exception of unused/unknown keys are passed
        throw std::runtime_error("unused key(s) for '"s + name + "':\n" + _j.dump());
}

double ExternalAkesson::energy(Change &change) {
    if (not fixed)                    // pho(z) unconverged, keep sampling
        if (key == Energybase::OLD) { // only sample on accepted configs
            cnt++;
            if (cnt % nstep == 0)
                update_rho();
            if (cnt % nstep * nphi == 0)
                update_phi();
        }
    return ExternalPotential::energy(change);
}

ExternalAkesson::~ExternalAkesson() {
    // save only if still updating and if energy type is "OLD",
    // that is, accepted configurations (not trial)
    if (not fixed and key == Energybase::OLD)
        save();
}

void ExternalAkesson::to_json(json &j) const {
    j = {{"lB", lB},         {"dz", dz},       {"nphi", nphi},          {"epsr", epsr},
         {"file", filename}, {"nstep", nstep}, {"Nupdates", updatecnt}, {"fixed", fixed}};
    ExternalPotential::to_json(j);
    _roundjson(j, 5);
}

void ExternalAkesson::save() {
    std::ofstream f(filename);
    if (f) {
        f.precision(16);
        f << rho;
    } else
        throw std::runtime_error("cannot save file '"s + filename + "'");
}

void ExternalAkesson::load() {
    std::ifstream f(filename);
    if (f) {
        rho << f;
        update_phi();
    } else
        faunus_logger->warn("density file {} not loaded", filename);
}

double ExternalAkesson::phi_ext(double z, double a) const {
    double a2 = a * a, z2 = z * z;
    return -2 * pc::pi * z - 8 * a * std::log((std::sqrt(2 * a2 + z2) + a) / std::sqrt(a2 + z2)) +
           2 * z * (0.5 * pc::pi + std::asin((a2 * a2 - z2 * z2 - 2 * a2 * z2) / std::pow(a2 + z2, 2)));
}

void ExternalAkesson::sync(Energybase *basePtr, Change &) {
    if (not fixed) {
        auto other = dynamic_cast<decltype(this)>(basePtr);
        assert(other);
        // only trial energy (new) require sync
        if (other->key == Energybase::OLD)
            if (cnt != other->cnt) {
                assert(cnt < other->cnt && "trial cnt's must be smaller");
                cnt = other->cnt;
                rho = other->rho;
                phi = other->phi;
            }
    }
}

void ExternalAkesson::update_rho() {
    updatecnt++;
    Point L = space.geo.getLength();
    double area = L.x() * L.y();
    if (L.x() not_eq L.y() or 0.5 * L.z() != halfz)
        throw std::runtime_error("Requires box Lx=Ly and Lz=const.");

    Q.clear();
    for (auto &g : space.groups) // loop over all groups
        for (auto &p : g)      // ...and their active particles
            Q(p.pos.z()) += p.charge;
    for (double z = -halfz; z <= halfz; z += dz)
        rho(z) += Q(z) / area;
}

void ExternalAkesson::update_phi() {
    Point L = space.geo.getLength();
    double a = 0.5 * L.x();
    for (double z = -halfz; z <= halfz; z += dz) {
        double s = 0;
        for (double zn = -halfz; zn <= halfz; zn += dz)
            if (rho(zn).cnt > 0)
                s += rho(zn).avg() * phi_ext(std::fabs(z - zn), a); // Eq. 14 in Greberg paper
        phi(z) = lB * s;
    }
}

// ------------ createGouyChapman -------------

std::function<double(const Particle &)> createGouyChapmanPotential(const json &j, const Geometry::Chameleon &geo) {
    if (geo.boundaryConditions().direction.z() != Geometry::FIXED) {
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

// ------------ CustomExternal -------------

CustomExternal::CustomExternal(const json &j, Space &spc) : ExternalPotential(j, spc), json_input_backup(j) {
    name = "customexternal";
    auto &constants = json_input_backup["constants"];
    if (std::string function = j.at("function"); function == "gouychapman") {
        externalPotentialFunc = createGouyChapmanPotential(constants, spc.geo);
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
