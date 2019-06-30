#include "externalpotential.h"
#include "multipole.h"

namespace Faunus {
namespace Energy {

// ------------ Energybase -------------

void Energybase::to_json(json &) const {}

void Energybase::sync(Energybase *, Change &) {}

void Energybase::init() {}

void to_json(json &j, const Energybase &base) {
    assert(not base.name.empty());
    if (base.timer)
        j["relative time"] = base.timer.result();
    if (not base.cite.empty())
        j[base.name]["reference"] = base.cite;
    base.to_json(j[base.name]);
}

// ------------ ExternalPotential -------------

// this calculates the interaction of a whole group
// with the applied external potential.
double ExternalPotential::_energy(const Group<Particle> &g) const {
    double u = 0;
    if (molids.find(g.id) != molids.end()) {
        if (COM and g.atomic == false) { // apply only to center of mass
            Particle cm;                 // temp. particle representing molecule
            cm.charge = Faunus::monopoleMoment(g.begin(), g.end());
            cm.pos = g.cm;
            return func(cm);
        } else {
            for (auto &p : g) { // loop over active particles
                u += func(p);
                if (std::isnan(u))
                    break;
            }
        }
    }
    return u;
}

ExternalPotential::ExternalPotential(const json &j, Tspace &spc) : spc(spc) {
    name = "external";
    COM = j.value("com", false);
    _names = j.at("molecules").get<decltype(_names)>(); // molecule names
    auto _ids = names2ids(molecules, _names);           // names --> molids
    molids = std::set<int>(_ids.begin(), _ids.end());   // vector --> set
    if (molids.empty() || molids.size() != _names.size())
        throw std::runtime_error(name + ": molecule list is empty");
}
double ExternalPotential::energy(Change &change) {
    assert(func != nullptr);
    double u = 0;
    if (change.dV or change.all or change.dN) {
        for (auto &g : spc.groups) { // check all groups
            u += _energy(g);
            if (std::isnan(u))
                break;
        }
    } else
        for (auto &d : change.groups) {
            auto &g = spc.groups.at(d.index); // check specified groups
            if (d.all or COM)                 // check all atoms in group
                u += _energy(g);              // _energy also checks for molecule id
            else {                            // check only specified atoms in group
                if (molids.find(g.id) != molids.end())
                    for (auto i : d.atoms)
                        u += func(*(g.begin() + i));
            }
            if (std::isnan(u))
                break;
        }
    return u;
}
void ExternalPotential::to_json(json &j) const {
    j["molecules"] = _names;
    j["com"] = COM;
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
        func = [&radius = radius, origo = origo, k = k, dir = dir](const Particle &p) {
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
        func = [low = low, high = high, k = k](const Particle &p) {
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
    cite = "doi:10/dhb9mj";

    SingleUseJSON _j = j; // json variant where items are deleted after access
    _j.erase("com");
    _j.erase("molecules");

    nstep = _j.at("nstep").get<unsigned int>();
    epsr = _j.at("epsr").get<double>();
    fixed = _j.value("fixed", false);
    nphi = _j.value("nphi", 10);

    halfz = 0.5 * spc.geo.getLength().z();
    lB = pc::lB(epsr);

    dz = _j.value("dz", 0.2); // read z resolution
    Q.setResolution(dz, -halfz, halfz);
    rho.setResolution(dz, -halfz, halfz);
    phi.setResolution(dz, -halfz, halfz);

    filename = _j.value("file", "mfcorr.dat"s);
    load();

    func = [&phi = phi](const typename Tspace::Tparticle &p) { return p.charge * phi(p.pos.z()); };

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
        std::cerr << "density file '" << filename << "' not loaded." << endl;
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
    Point L = spc.geo.getLength();
    double area = L.x() * L.y();
    if (L.x() not_eq L.y() or 0.5 * L.z() != halfz)
        throw std::runtime_error("Requires box Lx=Ly and Lz=const.");

    Q.clear();
    for (auto &g : spc.groups) // loop over all groups
        for (auto &p : g)      // ...and their active particles
            Q(p.pos.z()) += p.charge;
    for (double z = -halfz; z <= halfz; z += dz)
        rho(z) += Q(z) / area;
}

void ExternalAkesson::update_phi() {
    Point L = spc.geo.getLength();
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

std::function<double(const Particle &)> createGouyChapmanPotential(const json &j) {
    double rho;
    double c0 = j.at("ionicstrength").get<double>() * 1.0_molar; // assuming 1:1 salt, so c0=I
    double lB = pc::lB(j.at("epsr").get<double>());
    double k = 1 / (3.04 / sqrt(c0));   // hack!
    double phi0 = j.value("phi0", 0.0); // Unitless potential = beta*e*phi0
    if (std::fabs(phi0) > 1e-6)
        rho = sqrt(2 * c0 / (pc::pi * lB)) * sinh(.5 * phi0); // Evans&Wennerstrom,Colloidal Domain p
    // 138-140
    else {
        rho = 1.0 / j.value("qarea", 0.0);
        if (rho > 1e9)
            rho = j.at("rho");
        phi0 = 2. * std::asinh(rho * std::sqrt(0.5 * lB * pc::pi / c0)); //[Evans..]
    }
    double gamma0 = std::tanh(phi0 / 4); // assuming z=1  [Evans..]
    double surface_z_pos = j.value("zpos", 0.0);
    bool linearize = j.value("linearize", false);

    // return gamma function for calculation of GC potential on single particle.
    return [=](const Particle &p) {
        if (p.charge != 0) {
            double x = std::exp(-k * std::fabs(surface_z_pos - p.pos.z()));
            if (linearize)
                return p.charge * phi0 * x;
            else {
                x = gamma0 * x;
                return 2 * p.charge * std::log((1 + x) / (1 - x));
            }
        }
        return 0.0;
    };
}

// ------------ CustomExternal -------------

CustomExternal::CustomExternal(const json &j, Tspace &spc) : ExternalPotential(j, spc) {
    name = "customexternal";
    jin = j;
    auto &_j = jin["constants"];
    if (_j == nullptr)
        _j = json::object();
    _j["e0"] = pc::e0;
    _j["kB"] = pc::kB;
    _j["kT"] = pc::kT();
    _j["Nav"] = pc::Nav;
    _j["T"] = pc::temperature;
    std::string name = jin.at("function");

    // check of the custom potential match a name with a
    // predefined meaning.
    if (name == "gouychapman")
        func = createGouyChapmanPotential(_j);
    else if (name == "something") {
        // add additional potential here
        // base::func = createSomeOtherPotential(_j);
    } else {
        // if nothing found above, it is assumed that `function`
        // is a valid expression.
        expr.set(jin, {{"q", &d.q}, {"x", &d.x}, {"y", &d.y}, {"z", &d.z}});
        func = [&](const Particle &a) {
            d.x = a.pos.x();
            d.y = a.pos.y();
            d.z = a.pos.z();
            d.q = a.charge;
            return expr();
        };
    }
}
void CustomExternal::to_json(json &j) const {
    j = jin;
    ExternalPotential::to_json(j);
}

// ------------- ParticleSelfEnergy ---------------

/*
 * Upon construction, make sure the ExternalPotential base class loop
 * over all groups and particles (com=false)
 */
ParticleSelfEnergy::ParticleSelfEnergy(Space &spc, std::function<double(const Particle &)> selfEnergy)
    : ExternalPotential({{"molecules", "*"}, {"com", false}}, spc) {
    func = selfEnergy;
#ifndef NDEBUG
    // test if self energy can be called
    Particle myparticle;
    if (this->func)
        this->func(myparticle);
#endif
    name = "selfenergy";
}

} // namespace Energy
} // namespace Faunus
