#include "bonds.h"
#include "geometry.h"

namespace Faunus {
namespace Potential {

void to_json(Faunus::json &j, const std::shared_ptr<const BondData> &b) {
    to_json(j, *b);
}

void to_json(Faunus::json &j, const BondData &b) {
    json val;
    b.to_json(val);
    val["index"] = b.index;
    j = {{b.name(), val}};
}

void from_json(const Faunus::json &j, std::shared_ptr<BondData> &b) {
    if (j.is_object())
        if (j.size() == 1) {
            const auto &key = j.begin().key();
            const auto &val = j.begin().value();
            if (key == HarmonicBond().name())
                b = std::make_shared<HarmonicBond>();
            else if (key == FENEBond().name())
                b = std::make_shared<FENEBond>();
            else if (key == FENEWCABond().name())
                b = std::make_shared<FENEWCABond>();
            else if (key == HarmonicTorsion().name())
                b = std::make_shared<HarmonicTorsion>();
            else if (key == GromosTorsion().name())
                b = std::make_shared<GromosTorsion>();
            else if (key == PeriodicDihedral().name())
                b = std::make_shared<PeriodicDihedral>();
            else
                throw std::runtime_error("unknown bond type: " + key);
            try {
                b->from_json(val);
                b->index = val.at("index").get<decltype(b->index)>();
                if (b->index.size() != b->numindex())
                    throw std::runtime_error("exactly " + std::to_string(b->numindex()) + " indices required for " +
                                             b->name());
            } catch (std::exception &e) {
                throw std::runtime_error(e.what() + usageTip[key]);
            }
            return;
        }
    throw std::runtime_error("invalid bond data");
}

void setBondEnergyFunction(std::shared_ptr<BondData> &b, const ParticleVector &p) {
    if (b->type() == BondData::HARMONIC)
        std::dynamic_pointer_cast<HarmonicBond>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::FENE)
        std::dynamic_pointer_cast<FENEBond>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::FENEWCA)
        std::dynamic_pointer_cast<FENEWCABond>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::HARMONIC_TORSION)
        std::dynamic_pointer_cast<HarmonicTorsion>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::GROMOS_TORSION)
        std::dynamic_pointer_cast<GromosTorsion>(b)->setEnergyFunction(p);
    else if (b->type() == BondData::PERIODIC_DIHEDRAL)
        std::dynamic_pointer_cast<PeriodicDihedral>(b)->setEnergyFunction(p);
    else {
        assert(false); // we should never reach here
    }
}

void BondData::shift(int offset) {
    for (auto &i : index)
        i += offset;
}

bool BondData::hasEnergyFunction() const { return energy != nullptr; }

BondData::BondData(const std::vector<int> &index) : index(index) {}

HarmonicBond::HarmonicBond(double k, double req, const std::vector<int> &index)
    : StretchData(index), k_half(k / 2), req(req) {}

void HarmonicBond::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2; // k
    req = j.at("req").get<double>() * 1.0_angstrom;                               // req
}

void HarmonicBond::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / 1.0_kJmol * 1.0_angstrom * 1.0_angstrom}, {"req", req / 1.0_angstrom}};
}

std::string HarmonicBond::name() const { return "harmonic"; }

std::shared_ptr<BondData> HarmonicBond::clone() const { return std::make_shared<HarmonicBond>(*this); }

BondData::Variant HarmonicBond::type() const { return BondData::HARMONIC; }

void HarmonicBond::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        double d = req - dist(p[index[0]].pos, p[index[1]].pos).norm();
        return k_half * d * d;
    };
}

FENEBond::FENEBond(double k, double rmax, const std::vector<int> &index)
    : StretchData(index), k_half(k / 2), rmax_squared(rmax * rmax) {}

std::shared_ptr<BondData> FENEBond::clone() const { return std::make_shared<FENEBond>(*this); }

BondData::Variant FENEBond::type() const { return BondData::FENE; }

void FENEBond::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2; // k
    rmax_squared = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2); // rmax
}

void FENEBond::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / (1.0_kJmol / std::pow(1.0_angstrom, 2))}, {"rmax", std::sqrt(rmax_squared) / 1.0_angstrom}};
}

std::string FENEBond::name() const { return "fene"; }

void FENEBond::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        double r_squared = dist(p[index[0]].pos, p[index[1]].pos).squaredNorm();
        return (r_squared >= rmax_squared) ? pc::infty
                                           : -k_half * rmax_squared * std::log(1 - r_squared / rmax_squared);
    };
}

FENEWCABond::FENEWCABond(double k, double rmax, double epsilon, double sigma, const std::vector<int> &index)
    : StretchData(index), k_half(k / 2), rmax_squared(rmax * rmax), epsilon(epsilon), sigma_squared(sigma * sigma) {}

std::shared_ptr<BondData> FENEWCABond::clone() const { return std::make_shared<FENEWCABond>(*this); }

BondData::Variant FENEWCABond::type() const { return BondData::FENEWCA; }

void FENEWCABond::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2;
    rmax_squared = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2);
    epsilon = j.at("eps").get<double>() * 1.0_kJmol;
    sigma_squared = std::pow(j.at("sigma").get<double>() * 1.0_angstrom, 2);
}

void FENEWCABond::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / (1.0_kJmol / std::pow(1.0_angstrom, 2))},
         {"rmax", std::sqrt(rmax_squared) / 1.0_angstrom},
         {"eps", epsilon / 1.0_kJmol},
         {"sigma", std::sqrt(sigma_squared) / 1.0_angstrom}};
}

std::string FENEWCABond::name() const { return "fene+wca"; }
void FENEWCABond::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        double r_squared = dist(p[index[0]].pos, p[index[1]].pos).squaredNorm();
        double wca = 0;
        double x = sigma_squared;
        if (r_squared <= x * 1.2599210498948732) {
            x = x / r_squared;
            x = x * x * x;
            wca = epsilon * (x * x - x + 0.25);
        }
        return (r_squared > rmax_squared) ? pc::infty
                                          : -k_half * rmax_squared * std::log(1 - r_squared / rmax_squared) + wca;
    };
}

void HarmonicTorsion::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_rad, 2) / 2;
    aeq = j.at("aeq").get<double>() * 1.0_deg;
}

void HarmonicTorsion::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / (1.0_kJmol / std::pow(1.0_rad, 2))}, {"aeq", aeq / 1.0_deg}};
    _roundjson(j, 6);
}

HarmonicTorsion::HarmonicTorsion(double k, double aeq, const std::vector<int> &index)
    : TorsionData(index), k_half(k / 2), aeq(aeq) {}

BondData::Variant HarmonicTorsion::type() const { return BondData::HARMONIC_TORSION; }

std::string HarmonicTorsion::name() const { return "harmonic_torsion"; }

std::shared_ptr<BondData> HarmonicTorsion::clone() const { return std::make_shared<HarmonicTorsion>(*this); }

void HarmonicTorsion::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        Point ray1 = dist(p[index[0]].pos, p[index[1]].pos);
        Point ray2 = dist(p[index[2]].pos, p[index[1]].pos);
        double angle = std::acos(ray1.dot(ray2) / ray1.norm() / ray2.norm());
        return k_half * (angle - aeq) * (angle - aeq);
    };
}

void GromosTorsion::from_json(const Faunus::json &j) {
    k_half = j.at("k").get<double>() * 1.0_kJmol / 2;   // k
    cos_aeq = cos(j.at("aeq").get<double>() * 1.0_deg); // cos(angle)
}

void GromosTorsion::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k_half / 1.0_kJmol}, {"aeq", std::acos(cos_aeq) / 1.0_deg}};
}

GromosTorsion::GromosTorsion(double k, double cos_aeq, const std::vector<int> &index)
    : TorsionData(index), k_half(k / 2), cos_aeq(cos_aeq) {}

BondData::Variant GromosTorsion::type() const { return BondData::GROMOS_TORSION; }

std::string GromosTorsion::name() const { return "gromos_torsion"; }

std::shared_ptr<BondData> GromosTorsion::clone() const { return std::make_shared<GromosTorsion>(*this); }

void GromosTorsion::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        Point ray1 = dist(p[index[0]].pos, p[index[1]].pos);
        Point ray2 = dist(p[index[2]].pos, p[index[1]].pos);
        double dcos = cos_aeq - ray1.dot(ray2) / (ray1.norm() * ray2.norm());
        return k_half * dcos * dcos;
    };
}

int PeriodicDihedral::numindex() const { return 4; }

std::shared_ptr<BondData> PeriodicDihedral::clone() const { return std::make_shared<PeriodicDihedral>(*this); }

void PeriodicDihedral::from_json(const Faunus::json &j) {
    k = j.at("k").get<double>() * 1.0_kJmol;   // k
    n = j.at("n").get<double>();               // multiplicity/periodicity n
    phi = j.at("phi").get<double>() * 1.0_deg; // angle
}

void PeriodicDihedral::to_json(Faunus::json &j) const {
    j = {{"k", k / 1.0_kJmol}, {"n", n}, {"phi", phi / 1.0_deg}};
}

PeriodicDihedral::PeriodicDihedral(double k, double phi, double n, const std::vector<int> &index)
    : BondData(index), k(k), phi(phi), n(n) {}

BondData::Variant PeriodicDihedral::type() const { return BondData::PERIODIC_DIHEDRAL; }

std::string PeriodicDihedral::name() const { return "periodic_dihedral"; }

void PeriodicDihedral::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        Point vec1 = dist(p[index[1]].pos, p[index[0]].pos);
        Point vec2 = dist(p[index[2]].pos, p[index[1]].pos);
        Point vec3 = dist(p[index[3]].pos, p[index[2]].pos);
        Point norm1 = vec1.cross(vec2);
        Point norm2 = vec2.cross(vec3);
        // atan2( [v1×v2]×[v2×v3]⋅[v2/|v2|], [v1×v2]⋅[v2×v3] )
        double angle = atan2((norm1.cross(norm2)).dot(vec2) / vec2.norm(), norm1.dot(norm2));
        return k * (1 + cos(n * angle - phi));
    };
}

} // namespace Potential
} // namespace Faunus
