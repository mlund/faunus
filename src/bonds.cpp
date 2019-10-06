#include "bonds.h"
#include "geometry.h"

namespace Faunus {
namespace Potential {
void to_json(Faunus::json &j, const std::shared_ptr<BondData> &b) {
    json val;
    b->to_json(val);
    val["index"] = b->index;
    j = {{b->name(), val}};
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
    else if (b->type() == BondData::G96_TORSION)
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

BondData::~BondData() {}

void HarmonicBond::from_json(const Faunus::json &j) {
    k = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2) / 2; // k
    req = j.at("req").get<double>() * 1.0_angstrom;                          // req
}

void HarmonicBond::to_json(Faunus::json &j) const {
    j = {{"k", 2 * k / 1.0_kJmol * 1.0_angstrom * 1.0_angstrom}, {"req", req / 1.0_angstrom}};
}

std::string HarmonicBond::name() const { return "harmonic"; }

std::shared_ptr<BondData> HarmonicBond::clone() const { return std::make_shared<HarmonicBond>(*this); }

int HarmonicBond::numindex() const { return 2; }

BondData::Variant HarmonicBond::type() const { return BondData::HARMONIC; }

void HarmonicBond::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        double d = req - dist(p[index[0]].pos, p[index[1]].pos).norm();
        return k * d * d;
    };
}

std::shared_ptr<BondData> FENEBond::clone() const { return std::make_shared<FENEBond>(*this); }

int FENEBond::numindex() const { return 2; }

BondData::Variant FENEBond::type() const { return BondData::FENE; }

void FENEBond::from_json(const Faunus::json &j) {
    k[0] = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2);
    k[1] = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2);
}

void FENEBond::to_json(Faunus::json &j) const {
    j = {{"k", k[0] / (1.0_kJmol / std::pow(1.0_angstrom, 2))}, {"rmax", std::sqrt(k[1]) / 1.0_angstrom}};
}

std::string FENEBond::name() const { return "fene"; }
void FENEBond::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        double d = dist(p[index[0]].pos, p[index[1]].pos).squaredNorm();
        return (d > k[1]) ? pc::infty : -0.5 * k[0] * k[1] * std::log(1 - d / k[1]);
    };
}

std::shared_ptr<BondData> FENEWCABond::clone() const { return std::make_shared<FENEWCABond>(*this); }

int FENEWCABond::numindex() const { return 2; }

BondData::Variant FENEWCABond::type() const { return BondData::FENEWCA; }

void FENEWCABond::from_json(const Faunus::json &j) {
    k[0] = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2);
    k[1] = std::pow(j.at("rmax").get<double>() * 1.0_angstrom, 2);
    k[2] = j.at("eps").get<double>() * 1.0_kJmol;
    k[3] = std::pow(j.at("sigma").get<double>() * 1.0_angstrom, 2);
}

void FENEWCABond::to_json(Faunus::json &j) const {
    j = {{"k", k[0] / (1.0_kJmol / std::pow(1.0_angstrom, 2))},
         {"rmax", std::sqrt(k[1]) / 1.0_angstrom},
         {"eps", k[2] / 1.0_kJmol},
         {"sigma", std::sqrt(k[3]) / 1.0_angstrom}};
}

std::string FENEWCABond::name() const { return "fene+wca"; }
void FENEWCABond::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        double wca = 0, d = dist(p[index[0]].pos, p[index[1]].pos).squaredNorm();
        double x = k[3];
        if (d <= x * 1.2599210498948732) {
            x = x / d;
            x = x * x * x;
            wca = k[2] * (x * x - x + 0.25);
        }
        return (d > k[1]) ? pc::infty : -0.5 * k[0] * k[1] * std::log(1 - d / k[1]) + wca;
    };
}

int HarmonicTorsion::numindex() const { return 3; }

void HarmonicTorsion::from_json(const Faunus::json &j) {
    k = j.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_rad, 2);
    aeq = j.at("aeq").get<double>() * 1.0_deg;
}

void HarmonicTorsion::to_json(Faunus::json &j) const {
    j = {{"k", k / (1.0_kJmol / std::pow(1.0_rad, 2))}, {"aeq", aeq / 1.0_deg}};
    _roundjson(j, 6);
}

BondData::Variant HarmonicTorsion::type() const { return BondData::HARMONIC_TORSION; }

std::string HarmonicTorsion::name() const { return "harmonic_torsion"; }

std::shared_ptr<BondData> HarmonicTorsion::clone() const { return std::make_shared<HarmonicTorsion>(*this); }
void HarmonicTorsion::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        Point ray1 = dist(p[index[0]].pos, p[index[1]].pos);
        Point ray2 = dist(p[index[2]].pos, p[index[1]].pos);
        double angle = std::acos(ray1.dot(ray2) / ray1.norm() / ray2.norm());
        return 0.5 * k * (angle - aeq) * (angle - aeq);
    };
}

int GromosTorsion::numindex() const { return 3; }

void GromosTorsion::from_json(const Faunus::json &j) {
    k = j.at("k").get<double>() * 0.5 * 1.0_kJmol;  // k
    aeq = cos(j.at("aeq").get<double>() * 1.0_deg); // cos(angle)
}

void GromosTorsion::to_json(Faunus::json &j) const {
    j = {{"k", 2.0 * k / 1.0_kJmol}, {"aeq", std::acos(aeq) / 1.0_deg}};
}

BondData::Variant GromosTorsion::type() const { return BondData::G96_TORSION; }

std::string GromosTorsion::name() const { return "gromos_torsion"; }

std::shared_ptr<BondData> GromosTorsion::clone() const { return std::make_shared<GromosTorsion>(*this); }
void GromosTorsion::setEnergyFunction(const ParticleVector &p) {
    energy = [&](Geometry::DistanceFunction dist) {
        Point ray1 = dist(p[index[0]].pos, p[index[1]].pos);
        Point ray2 = dist(p[index[2]].pos, p[index[1]].pos);
        double dangle = aeq - std::acos(ray1.dot(ray2) / ray1.norm() / ray2.norm());
        return k * dangle * dangle;
    };
}

int PeriodicDihedral::numindex() const { return 4; }

std::shared_ptr<BondData> PeriodicDihedral::clone() const { return std::make_shared<PeriodicDihedral>(*this); }

void PeriodicDihedral::from_json(const Faunus::json &j) {
    k[0] = j.at("k").get<double>() * 1.0_kJmol; // k
    k[1] = j.at("n").get<double>();             // multiplicity/periodicity n
    k[2] = j.at("phi").get<double>() * 1.0_deg; // angle
}

void PeriodicDihedral::to_json(Faunus::json &j) const {
    j = {{"k", k[0] / 1.0_kJmol}, {"n", k[1]}, {"phi", k[2] / 1.0_deg}};
}

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
        return k[0] * (1 + cos(k[1] * angle - k[2]));
    };
}

} // namespace Potential
} // namespace Faunus
