#include "particle.h"
#include <Eigen/Geometry>

namespace Faunus {
void ParticlePropertyBase::rotate(const Eigen::Quaterniond &, const Eigen::Matrix3d &) {}

void Radius::to_json(json &j) const { j["r"] = radius; }

void Radius::from_json(const json &j) { radius = j.value("r", 0.0); }

void Charge::to_json(json &j) const { j["q"] = charge; }

void Charge::from_json(const json &j) { charge = j.value("q", 0.0); }

void Dipole::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) { mu = q * mu; }

void Dipole::to_json(json &j) const {
    j["mu"] = mu;
    j["mulen"] = mulen;
}

void Dipole::from_json(const json &j) {
    mu = j.value("mu", Point(1, 0, 0));
    mulen = j.value("mulen", mulen);
}

void Polarizable::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) {
    mui = q * mui;
    alpha.rotate(m);
}

void Polarizable::to_json(json &j) const {
    j["alpha"] = alpha;
    j["mui"] = mui;
    j["muilen"] = muilen;
}

void Polarizable::from_json(const json &j) {
    alpha = j.value("alpha", alpha);
    mui = j.value("mui", Point(1, 0, 0));
    muilen = j.value("mulen", muilen);
}

void Quadrupole::rotate(const Eigen::Quaterniond &, const Eigen::Matrix3d &m) { Q.rotate(m); }

void Quadrupole::to_json(json &j) const { j["Q"] = Q; }

void Quadrupole::from_json(const json &j) { Q = j.value("Q", Q); }

void Cigar::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) { scdir = q * scdir; }

void Cigar::to_json(json &j) const {
    j["scdir"] = scdir;
    j["sclen"] = sclen;
}
void Cigar::from_json(const json &j) {
    scdir = j.value("scdir", scdir);
    sclen = j.value("sclen", sclen);
}

const AtomData &Particle::traits() {
    assert(id >= 0 and id < atoms.size());
    return atoms.at(id);
}

/**
 * @warning Performance is sub-optimal as conversion is done through a json object
 */
Particle::Particle(const AtomData &a) { *this = json(a).front(); }

// copy constructor
Particle::Particle(const Particle &p) : id(p.id), charge(p.charge), pos(p.pos) {
    if (p.ext != nullptr)
        ext = std::make_shared<Particle::ParticleExtension>(*p.ext); // deep copy
}

// assignment operator
Particle &Particle::operator=(const Particle &p) {
    if (&p != this) {
        charge = p.charge;
        pos = p.pos;
        id = p.id;
        if (p.ext != nullptr) { // if particle has
            if (ext != nullptr) // extension, then
                *ext = *p.ext;  // deep copy
            else                // else if *this is empty, create new based on p
                ext = std::make_shared<Particle::ParticleExtension>(*p.ext); // create new
        } else                                                               // p doesn't have extended properties
            ext = nullptr;
    }
    return *this;
}

void Particle::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &m) {
    if (ext != nullptr)
        ext->rotate(q, m);
}

bool Particle::hasExtension() const { return ext != nullptr; }

Particle::ParticleExtension &Particle::createExtension() {
    assert(ext == nullptr && "extension already created");
    ext = std::make_shared<ParticleExtension>();
    return *ext;
}

void from_json(const json &j, Particle &p) {
    p.id = j.value("id", -1);
    p.pos = j.value("pos", Point(0, 0, 0));
    p.charge = j.value("q", 0.0);

    p.ext = std::make_shared<Particle::ParticleExtension>();
    from_json(j, *p.ext);
    Particle::ParticleExtension empty_extended_particle;
    // why can't we compare ParticleExtension directly?!
    // (slow and ugly)
    if (json(*p.ext) == json(empty_extended_particle))
        p.ext = nullptr; // no extended features found in json
}
void to_json(json &j, const Particle &p) {
    if (p.ext != nullptr)
        to_json(j, *p.ext);
    j["id"] = p.id;
    j["pos"] = p.pos;
    j["q"] = p.charge;
}

} // namespace Faunus
