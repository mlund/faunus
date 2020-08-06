#include "particle.h"
#include <Eigen/Geometry>
#include "aux/eigensupport.h"

namespace Faunus {

/**
 * @todo Inefficient to pass both Quaternion and rotation matrix. This is currently
 * done since tensors are rotated using the latter. Could these be rotated using
 * quaternions?
 */
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

const AtomData &Particle::traits() { return atoms.at(id); }

/**
 * @warning Performance is sub-optimal as conversion is done through a json object
 */
Particle::Particle(const AtomData &a) { *this = json(a).front(); }
Particle::Particle(const AtomData &a, const Point &pos) : Particle(a) { this->pos = pos; }
Particle::Particle(const Particle &other) : id(other.id), charge(other.charge), pos(other.pos) {
    if (other.ext) {
        ext = std::make_shared<Particle::ParticleExtension>(*other.ext); // deep copy
    }
}

Particle &Particle::operator=(const Particle &other) {
    if (&other != this) {
        charge = other.charge;
        pos = other.pos;
        id = other.id;
        if (other.ext) {           // if particle has
            if (ext) {             // extension, then
                *ext = *other.ext; // deep copy
            } else {               // else if *this is empty, create new based on p
                ext = std::make_shared<Particle::ParticleExtension>(*other.ext); // create new
            }
        } else { // other doesn't have extended properties
            ext = nullptr;
        }
    }
    return *this;
}

/**
 * @param quaternion Quaternion used to rotate points
 * @param rotation_matrix Rotation matrix used to rotate tensors
 * @todo Only one of the above should be enough
 *
 * This has effect only on anisotropic particles and does no
 * positional rotation, only internal rotation
 */
void Particle::rotate(const Eigen::Quaterniond &quaternion, const Eigen::Matrix3d &rotation_matrix) {
    if (hasExtension()) {
        ext->rotate(quaternion, rotation_matrix);
    }
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
    if (p.ext) {
        to_json(j, *p.ext);
    }
    j["id"] = p.id;
    j["pos"] = p.pos;
    j["q"] = p.charge;
}

} // namespace Faunus
