#include <doctest/doctest.h>
#include "particle.h"
#include "rotate.h"
#include "units.h"
#include "aux/eigensupport.h"
#include <fstream>
#include <Eigen/Geometry>
#include <cereal/archives/binary.hpp>

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

const AtomData &Particle::traits() const { return atoms.at(id); }

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

TEST_SUITE_BEGIN("Particle");

// convert test to use `Particle::shape`
TEST_CASE("[Faunus] Particle") {
    using doctest::Approx;
    Particle p1, p2;
    p1.id = 100;
    p1.pos = {1, 2, 3};
    p1.charge = -0.8;

    CHECK(p1.hasExtension() == false);

    p1.createExtension();
    CHECK(p1.hasExtension() == true);

    p2.createExtension();
    CHECK(p2.hasExtension() == true);

    p1.getExt().mu = {0, 0, 1};
    p1.getExt().mulen = 2.8;
    p1.getExt().scdir = {-0.1, 0.3, 1.9};
    p1.getExt().sclen = 0.5;
    p1.getExt().Q = Tensor(1, 2, 3, 4, 5, 6);

    p2 = json(p1); // p1 --> json --> p2
    CHECK(p2.hasExtension() == true);

    CHECK(json(p1) == json(p2)); // p1 --> json == json <-- p2 ?

    CHECK(p2.id == 100);
    CHECK(p2.pos == Point(1, 2, 3));
    CHECK(p2.charge == -0.8);
    CHECK(p2.getExt().mu == Point(0, 0, 1));
    CHECK(p2.getExt().mulen == 2.8);
    CHECK(p2.getExt().scdir == Point(-0.1, 0.3, 1.9));
    CHECK(p2.getExt().sclen == 0.5);
    CHECK(p2.getExt().Q == Tensor(1, 2, 3, 4, 5, 6));

    // check of all properties are rotated
    QuaternionRotate qrot(pc::pi / 2, {0, 1, 0});
    p1.getExt().mu = p1.getExt().scdir = {1, 0, 0};
    p1.rotate(qrot.getQuaternion(), qrot.getRotationMatrix());

    CHECK(p1.getExt().mu.x() == Approx(0));
    CHECK(p1.getExt().mu.z() == Approx(-1));
    CHECK(p1.getExt().scdir.x() == Approx(0));
    CHECK(p1.getExt().scdir.z() == Approx(-1));

    CHECK(p1.getExt().Q(0, 0) == Approx(6));
    CHECK(p1.getExt().Q(0, 1) == Approx(5));
    CHECK(p1.getExt().Q(0, 2) == Approx(-3));
    CHECK(p1.getExt().Q(1, 1) == Approx(4));
    CHECK(p1.getExt().Q(1, 2) == Approx(-2));
    CHECK(p1.getExt().Q(2, 2) == Approx(1));

    SUBCASE("Cereal serialisation") {
        Particle p;
        p.pos = {10, 20, 30};
        p.charge = -1;
        p.id = 8;
        p.ext = std::make_shared<Particle::ParticleExtension>();
        p.getExt().mu = {0.1, 0.2, 0.3};
        p.getExt().mulen = 104;

        { // write
            std::ofstream os("out.cereal", std::ios::binary);
            cereal::BinaryOutputArchive archive(os);
            archive(p);
        }

        // set to zero
        p.pos = {0, 0, 0};
        p.charge = 0;
        p.id = 0;
        p.getExt().mu.setZero();
        p.getExt().mulen = 0;

        { // read
            std::ifstream is("out.cereal", std::ios::binary);
            cereal::BinaryInputArchive archive(is);
            archive(p);
            CHECK(p.id == 8);
            CHECK(p.charge == -1);
            CHECK(p.pos.x() == 10);
            CHECK(p.pos.y() == 20);
            CHECK(p.pos.z() == 30);
            CHECK(p.getExt().mu.x() == 0.1);
            CHECK(p.getExt().mu.y() == 0.2);
            CHECK(p.getExt().mu.z() == 0.3);
            CHECK(p.getExt().mulen == 104);
        }
    }
}

TEST_SUITE_END();

} // namespace Faunus
