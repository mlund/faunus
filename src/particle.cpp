#include <doctest/doctest.h>
#include "particle.h"
#include "rotate.h"
#include "units.h"
#include "aux/eigensupport.h"
#include <fstream>
#include <Eigen/Geometry>
#include <cereal/types/memory.hpp>
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
    mu = j.value("mu", Point::Zero().eval());
    mulen = j.value("mulen", 0.0);
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

void Quadrupole::from_json(const json& j) { Q = j.value("Q", Tensor(Tensor::Zero()));}

void Cigar::rotate(const Eigen::Quaterniond &q, const Eigen::Matrix3d &) {
    scdir = q * scdir;
    patchdir = q * patchdir;
    chirality_direction = q * patchdir;
    patchsides.at(0) = q * patchsides.at(0);
    patchsides.at(1) = q * patchsides.at(1);
    // here we may want to run `initialize()` to reduce numerical rounding issues
    // after many rotations...will require access to AtomData.
}

void Cigar::to_json(json& j) const { j["scdir"] = scdir; }

void Cigar::from_json(const json& j) {
    const auto& psc = Faunus::atoms.at(j.at("id").get<int>()).sphero_cylinder;
    scdir = j.value("scdir", Point(1.0, 0.0, 0.0));
    initialize(psc);
}

/**
 * Calculates cosine of angles, patch direction including chirality
 * and vector corresponding to sides of patch that are used in
 * calculations of interactions.
 * This function must be called at the beginning of calculations and after changes
 * of patch properties.
 * It shall be also after a lot of moves to remove accumulated errors
 *
 * @note Largely from Robert Vacha's C code
 */
void Cigar::initialize(const SpheroCylinderData& psc) {
    constexpr auto very_small_number = 1e-9;
    half_length = 0.5 * psc.length;
    if (half_length > very_small_number) {
        Point vec;
        Eigen::Quaterniond Q;
        pcangl = std::cos(0.5 * psc.patch_angle);
        pcanglsw = std::cos(0.5 * psc.patch_angle + psc.patch_angle_switch);

        if (scdir.squaredNorm() < very_small_number) {
            scdir = {1, 0, 0};
        }
        if (patchdir.squaredNorm() < very_small_number) {
            patchdir = {0, 1, 0};
        }
        scdir.normalize();

        patchdir = patchdir - scdir * patchdir.dot(scdir); // perp. project
        patchdir.normalize();

        /* calculate patch sides */
        if (psc.chiral_angle < very_small_number) {
            vec = scdir;
        } else {
            chirality_direction = scdir;
            Q = Eigen::AngleAxisd(0.5 * psc.chiral_angle, patchdir);
            chirality_direction = Q * chirality_direction; // rotate
            vec = chirality_direction;
        }

        /* create side vector by rotating patch vector by half size of patch*/
        /* the first side */
        patchsides[0] = patchdir;
        Q = Eigen::AngleAxisd(0.5 * psc.patch_angle + psc.patch_angle_switch, vec);
        patchsides[0] = Q * patchsides[0]; // rotate
        patchsides[0].normalize();

        /* the second side */
        patchsides[1] = patchdir;
        Q = Eigen::AngleAxisd(-0.5 * psc.patch_angle - psc.patch_angle_switch, vec);
        patchsides[1] = Q * patchsides[1]; // rotate
        patchsides[1].normalize();

        if (patchsides[0].squaredNorm() < very_small_number) {
            throw std::runtime_error("Patch side vector has zero size.");
        }
    }
}

const AtomData &Particle::traits() const { return atoms[id]; }

/**
 * @warning Performance is sub-optimal as conversion is done through a json object
 */
Particle::Particle(const AtomData &a) {*this = json(a).front(); }
Particle::Particle(const AtomData &a, const Point &pos) : Particle(a) { this->pos = pos; }
Particle::Particle(const Particle &other) : id(other.id), charge(other.charge), pos(other.pos) {
    if (other.hasExtension()) {
        ext = std::make_unique<Particle::ParticleExtension>(other.getExt()); // deep copy
    }
}

Particle &Particle::operator=(const Particle &other) {
    if (&other != this) {
        charge = other.charge;
        pos = other.pos;
        id = other.id;
        if (other.hasExtension()) {    // if particle has
            if (ext) {                 // extension, then
                *ext = other.getExt(); // deep copy
            } else {                   // else if *this is empty, create new based on p
                ext = std::make_unique<Particle::ParticleExtension>(other.getExt()); // create and deep copy
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

Particle::ParticleExtension& Particle::createExtension() {
    if (!ext) {
        ext = std::make_unique<ParticleExtension>();
    }
    return *ext;
}

void from_json(const json &j, Particle &particle) {
    particle.id = j.value("id", -1);
    particle.pos = j.value("pos", Point::Zero().eval());
    particle.charge = j.value("q", 0.0);
    particle.ext = std::make_unique<Particle::ParticleExtension>(j);
    // Disable extended features if unused. Slow and ugly check:
    const auto empty_extended_particle = json(Particle::ParticleExtension());
    if (json(*particle.ext) == empty_extended_particle) {
        particle.ext = nullptr; // no extended features found in json
    }
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

TEST_CASE("[Faunus] Particle") {
    using doctest::Approx;

    json j = R"({ "atomlist" : [
             { "A": { "sigma": 2.5, "pactivity":2, "eps_custom": 0.1 } },
             { "B": { "psc": {"length": 9, "chiral_angle": 5.0} } }
             ]})"_json;

    pc::temperature = 298.15_K;
    atoms = j["atomlist"].get<decltype(atoms)>();

    Particle p1, p2;
    p1.id = 1;
    p1.pos = {1, 2, 3};
    p1.charge = -0.8;

    CHECK(p1.hasExtension() == false);
    CHECK(p2.hasExtension() == false);

    p1.createExtension();
    [[maybe_unused]] auto& newly_created_extension = p2.getExt();
    CHECK(p1.hasExtension() == true);
    CHECK(p2.hasExtension() == true);

    p1.getExt().mu = {0, 0, 1};
    p1.getExt().mulen = 2.8;
    p1.getExt().scdir = Point(-0.1, 0.3, 1.9).normalized();
    p1.getExt().Q = Tensor(1, 2, 3, 4, 5, 6);

    p2 = json(p1); // p1 --> json --> p2
    CHECK(p2.hasExtension() == true);

    CHECK(json(p1) == json(p2)); // p1 --> json == json <-- p2 ?

    CHECK(p2.id == 1);
    CHECK(p2.pos == Point(1, 2, 3));
    CHECK(p2.charge == -0.8);
    CHECK(p2.getExt().mu == Point(0, 0, 1));
    CHECK(p2.getExt().mulen == 2.8);
    CHECK(p2.getExt().Q == Tensor(1, 2, 3, 4, 5, 6));

    SUBCASE("Cigar") {
        CHECK(p2.getExt().scdir == Point(-0.1, 0.3, 1.9).normalized());
        CHECK(p2.getExt().scdir.x() == Approx(-0.0519174));
        CHECK(p2.getExt().scdir.y() == Approx(0.155752));
        CHECK(p2.getExt().scdir.z() == Approx(0.986431));
        CHECK(p2.getExt().chirality_direction.x() == Approx(-0.0083089));
        CHECK(p2.getExt().chirality_direction.y() == Approx(0.155604));
        CHECK(p2.getExt().chirality_direction.z() == Approx(0.987785));
        CHECK(p2.getExt().patchdir.x() == Approx(0.008186156));
        CHECK(p2.getExt().patchdir.y() == Approx(0.987796153));
        CHECK(p2.getExt().patchdir.z() == Approx(-0.1555369633));
        CHECK(p2.getExt().half_length == 4.5);
    }

    SUBCASE("rotate") {
        // check if all properties are rotated
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
    }

    SUBCASE("Cereal serialisation") {
        Particle p;
        p.pos = {10, 20, 30};
        p.charge = -1;
        p.id = 8;
        p.ext = std::make_unique<Particle::ParticleExtension>();
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
