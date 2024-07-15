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
void ParticlePropertyBase::rotate(const Eigen::Quaterniond&, const Eigen::Matrix3d&) {}

void Radius::to_json(json& j) const
{
    j["r"] = radius;
}

void Radius::from_json(const json& j)
{
    radius = j.value("r", 0.0);
}

void Charge::to_json(json& j) const
{
    j["q"] = charge;
}

void Charge::from_json(const json& j)
{
    charge = j.value("q", 0.0);
}

void Dipole::rotate(const Eigen::Quaterniond& q, const Eigen::Matrix3d&)
{
    mu = q * mu;
}

void Dipole::to_json(json& j) const
{
    j["mu"] = mu;
    j["mulen"] = mulen;
}

void Dipole::from_json(const json& j)
{
    mu = j.value("mu", Point::Zero().eval());
    mulen = j.value("mulen", 0.0);
}

bool Dipole::isDipolar() const
{
    return mulen > pc::epsilon_dbl;
}

void Polarizable::rotate(const Eigen::Quaterniond& q, const Eigen::Matrix3d& m)
{
    mui = q * mui;
    alpha.rotate(m);
}

void Polarizable::to_json(json& j) const
{
    j["alpha"] = alpha;
    j["mui"] = mui;
    j["muilen"] = muilen;
}

void Polarizable::from_json(const json& j)
{
    alpha = j.value("alpha", alpha);
    mui = j.value("mui", Point(1, 0, 0));
    muilen = j.value("mulen", muilen);
}

bool Polarizable::isPolarizable() const
{
    return alpha.cwiseAbs().sum() > pc::epsilon_dbl;
}

void Quadrupole::rotate(const Eigen::Quaterniond&, const Eigen::Matrix3d& m)
{
    Q.rotate(m);
}

void Quadrupole::to_json(json& j) const
{
    j["Q"] = Q;
}

void Quadrupole::from_json(const json& j)
{
    Q = j.value("Q", Tensor(Tensor::Zero()));
}

bool Quadrupole::isQuadrupolar() const
{
    return Q.cwiseAbs().sum() > pc::epsilon_dbl;
}

void Cigar::rotate(const Eigen::Quaterniond& quaternion,
                   [[maybe_unused]] const Eigen::Matrix3d& rotation_matrix)
{
    scdir = quaternion * scdir;
    patchdir = quaternion * patchdir;
    patchsides.at(0) = quaternion * patchsides.at(0);
    patchsides.at(1) = quaternion * patchsides.at(1);
    // here we may want to run `initialize()` to reduce numerical rounding issues
    // after many rotations...will require access to AtomData.
}

void Cigar::to_json(json& j) const
{
    j["scdir"] = scdir;
    j["pdir"] = patchdir;
}

void Cigar::from_json(const json& j)
{
    assert(j.contains("id"));
    const auto& psc = Faunus::atoms.at(j.at("id").get<int>()).sphero_cylinder;
    setDirections(psc, j.value("scdir", Point(1.0, 0.0, 0.0)),
                  j.value("pdir", Point(0.0, 1.0, 0.0)));
}

/**
 * Calculates cosine of angles, patch direction including chirality
 * and vector corresponding to sides of patch that are used in
 * calculations of interactions.
 * This function must be called at the beginning of calculations and after changes
 * of patch properties.
 * It shall be also after a lot of moves to remove accumulated errors
 */
void Cigar::setDirections(const SpheroCylinderData& psc_data, const Point& new_direction,
                          const Point& new_patch_direction)
{
    constexpr auto very_small_number = 1e-9;
    half_length = 0.5 * psc_data.length;
    if (half_length < very_small_number) {
        return;
    }

    pcangl = std::cos(0.5 * psc_data.patch_angle);
    pcanglsw = std::cos(0.5 * psc_data.patch_angle + psc_data.patch_angle_switch);

    scdir = scdir.squaredNorm() < very_small_number ? Point(1, 0, 0) : new_direction.normalized();
    patchdir = patchdir.squaredNorm() < very_small_number ? Point(0, 1, 0)
                                                          : new_patch_direction.normalized();
    if (fabs(patchdir.dot(scdir)) > very_small_number) { // must be perpendicular
        faunus_logger->trace("straigthening patch direction");
        patchdir = (patchdir - scdir * patchdir.dot(scdir)).normalized(); // perp. project
    }

    /* calculate patch sides */
    Point patch_length_axis = scdir;
    Eigen::Quaterniond quaternion;
    if (psc_data.chiral_angle > very_small_number) {
        quaternion = Eigen::AngleAxisd(0.5 * psc_data.chiral_angle, patchdir);
        patch_length_axis = quaternion * patch_length_axis;
    }

    /* create side vector by rotating patch vector by half size of patch*/
    const auto half_angle = 0.5 * psc_data.patch_angle + psc_data.patch_angle_switch;
    quaternion = Eigen::AngleAxisd(half_angle, patch_length_axis);
    patchsides.at(0) = (quaternion * patchdir).normalized();
    quaternion = Eigen::AngleAxisd(-half_angle, patch_length_axis);
    patchsides.at(1) = (quaternion * patchdir).normalized();

    if (patchsides.at(0).squaredNorm() < very_small_number) {
        throw std::runtime_error("patch side vector has zero size");
    }
}

bool Cigar::isCylindrical() const
{
    return half_length > pc::epsilon_dbl;
}

const AtomData& Particle::traits() const
{
    return atoms[id];
}

/**
 * @warning Performance is sub-optimal as conversion is done through a json object
 */
Particle::Particle(const AtomData& atomdata)
{
    *this = static_cast<json>(atomdata).front();
}

Particle::Particle(const AtomData& atomdata, const Point& pos)
    : Particle(atomdata)
{
    this->pos = pos;
}

Particle::Particle(const Particle& other)
    : id(other.id)
    , charge(other.charge)
    , pos(other.pos)
{
    if (other.hasExtension()) {
        ext = std::make_unique<Particle::ParticleExtension>(other.getExt()); // deep copy
    }
}

Particle& Particle::operator=(const Particle& other)
{
    if (&other != this) {
        charge = other.charge;
        pos = other.pos;
        id = other.id;
        if (other.hasExtension()) {    // if particle has
            if (ext) {                 // extension, then
                *ext = other.getExt(); // deep copy
            }
            else { // else if *this is empty, create new based on p
                ext = std::make_unique<Particle::ParticleExtension>(
                    other.getExt()); // create and deep copy
            }
        }
        else { // other doesn't have extended properties
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
void Particle::rotate(const Eigen::Quaterniond& quaternion, const Eigen::Matrix3d& rotation_matrix)
{
    if (hasExtension()) {
        ext->rotate(quaternion, rotation_matrix);
    }
}

bool Particle::hasExtension() const
{
    return ext != nullptr;
}

Particle::ParticleExtension& Particle::createExtension()
{
    if (!ext) {
        ext = std::make_unique<ParticleExtension>();
    }
    return *ext;
}

void from_json(const json& j, Particle& particle)
{
    assert(j.contains("id"));
    particle.id = j.at("id").get<int>();
    particle.pos = j.value("pos", Point::Zero().eval());
    particle.charge = j.value("q", 0.0);
    if (j.contains("psc") || j.contains("mu") || j.contains("Q") || j.contains("scdir")) {
        particle.ext = std::make_unique<Particle::ParticleExtension>(j);
        // delete extension if not in use:
        if (!particle.ext->isCylindrical() && !particle.ext->isDipolar() &&
            !particle.ext->isQuadrupolar() && !particle.ext->isQuadrupolar()) {
            particle.ext = nullptr;
        }
    }
}

void to_json(json& j, const Particle& particle)
{
    if (particle.hasExtension()) {
        to_json(j, particle.getExt());
    }
    j["id"] = particle.id;
    j["pos"] = particle.pos;
    j["q"] = particle.charge;
}

TEST_SUITE_BEGIN("Particle");

TEST_CASE("[Faunus] Particle")
{
    using doctest::Approx;

    json j = R"({ "atomlist" : [
             { "A": { "sigma": 2.5, "pactivity":2, "eps_custom": 0.1 } },
             { "B": { "psc": {"length": 9, "chiral_angle": 5.0, "patch_angle": 30} } }
             ]})"_json;

    pc::temperature = 298.15_K;
    atoms = j["atomlist"].get<decltype(atoms)>();

    Particle p1, p2;
    p1.id = 1;
    p1.pos = {1, 2, 3};
    p1.charge = -0.8;

    CHECK_EQ(p1.hasExtension(), false);
    CHECK_EQ(p2.hasExtension(), false);

    p1.createExtension();
    [[maybe_unused]] auto& newly_created_extension = p2.getExt();
    CHECK_EQ(p1.hasExtension(), true);
    CHECK_EQ(p2.hasExtension(), true);

    p1.getExt().mu = {0, 0, 1};
    p1.getExt().mulen = 2.8;

    p1.getExt().scdir = Point(-0.1, 0.3, 1.9).normalized();
    p1.getExt().patchdir = Point::Ones().normalized();
    p1.getExt().patchdir =
        (p1.getExt().patchdir - p1.getExt().scdir * p1.getExt().patchdir.dot(p1.getExt().scdir))
            .normalized();
    CHECK(p1.getExt().patchdir.dot(p1.getExt().scdir) < 1e-9); // must be perpendicular

    p1.getExt().Q = Tensor(1, 2, 3, 4, 5, 6);

    p2 = json(p1); // p1 --> json --> p2
    CHECK_EQ(p2.hasExtension(), true);

    CHECK_EQ(json(p1), json(p2)); // p1 --> json == json <-- p2 ?

    CHECK_EQ(p2.id, 1);
    CHECK_EQ(p2.pos, Point(1, 2, 3));
    CHECK_EQ(p2.charge, -0.8);
    CHECK_EQ(p2.getExt().mu, Point(0, 0, 1));
    CHECK_EQ(p2.getExt().mulen, 2.8);
    CHECK_EQ(p2.getExt().Q, Tensor(1, 2, 3, 4, 5, 6));

    SUBCASE("Cigar")
    {
        const auto& cigar = p2.getExt();
        CHECK_EQ(cigar.scdir, Point(-0.1, 0.3, 1.9).normalized());
        CHECK_EQ(cigar.scdir.x(), Approx(-0.0519174));
        CHECK_EQ(cigar.scdir.y(), Approx(0.155752));
        CHECK_EQ(cigar.scdir.z(), Approx(0.986431));
        CHECK(cigar.scdir.dot(p2.getExt().patchdir) < 1e-9);
        CHECK_EQ(cigar.patchdir.x(), Approx(0.7850810157));
        CHECK_EQ(cigar.patchdir.y(), Approx(0.6168493695));
        CHECK_EQ(cigar.patchdir.z(), Approx(-0.0560772154));
        CHECK_EQ(cigar.half_length, Approx(4.5));
        CHECK_EQ(cigar.patchsides.at(0).x(), Approx(0.5981493663));
        CHECK_EQ(cigar.patchsides.at(0).y(), Approx(0.7970822805));
        CHECK_EQ(cigar.patchsides.at(0).z(), Approx(-0.0829287266));
        CHECK_EQ(cigar.patchsides.at(1).x(), Approx(0.9185106913));
        CHECK_EQ(cigar.patchsides.at(1).y(), Approx(0.3945791934));
        CHECK_EQ(cigar.patchsides.at(1).z(), Approx(-0.0254041347));
        CHECK_EQ(p2.traits().sphero_cylinder.patch_angle, 30.0 * 1.0_deg);
        CHECK_EQ(p2.traits().sphero_cylinder.chiral_angle, 5.0 * 1.0_deg);
    }

    SUBCASE("rotate")
    {
        // check if all properties are rotated
        QuaternionRotate qrot(pc::pi / 2, {0, 1, 0});
        p1.getExt().mu = p1.getExt().scdir = {1, 0, 0};
        p1.rotate(qrot.getQuaternion(), qrot.getRotationMatrix());

        CHECK_EQ(p1.getExt().mu.x(), Approx(0));
        CHECK_EQ(p1.getExt().mu.z(), Approx(-1));
        CHECK_EQ(p1.getExt().scdir.x(), Approx(0));
        CHECK_EQ(p1.getExt().scdir.z(), Approx(-1));

        CHECK_EQ(p1.getExt().Q(0, 0), Approx(6));
        CHECK_EQ(p1.getExt().Q(0, 1), Approx(5));
        CHECK_EQ(p1.getExt().Q(0, 2), Approx(-3));
        CHECK_EQ(p1.getExt().Q(1, 1), Approx(4));
        CHECK_EQ(p1.getExt().Q(1, 2), Approx(-2));
        CHECK_EQ(p1.getExt().Q(2, 2), Approx(1));
    }

    SUBCASE("Cereal serialisation")
    {
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
            CHECK_EQ(p.id, 8);
            CHECK_EQ(p.charge, -1);
            CHECK_EQ(p.pos.x(), 10);
            CHECK_EQ(p.pos.y(), 20);
            CHECK_EQ(p.pos.z(), 30);
            CHECK_EQ(p.getExt().mu.x(), 0.1);
            CHECK_EQ(p.getExt().mu.y(), 0.2);
            CHECK_EQ(p.getExt().mu.z(), 0.3);
            CHECK_EQ(p.getExt().mulen, 104);
        }
    }
}

TEST_SUITE_END();

} // namespace Faunus
