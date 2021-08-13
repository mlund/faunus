#include "geometry.h"
#include "atomdata.h"
#include "random.h"
#include "average.h"
#include "aux/eigensupport.h"
#include <spdlog/spdlog.h>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>
#include <Eigen/Eigenvalues>

namespace Faunus::Geometry {

GeometryBase::~GeometryBase() = default;

/**
 * @return Boolean matrix with true for each periodic direction
 */
Eigen::Matrix<bool, 3, 1> BoundaryCondition::isPeriodic() const {
    return direction.cwiseEqual(
        BoundaryXYZ(Boundary::PERIODIC, Boundary::PERIODIC, Boundary::PERIODIC));
}

GeometryImplementation::~GeometryImplementation() = default;

Cuboid::Cuboid(const Point &side_length) {
    boundary_conditions =
        BoundaryCondition(Coordinates::ORTHOGONAL, {Boundary::PERIODIC, Boundary::PERIODIC, Boundary::PERIODIC});
    setLength(side_length);
}

Cuboid::Cuboid() : Cuboid(Point::Zero()) {}

Point Cuboid::getLength() const { return box; }

void Cuboid::setLength(const Point &len) {
    box = len;
    box_half = 0.5 * box;
    box_inv = box.cwiseInverse();
}

double Cuboid::getVolume(int) const { return box.x() * box.y() * box.z(); }

Point Cuboid::setVolume(double volume, const VolumeMethod method) {
    const double old_volume = getVolume();
    double alpha;
    Point new_box;
    Point box_scaling;
    switch (method) {
    case VolumeMethod::ISOTROPIC:
        alpha = std::cbrt(volume / old_volume);
        box_scaling = {alpha, alpha, alpha};
        break;
    case VolumeMethod::XY:
        alpha = std::sqrt(volume / old_volume);
        box_scaling = {alpha, alpha, 1.0};
        break;
    case VolumeMethod::Z:
        alpha = volume / old_volume;
        box_scaling = {1.0, 1.0, alpha};
        break;
    case VolumeMethod::ISOCHORIC:
        // z is scaled by 1/alpha/alpha, x and y are scaled by alpha
        alpha = std::cbrt(volume / old_volume);
        box_scaling = {alpha, alpha, 1 / (alpha * alpha)};
        volume = old_volume;
        break;
    default:
        throw std::invalid_argument("unsupported volume scaling method for the cuboid geometry");
    }
    new_box = box.cwiseProduct(box_scaling);
    setLength(new_box);
    assert(fabs(getVolume() - volume) < 1e-6);
    return box_scaling; // this will scale any point to new volume
}

void Cuboid::boundary(Point &a) const {
    if (boundary_conditions.direction.x() == Boundary::PERIODIC) {
        if (std::fabs(a.x()) > box_half.x())
            a.x() -= box.x() * anint(a.x() * box_inv.x());
    }
    if (boundary_conditions.direction.y() == Boundary::PERIODIC) {
        if (std::fabs(a.y()) > box_half.y())
            a.y() -= box.y() * anint(a.y() * box_inv.y());
    }
    if (boundary_conditions.direction.z() == Boundary::PERIODIC) {
        if (std::fabs(a.z()) > box_half.z())
            a.z() -= box.z() * anint(a.z() * box_inv.z());
    }
}

Point Cuboid::vdist(const Point &a, const Point &b) const {
    Point distance(a - b);
    if (boundary_conditions.direction.x() == Boundary::PERIODIC) {
        if (distance.x() > box_half.x())
            distance.x() -= box.x();
        else if (distance.x() < -box_half.x())
            distance.x() += box.x();
    }
    if (boundary_conditions.direction.y() == Boundary::PERIODIC) {
        if (distance.y() > box_half.y())
            distance.y() -= box.y();
        else if (distance.y() < -box_half.y())
            distance.y() += box.y();
    }
    if (boundary_conditions.direction.z() == Boundary::PERIODIC) {
        if (distance.z() > box_half.z())
            distance.z() -= box.z();
        else if (distance.z() < -box_half.z())
            distance.z() += box.z();
    }
    return distance;
}

void Cuboid::randompos(Point &m, Random &rand) const {
    m.x() = (rand() - 0.5) * box.x();
    m.y() = (rand() - 0.5) * box.y();
    m.z() = (rand() - 0.5) * box.z();
}

bool Cuboid::collision(const Point &a) const {
    return (std::fabs(a.x()) > box_half.x() || std::fabs(a.y()) > box_half.y() || std::fabs(a.z()) > box_half.z());
}

void Cuboid::from_json(const json &j) {
    box.setZero();
    if (auto length_ = j.at("length"); length_.is_number()) {
        auto l = length_.get<double>();
        setLength({l, l, l});
        return;
    } else if (length_.is_array()) {
        if (length_.size() == 3) {
            setLength(length_.get<Point>());
            return;
        }
    }
    throw ConfigurationError("sidelength syntax error");
}

void Cuboid::to_json(json &j) const { j = {{"length", box}}; }

// =============== Slit ===============

Slit::Slit(const Point &p) : Tbase(p) {
    boundary_conditions =
        BoundaryCondition(Coordinates::ORTHOGONAL, {Boundary::PERIODIC, Boundary::PERIODIC, Boundary::FIXED});
}

Slit::Slit(double x, double y, double z) : Slit(Point(x, y, z)) {}
Slit::Slit(double x) : Slit(x, x, x) {}

// =============== Sphere ===============

Sphere::Sphere(double radius) : radius(radius) {
    boundary_conditions =
        BoundaryCondition(Coordinates::ORTHOGONAL, {Boundary::FIXED, Boundary::FIXED, Boundary::FIXED});
}

Point Sphere::getLength() const { return {2.0 * radius, 2.0 * radius, 2.0 * radius}; }

double Sphere::getVolume(int dim) const {
    double result;
    switch (dim) {
    case 3:
        result = 4.0 / 3.0 * pc::pi * radius * radius * radius; // volume
        break;
    case 2:
        result = 4.0 * pc::pi * radius * radius; // surface area
        break;
    case 1:
        result = 2.0 * radius; // diameter
        break;
    default:
        throw std::invalid_argument("unsupported volume dimension for the sphere: " + std::to_string(dim));
    }
    return result;
}

Point Sphere::setVolume(double volume, const VolumeMethod method) {
    const double old_radius = radius;
    Point box_scaling;
    if (method == VolumeMethod::ISOTROPIC) {
        radius = std::cbrt(volume / (4.0 / 3.0 * pc::pi));
        assert(std::fabs(getVolume() - volume) < 1e-6 && "error setting sphere volume");
    } else {
        throw std::invalid_argument("unsupported volume scaling method for the spherical geometry");
    }
    box_scaling.setConstant(radius / old_radius);
    return box_scaling;
}

void Sphere::boundary(Point &) const {} // no PBC

bool Sphere::collision(const Point &point) const { return point.squaredNorm() > radius * radius; }

Point Sphere::vdist(const Point &a, const Point &b) const {
    // no pbc; shall we check points coordinates?
    Point distance(a - b);
    return distance;
}

void Sphere::randompos(Point &m, Random &rand) const {
    double r2 = radius * radius, d = 2 * radius;
    do {
        m.x() = (rand() - 0.5) * d;
        m.y() = (rand() - 0.5) * d;
        m.z() = (rand() - 0.5) * d;
    } while (m.squaredNorm() > r2);
}

void Sphere::from_json(const json &j) { radius = j.at("radius").get<double>(); }

void Sphere::to_json(json &j) const { j = {{"radius", radius}}; }

double Sphere::getRadius() const { return radius; }

// =============== Hypersphere 2D ===============

Hypersphere2d::Hypersphere2d(double radius) : Sphere(radius) {
    boundary_conditions = BoundaryCondition(Coordinates::NON3D);
}

Point Hypersphere2d::vdist(const Point &a, const Point &b) const {
    // ugly but works, needs fixing though...
    Point distance3d(a - b);
    double angle = std::acos(a.dot(b) / (radius * radius));
    double distance = radius * angle;
    return distance3d / distance3d.norm() * distance;
}

void Hypersphere2d::randompos(Point &m, Random &rand) const {
    Sphere::randompos(m, rand);
    m = m / m.norm() * radius;
}

bool Hypersphere2d::collision(const Point &a) const {
    bool collision = std::fabs(a.norm() - radius) > 1e-6;
    return collision;
}

// =============== Hexagonal Prism ===============

const Eigen::Matrix3d HexagonalPrism::rhombic2cartesian =
    (Eigen::Matrix3d() << std::cos(-pc::pi / 6.0), 0.0, 0.0, std::sin(-pc::pi / 6.0), 1.0, 0.0, 0.0, 0.0, 1.0)
        .finished();

const Eigen::Matrix3d HexagonalPrism::cartesian2rhombic = rhombic2cartesian.inverse();

HexagonalPrism::HexagonalPrism(double side, double height) {
    // the current implementation is hardcoded as bellow and ignores other periodic_directions settings
    boundary_conditions =
        BoundaryCondition(Coordinates::ORTHOHEXAGONAL, {Boundary::PERIODIC, Boundary::PERIODIC, Boundary::FIXED});
    set_box(side, height);
}

Point HexagonalPrism::getLength() const { return box; }

double HexagonalPrism::getVolume(int) const {
    return 3.0 / 4.0 * box.x() * box.y() * box.z(); // 3 * inner_radius * outer_radius * height
}

void HexagonalPrism::set_box(double side, double height) { box = {std::sqrt(3.0) * side, 2.0 * side, height}; }

Point HexagonalPrism::setVolume(double volume, const VolumeMethod method) {
    const double old_volume = getVolume();
    double alpha;
    Point box_scaling;

    switch (method) {
    case VolumeMethod::ISOTROPIC:
        alpha = std::cbrt(volume / old_volume);
        box_scaling = {alpha, alpha, alpha};
        break;
    case VolumeMethod::XY:
        alpha = std::sqrt(volume / old_volume);
        box_scaling = {alpha, alpha, 1.0};
        break;
    case VolumeMethod::Z:
        alpha = volume / old_volume;
        box_scaling = {1.0, 1.0, alpha};
        break;
    case VolumeMethod::ISOCHORIC:
        // radius is scaled by alpha, z is scaled by 1/alpha/alpha
        alpha = std::cbrt(volume / old_volume);
        box_scaling = {alpha, alpha, 1.0 / (alpha * alpha)};
        volume = old_volume;
        break;
    default:
        throw std::invalid_argument("unsupported volume scaling method for the hexagonal-prism geometry");
    }
    box = box.cwiseProduct(box_scaling);
    assert(fabs(getVolume() - volume) < 1e-6);
    return box_scaling;
}

Point HexagonalPrism::vdist(const Point &a, const Point &b) const {
    Point distance(a - b);
    boundary(distance);
    return distance;
}

bool HexagonalPrism::collision(const Point &a) const {
    const double height = box.z();
    const double outer_radius = 0.5 * box.y();

    // Hexagon can be divided into three rhombuses. Using natural rhombic coordinates, it can be
    // straightforwardly determined if the particular point is inside the rhombus. If the point is mirrored to
    // the first quadrant taking the absolute values of coordinates, only one rhombus needs to be evaluated.
    Point b = cartesian2rhombic * a.cwiseAbs();
    bool collision = b.z() > 0.5 * height || b.x() > outer_radius || b.y() > outer_radius;
    return collision;
}

void HexagonalPrism::boundary(Point &a) const {
    // TODO optimise and add documentation
    const double sqrtThreeByTwo = sqrt(3.0) / 2.0;
    const Point unitvX = {1.0, 0.0, 0.0};
    const Point unitvY = {0.5, sqrtThreeByTwo, 0.0};
    const Point unitvZ = {-0.5, sqrtThreeByTwo, 0.0};

    double tmp = a.dot(unitvX);
    if (std::fabs(tmp) > 0.5 * box.x())
        a -= box.x() * anint(tmp / box.x()) * unitvX;

    if (a.dot(unitvY) > 0.5 * box.x()) {
        a -= box.x() * unitvY;
        if (a.dot(unitvX) < -0.5 * box.x()) // Check that point did not get past x-limit
            a += box.x() * unitvX;
    }
    if (a.dot(unitvY) < -0.5 * box.x()) {
        a = a + box.x() * unitvY;
        if (a.dot(unitvX) > 0.5 * box.x()) // Check that point did not get past x-limit
            a = a - box.x() * unitvX;
    }

    tmp = a.dot(unitvZ);
    if (std::fabs(tmp) > 0.5 * box.x())
        a -= box.x() * anint(tmp / box.x()) * unitvZ;
    if (std::fabs(a.z()) > 0.5 * box.z())
        a.z() -= box.z() * anint(a.z() / box.z());
}

void HexagonalPrism::randompos(Point &m, Random &rand) const {
    // Generating random points in hexagonal-prism coordinates is not feasible as the space coverage
    // will not be homogeneous back in the cartesian coordinates.
    m.z() = (rand() - 0.5) * box.z();
    do {
        // x,y shall be always generated as a pair; 75% chance to hit
        m.x() = (rand() - 0.5) * box.x();
        m.y() = (rand() - 0.5) * box.y();
    } while (collision(m));
}

void HexagonalPrism::from_json(const json &j) {
    auto radius = j.at("radius").get<double>(); // inner radius
    auto height = j.at("length").get<double>();
    auto edge = 2.0 / std::sqrt(3.0) * radius;
    set_box(edge, height);
}

void HexagonalPrism::to_json(json &j) const { j = {{"radius", 0.5 * box.x()}, {"length", box.z()}}; }

double HexagonalPrism::innerRadius() const { return 0.5 * box.x(); }
double HexagonalPrism::outerRadius() const { return 0.5 * box.y(); }
double HexagonalPrism::height() const { return box.z(); }

// =============== Cylinder ===============

Cylinder::Cylinder(double radius, double height) : radius(radius), height(height) {
    boundary_conditions =
        BoundaryCondition(Coordinates::ORTHOGONAL, {Boundary::FIXED, Boundary::FIXED, Boundary::PERIODIC});
}

Point Cylinder::getLength() const { return {2.0 * radius, 2.0 * radius, height}; }

double Cylinder::getVolume(int) const { return pc::pi * radius * radius * height; }

Point Cylinder::setVolume(double volume, const VolumeMethod method) {
    const double old_volume = getVolume();
    double alpha;
    Point box_scaling;

    switch (method) {
    case VolumeMethod::ISOTROPIC:
        alpha = std::cbrt(volume / old_volume);
        radius *= alpha;
        height *= alpha;
        box_scaling = {alpha, alpha, alpha};
        break;
    case VolumeMethod::XY: // earlier wrongly named ISOTROPIC!
        alpha = std::sqrt(volume / old_volume);
        radius *= alpha;
        box_scaling = {alpha, alpha, 1.0};
        break;
    case VolumeMethod::Z:
        alpha = volume / old_volume;
        height *= alpha;
        box_scaling = {1.0, 1.0, alpha};
        break;
    case VolumeMethod::ISOCHORIC:
        // height is scaled by 1/alpha/alpha, radius is scaled by alpha
        alpha = std::cbrt(volume / old_volume);
        radius *= alpha;
        height /= (alpha * alpha);
        box_scaling = {alpha, alpha, 1.0 / (alpha * alpha)};
        volume = old_volume;
        break;
    default:
        throw std::invalid_argument("unsupported volume scaling method for the cylindrical geometry");
    }
    assert(std::fabs(getVolume() - volume) < 1e-6 && "error setting sphere volume");
    return box_scaling;
}

void Cylinder::boundary(Point &a) const {
    // z-pbc
    if (std::fabs(a.z()) > 0.5 * height)
        a.z() -= height * anint(a.z() / height);
}

bool Cylinder::collision(const Point &a) const {
    return std::fabs(a.z()) > 0.5 * height || a.x() * a.x() + a.y() * a.y() > radius * radius;
}

Point Cylinder::vdist(const Point &a, const Point &b) const {
    Point distance(a - b);
    if (distance.z() > 0.5 * height)
        distance.z() -= height;
    else if (distance.z() < -0.5 * height)
        distance.z() += height;
    return distance;
}

void Cylinder::randompos(Point &m, Random &rand) const {
    double r2 = radius * radius, d = 2.0 * radius;
    m.z() = (rand() - 0.5) * height;
    do {
        // x,y shall always generated as a pair; 78.5% chance to hit
        m.x() = (rand() - 0.5) * d;
        m.y() = (rand() - 0.5) * d;
    } while (m.x() * m.x() + m.y() * m.y() > r2);
}

void Cylinder::from_json(const json &j) {
    radius = j.at("radius").get<double>();
    height = j.at("length").get<double>();
}

void Cylinder::to_json(json &j) const { j = {{"radius", radius}, {"length", height}}; }

// =============== Truncated Octahedron ===============

TruncatedOctahedron::TruncatedOctahedron(double side) : side(side) {
    // the current implementation is hardcoded as bellow and ignores other periodic_directions settings
    boundary_conditions =
        BoundaryCondition(Coordinates::TRUNC_OCTAHEDRAL, {Boundary::PERIODIC, Boundary::PERIODIC, Boundary::PERIODIC});
}

Point TruncatedOctahedron::getLength() const {
    // todo check orientation in xyz
    return Point::Constant(2.0 * std::sqrt(2.0) * side); // distance between opposite square faces
}

double TruncatedOctahedron::getVolume(int) const { return std::sqrt(128.0) * side * side * side; }

Point TruncatedOctahedron::setVolume(double volume, const VolumeMethod method) {
    const double old_side = side;
    Point box_scaling;

    if (method == VolumeMethod::ISOTROPIC) {
        side = std::cbrt(volume / std::sqrt(128.0));
        assert(std::fabs(getVolume() - volume) < 1e-6 && "error setting sphere volume");
    } else {
        throw std::invalid_argument("unsupported volume scaling method for the truncated-octahedral geometry");
    }
    box_scaling.setConstant(side / old_side);
    assert(fabs(getVolume() - volume) < 1e-6);
    return box_scaling;
}

Point TruncatedOctahedron::vdist(const Point &a, const Point &b) const {
    Point distance(a - b);
    boundary(distance);
    return distance;
}

bool TruncatedOctahedron::collision(const Point &a) const {
    const double sqrtThreeI = 1.0 / std::sqrt(3.0);
    const double origin_to_square_face = std::sqrt(2.0) * side;
    const double origin_to_hexagonal_face = std::sqrt(1.5) * side;

    // ugly
    if (std::fabs(a.dot(Point(1.0, 0.0, 0.0))) > origin_to_square_face)
        return true;
    if (std::fabs(a.dot(Point(0.0, 1.0, 0.0))) > origin_to_square_face)
        return true;
    if (std::fabs(a.dot(Point(0.0, 0.0, 1.0))) > origin_to_square_face)
        return true;
    if (std::fabs(a.dot(Point(1.0, 1.0, 1.0) * sqrtThreeI)) > origin_to_hexagonal_face)
        return true;
    if (std::fabs(a.dot(Point(1.0, 1.0, -1.0) * sqrtThreeI)) > origin_to_hexagonal_face)
        return true;
    if (std::fabs(a.dot(Point(1.0, -1.0, -1.0) * sqrtThreeI)) > origin_to_hexagonal_face)
        return true;
    if (std::fabs(a.dot(Point(1.0, -1.0, 1.0) * sqrtThreeI)) > origin_to_hexagonal_face)
        return true;
    return false;
}

void TruncatedOctahedron::boundary(Point &a) const {
    const double sqrtThreeI = 1.0 / std::sqrt(3.0);
    const double square_face_distance = std::sqrt(8.0) * side;
    const double hexagonal_face_distance = std::sqrt(6.0) * side;
    const Point unitvXYZ = Point(1, 1, 1) * sqrtThreeI;
    const Point unitvXiYZ = Point(1, 1, -1) * sqrtThreeI;
    const Point unitvXYiZ = Point(1, -1, -1) * sqrtThreeI;
    const Point unitvXYZi = Point(1, -1, 1) * sqrtThreeI;

    // todo improve
    bool outside = false;
    do {
        outside = false;
        double tmp = a.dot(unitvXYZ);
        if (std::fabs(tmp) > hexagonal_face_distance * 0.5) {
            a -= hexagonal_face_distance * anint(tmp / hexagonal_face_distance) * unitvXYZ;
            outside = true;
        }
        tmp = a.dot(unitvXiYZ);
        if (std::fabs(tmp) > hexagonal_face_distance * 0.5) {
            a -= hexagonal_face_distance * anint(tmp / hexagonal_face_distance) * unitvXiYZ;
            outside = true;
        }
        tmp = a.dot(unitvXYiZ);
        if (std::fabs(tmp) > hexagonal_face_distance * 0.5) {
            a -= hexagonal_face_distance * anint(tmp / hexagonal_face_distance) * unitvXYiZ;
            outside = true;
        }
        tmp = a.dot(unitvXYZi);
        if (std::fabs(tmp) > hexagonal_face_distance * 0.5) {
            a -= hexagonal_face_distance * anint(tmp / hexagonal_face_distance) * unitvXYZi;
            outside = true;
        }
    } while (outside);

    if (std::fabs(a.x()) > square_face_distance * 0.5)
        a.x() -= square_face_distance * anint(a.x() / square_face_distance);

    if (std::fabs(a.y()) > square_face_distance * 0.5)
        a.y() -= square_face_distance * anint(a.y() / square_face_distance);

    if (std::fabs(a.z()) > square_face_distance * 0.5)
        a.z() -= square_face_distance * anint(a.z() / square_face_distance);
}

void TruncatedOctahedron::randompos(Point &pos, Random &rand) const {
    const double d = std::sqrt(10.0) * side; // use circumdiameter
    const double r2 = d * d / 4.0;
    do {
        do {
            pos.x() = (rand() - 0.5) * d;
            pos.y() = (rand() - 0.5) * d;
            pos.z() = (rand() - 0.5) * d;
        } while (pos.squaredNorm() > r2);
    } while (collision(pos));
}

void TruncatedOctahedron::from_json(const json &j) { side = j.at("radius").get<double>(); }

void TruncatedOctahedron::to_json(json &j) const { j = {{"radius", side}}; }

// =============== Chameleon==============

const std::map<std::string, Variant> Chameleon::names = {{{"cuboid", Variant::CUBOID},
                                                          {"cylinder", Variant::CYLINDER},
                                                          {"slit", Variant::SLIT},
                                                          {"sphere", Variant::SPHERE},
                                                          {"hexagonal", Variant::HEXAGONAL},
                                                          {"octahedron", Variant::OCTAHEDRON},
                                                          {"hypersphere2d", Variant::HYPERSPHERE2D}}};

void from_json(const json &j, Chameleon &g) {
    try {
        g.from_json(j);
    }
    catch (std::exception &e) {
        usageTip.pick("geometry");
        throw ConfigurationError("geometry construction error: {}", e.what());
    }
}

void to_json(json &j, const Chameleon &g) { g.to_json(j); }

ParticleVector mapParticlesOnSphere(const ParticleVector &source) {
    assert(source.size() > 1);
    Average<double> radius; // average radial distance
    ParticleVector destination = source;
    Point COM = massCenter(source.begin() + 1, source.end());
    for (size_t i = 1; i < source.size(); i++) {
        destination[i].pos = source[i].pos - COM; // make COM origin
        double r = destination[i].pos.norm();     // radial distance
        destination[i].pos /= r;                  // normalize to unit vector
        radius += r;                         // radius is the average r
    }
    destination[0].pos.setZero();
    for (auto &i : destination) // scale positions to surface of sphere
        i.pos = i.pos * radius.avg() + COM;

    // rmsd, skipping the first particle
    double _rmsd =
        rootMeanSquareDeviation(destination.begin() + 1, destination.end(), source.begin() + 1,
                                [](const Particle &a, const Particle &b) { return (a.pos - b.pos).squaredNorm(); });

    faunus_logger->info("{} particles mapped on sphere of radius {} with RMSD {} {}; the first particle ({}) is a "
                        "dummy and COM placeholder",
                        destination.size(), radius.avg(), _rmsd, u8::angstrom,
                        Faunus::atoms.at(destination.at(0).id).name);
    return destination;
}

ShapeDescriptors::ShapeDescriptors(const Tensor &gyration_tensor) {
    const auto principal_moment = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(gyration_tensor).eigenvalues();
    gyration_radius_squared = gyration_tensor.trace();
    asphericity = 3.0 / 2.0 * principal_moment.z() - gyration_radius_squared / 2.0;
    acylindricity = principal_moment.y() - principal_moment.x();
    relative_shape_anisotropy = (asphericity * asphericity + 3.0 / 4.0 * acylindricity * acylindricity) /
                                (gyration_radius_squared * gyration_radius_squared);
}

ShapeDescriptors &ShapeDescriptors::operator+=(const ShapeDescriptors &other) {
    gyration_radius_squared += other.gyration_radius_squared;
    asphericity += other.asphericity;
    acylindricity += other.acylindricity;
    relative_shape_anisotropy += other.relative_shape_anisotropy;
    return *this;
}

ShapeDescriptors ShapeDescriptors::operator*(const double scale) const {
    auto scaled = *this;
    scaled.gyration_radius_squared *= scale;
    scaled.asphericity *= scale;
    scaled.acylindricity *= scale;
    scaled.relative_shape_anisotropy *= scale;
    return scaled;
}

void to_json(json &j, const ShapeDescriptors &shape) {
    j = {{"Rg = √⟨s²⟩", std::sqrt(shape.gyration_radius_squared)},
         {"asphericity (b)", shape.asphericity},
         {"acylindricity (c)", shape.acylindricity},
         {"relative shape anisotropy (𝜅²)", shape.relative_shape_anisotropy}};
}

TEST_CASE("[Faunus] ShapeDescriptors") {
    using doctest::Approx;
    std::vector<Point> positions = {{0, 0, 0}, {1, 0, 0}};
    std::vector<double> masses = {1, 1};
    Point origin = {0, 0, 0};
    auto gyration_tensor = gyration(positions.begin(), positions.end(), masses.begin(), origin);
    CHECK(gyration_tensor(0, 0) == Approx(0.5));

    auto shape = ShapeDescriptors(gyration_tensor);
    CHECK(shape.relative_shape_anisotropy == Approx(1.0));

    positions = {{0, 1, 0}, {1, 0, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 1, 1}, {1, 0, 1}, {-1, 0, 1}, {0, -1, 1}};
    masses = {1, 1, 1, 1, 1, 1, 1, 1};
    gyration_tensor = gyration(positions.begin(), positions.end(), masses.begin(), origin);
    shape = ShapeDescriptors(gyration_tensor);
    CHECK(shape.relative_shape_anisotropy == Approx(0.0));
}

TEST_CASE("[Faunus] HexagonalPrismToCuboid") {
    using doctest::Approx;
    double radius = 2.0, height = 20.0;
    double side = 2.0 / std::sqrt(3.0) * radius;
    HexagonalPrism hexagonal_prism(side, height);
    ParticleVector p(6); // corners of a hexagon
    p[0].pos = {0, 1, 0};
    p[1].pos = {0.866, 0.5, 0};
    p[2].pos = {0.866, -0.5, 0};
    p[3].pos = {0, -1, 0};
    p[4].pos = {-0.866, -0.5, 0};
    p[5].pos = {-0.866, 0.5, 0};
    auto [cuboid, p_new] = HexagonalPrismToCuboid(hexagonal_prism, p);
    CHECK(p_new.size() == 12);
    CHECK(cuboid.getLength().x() == Approx(radius * 2.0));
    CHECK(cuboid.getLength().y() == Approx(side * 3.0));
    CHECK(cuboid.getLength().z() == Approx(height));
    CHECK(cuboid.getVolume() == Approx(2.0 * hexagonal_prism.getVolume()));

    std::vector<Point> positions = {{0, 1, 0},           {0.866, 0.5, 0},  {0.866, -0.5, 0},   {0, -1, 0},
                                    {-0.866, -0.5, 0},   {-0.866, 0.5, 0}, {2, -2.4641, 0},    {-1.134, -2.9641, 0},
                                    {-1.134, 2.9641, 0}, {2, 2.4641, 0},   {1.134, 2.9641, 0}, {1.134, -2.9641, 0}};
    size_t i = 0;
    for (auto &particle : p_new) { // compared actual positions w. expected positions
        CHECK(particle.pos.x() == Approx(positions[i].x()));
        CHECK(particle.pos.y() == Approx(positions[i].y()));
        CHECK(particle.pos.z() == Approx(positions[i].z()));
        i++;
    }
}

Chameleon::Chameleon(const Variant type) {
    makeGeometry(type);
    _setLength(geometry->getLength());
}

Chameleon::Chameleon(const GeometryImplementation &geo, const Variant type) : geometry(geo.clone()), _type(type) {
    _setLength(geometry->getLength());
}

Point Chameleon::getLength() const {
    assert(geometry);
    return geometry->getLength();
}

void Chameleon::setLength(const Point &l) {
    assert(geometry);
    _setLength(l);
    // ugly
    if (type == Variant::CUBOID) {
        Cuboid &cuboid = dynamic_cast<Cuboid &>(*geometry);
        cuboid.setLength(l);
    } else {
        throw std::runtime_error("setLength allowed only for the Cuboid geometry");
    }
}

void Chameleon::makeGeometry(const Variant type) {
    switch (type) {
    case Variant::CUBOID:
        geometry = std::make_unique<Cuboid>();
        break;
    case Variant::SLIT:
        geometry = std::make_unique<Slit>();
        break;
    case Variant::SPHERE:
        geometry = std::make_unique<Sphere>();
        break;
    case Variant::CYLINDER:
        geometry = std::make_unique<Cylinder>();
        break;
    case Variant::HEXAGONAL:
        geometry = std::make_unique<HexagonalPrism>();
        break;
    case Variant::OCTAHEDRON:
        geometry = std::make_unique<TruncatedOctahedron>();
        break;
    case Variant::HYPERSPHERE2D:
        geometry = std::make_unique<Hypersphere2d>();
        break;
    default:
        throw std::invalid_argument("unknown geometry");
    }
    _type = type;
}

void Chameleon::from_json(const json &j) {
    std::tie(_name, _type) = variantName(j);
    makeGeometry(_type);
    geometry->from_json(j);
    _setLength(geometry->getLength());
}

void Chameleon::to_json(json &j) const {
    assert(geometry);
    geometry->to_json(j);
    j["type"] = name;
}

void Chameleon::_setLength(const Point &l) {
    len = l;
    len_half = l * 0.5;
    len_inv = l.cwiseInverse();
    len_or_zero.setZero();

    // this is to speed up the `sqdist()` by avoiding branching when testing
    // for PBC in each direction. The variable `len_or_zero` either equals
    // `len` for PBC or zero if not.
    if (geometry->boundary_conditions.coordinates == Coordinates::ORTHOGONAL)
        for (size_t i = 0; i < 3; i++)
            len_or_zero[i] = len[i] * (geometry->boundary_conditions.direction[i] == Boundary::PERIODIC);
}

double Chameleon::getVolume(int dim) const {
    assert(geometry);
    return geometry->getVolume(dim);
}

Point Chameleon::setVolume(double V, VolumeMethod method) {
    auto scale = geometry->setVolume(V, method);
    _setLength(geometry->getLength());
    return scale;
}

Chameleon::VariantName Chameleon::variantName(const std::string &name) {
    if (auto it = names.find(name); it == names.end()) {
        throw std::runtime_error("unknown geometry: " + name);
    } else {
        return *it;
    }
}

Chameleon::VariantName Chameleon::variantName(const json &j) { return variantName(j.at("type").get<std::string>()); }

Chameleon &Chameleon::operator=(const Chameleon &geo) {
    if (&geo != this) {
        GeometryBase::operator=(geo);
        len = geo.len;
        len_half = geo.len_half;
        len_inv = geo.len_inv;
        len_or_zero = geo.len_or_zero;
        _type = geo._type;
        _name = geo._name;
        geometry = geo.geometry != nullptr ? geo.geometry->clone() : nullptr;
    }
    return *this;
}

const BoundaryCondition &Chameleon::boundaryConditions() const { return geometry->boundary_conditions; }

Chameleon::Chameleon(const Chameleon &geo)
    : GeometryBase(geo), len(geo.len), len_half(geo.len_half), len_inv(geo.len_inv),
      geometry(geo.geometry != nullptr ? geo.geometry->clone() : nullptr), _type(geo._type), _name(geo._name) {}

std::shared_ptr<GeometryImplementation> Chameleon::asSimpleGeometry() const { return geometry->clone(); }

TEST_CASE("[Faunus] spherical coordinates") {
    using doctest::Approx;

    Point sph1 = {2, 0.5, -0.3};
    auto pnt1 = rtp2xyz(sph1); // sph --> cart
    auto sph2 = xyz2rtp(pnt1); // cart --> sph

    CHECK(pnt1.norm() == Approx(2));
    CHECK(sph1.x() == Approx(sph2.x()));
    // CHECK( sph1.y() == Approx(sph2.y()));
    // CHECK( sph1.z() == Approx(sph2.z()));
}

TEST_CASE("[Faunus] Geometry") {
    using doctest::Approx;
    Random slump;

    SUBCASE("cuboid") {
        double x = 2, y = 3, z = 4;
        Cuboid geo({x, y, z});
        CHECK(geo.getVolume() == doctest::Approx(x * y * z));

        // check boundaries and pbc
        Point a(1.1, 1.5, -2.001);
        CHECK(geo.collision(a) == true);
        geo.getBoundaryFunc()(a);
        CHECK(geo.collision(a) == false);
        CHECK(a.x() == Approx(-0.9));  // x has been wrapped
        CHECK(a.y() == Approx(1.5));   // y is unchanged
        CHECK(a.z() == Approx(1.999)); // z has been wrapped
        a.y() = 1.51;                  // move y out of box
        geo.getBoundaryFunc()(a);      // wrap around boundary
        CHECK(a.y() == Approx(-1.49)); // check y-boundary
        a.y() = 1.5;                   // restore

        // check distances
        Point distance = geo.vdist({0.1, 0.5, -1.001}, a);
        CHECK(distance.x() == Approx(1.0));
        CHECK(distance.y() == Approx(-1.0));
        CHECK(distance.z() == Approx(1.0));
        CHECK(geo.vdist({1, 2, 3}, a) == geo.getDistanceFunc()({1, 2, 3}, a));

        // check that geometry is properly inscribed in a cuboid
        Point box = geo.getLength();
        CHECK(box.x() == Approx(x));
        CHECK(box.y() == Approx(y));
        CHECK(box.z() == Approx(z));

        // check random position
        Point c(x + 1, y + 1, z + 1); // out of the box
        bool container_overlap = false;
        for (int i = 0; i < 1e4; i++) {
            geo.randompos(c, slump);
            if (geo.collision(c))
                container_overlap = true;
        }
        CHECK(container_overlap == false);

        // volume scaling
        double sf = 2.;
        auto scaling = geo.setVolume(sf * sf * sf * x * y * z);
        CHECK(geo.getVolume() == doctest::Approx(sf * sf * sf * x * y * z));
        CHECK(geo.getLength().x() == Approx(sf * x));
        CHECK(geo.getLength().y() == Approx(sf * y));
        CHECK(geo.getLength().z() == Approx(sf * z));
        CHECK(scaling.x() == Approx(sf));
        CHECK(scaling.y() == Approx(sf));
        CHECK(scaling.z() == Approx(sf));

        // check json
        geo.from_json(R"( {"type": "cuboid", "length": [2.5,3.5,4.5]} )"_json);
        CHECK(geo.getVolume() == doctest::Approx(2.5 * 3.5 * 4.5));
    }

    SUBCASE("slit") {
        double x = 2, y = 4, z = 3;
        Slit geo(x, y, z);
        CHECK(geo.getVolume() == doctest::Approx(x * y * z));

        // check boundaries and pbc
        Point a(1.1, -2.001, 1.499);
        CHECK(geo.collision(a) == true);
        geo.getBoundaryFunc()(a);
        CHECK(geo.collision(a) == false);
        CHECK(a.x() == Approx(-0.9));
        CHECK(a.y() == Approx(1.999));
        CHECK(a.z() == Approx(1.499));
        Point b = a;
        geo.boundary(b);
        CHECK(a == b);

        Point outsize_z = {1.1, -2.001, 1.5001};
        CHECK(geo.collision(outsize_z) == true);
        geo.boundary(outsize_z);
        CHECK(geo.collision(outsize_z) == true);

        Point c(0, 0, -0.51 * z);

        // check distances
        Point distance = geo.vdist({0.1, -1.001, -0.501}, a);
        CHECK(distance.x() == Approx(1.0));
        CHECK(distance.y() == Approx(1.0));
        CHECK(distance.z() == Approx(-2.0));
        CHECK(geo.vdist({1, 2, 3}, a) == geo.getDistanceFunc()({1, 2, 3}, a));

        // check that geometry is properly enscribed in a cuboid
        Point box = geo.getLength();
        CHECK(box.x() == Approx(x));
        CHECK(box.y() == Approx(y));
        CHECK(box.z() == Approx(z));

        // check random position
        Point d(x + 1, y + 1, z + 1); // out of the box
        bool containerOverlap = false;
        for (int i = 0; i < 1e4; i++) {
            geo.randompos(d, slump);
            if (geo.collision(d))
                containerOverlap = true;
        }
        CHECK(containerOverlap == false);

        // volume scaling
        double sf = 2.;
        auto scaling = geo.setVolume(sf * sf * x * y * z, VolumeMethod::XY);
        CHECK(geo.getVolume() == doctest::Approx(sf * sf * x * y * z));
        CHECK(geo.getLength().x() == Approx(sf * x));
        CHECK(geo.getLength().y() == Approx(sf * y));
        CHECK(geo.getLength().z() == Approx(z));
        CHECK(scaling.x() == Approx(sf));
        CHECK(scaling.y() == Approx(sf));
        CHECK(scaling.z() == Approx(1.0));

        // check json
        geo.from_json(R"( {"type": "cuboid", "length": [2.5,3.5,4.5]} )"_json);
        CHECK(geo.getVolume() == doctest::Approx(2.5 * 3.5 * 4.5));
    }

    SUBCASE("sphere") {
        double radius = 5.;
        Sphere geo(radius);
        CHECK(geo.getVolume() == doctest::Approx(4. / 3. * pc::pi * radius * radius * radius));

        // check boundaries
        CHECK(geo.collision({5.01, 0, 0}) == true);
        CHECK(geo.collision({4.99, 0, 0}) == false);
        Point a(radius - 1, 0, -0.5 * radius);
        Point b = a;
        geo.boundary(a);
        CHECK(a == b);

        // check distances
        Point distance = geo.vdist({3.0, 1.0, -2.0}, {-3.0, -1.0, 2.0});
        CHECK(distance.x() == Approx(6.0));
        CHECK(distance.y() == Approx(2.0));
        CHECK(distance.z() == Approx(-4.0));

        // check that geometry is properly enscribed in a cuboid
        Point box = geo.getLength();
        CHECK(box.x() == Approx(10));
        CHECK(box.y() == Approx(10));
        CHECK(box.z() == Approx(10));

        // check random position
        Point c(radius + 1, radius + 1, radius + 1); // out of the box
        bool container_overlap = false;
        for (int i = 0; i < 1e4; i++) {
            geo.randompos(c, slump);
            if (geo.collision(c))
                container_overlap = true;
        }
        CHECK(container_overlap == false);

        // volume scaling
        geo.setVolume(123.4);
        CHECK(geo.getVolume() == Approx(123.4));
        CHECK_THROWS_AS(geo.setVolume(100., VolumeMethod::ISOCHORIC), std::invalid_argument);
        CHECK_THROWS_AS(geo.setVolume(100., VolumeMethod::XY), std::invalid_argument);

        // check json
        geo.from_json(R"( { "type": "sphere", "radius": 2.0 } )"_json);
        CHECK(geo.getVolume() == doctest::Approx(4. / 3. * pc::pi * 2.0 * 2.0 * 2.0));
    }

    SUBCASE("cylinder") {
        double radius = 1., volume = 1.;
        double height = volume / (pc::pi * radius * radius);
        Point box;
        Cylinder geo(radius, height);

        // check boundaries
        CHECK(geo.getVolume() == Approx(volume));
        CHECK(geo.collision({-1.01 * radius, 0, 0}) == true);
        CHECK(geo.collision({0.99 * radius, 0, 0}) == false);
        CHECK(geo.collision({-0.99 * radius, 0.15 * radius, 0}) == true);
        CHECK(geo.collision({0, 0, -0.51 * height}) == true);
        CHECK(geo.collision({0, 0, 0.49 * height}) == false);

        // check that geometry is properly enscribed in a cuboid
        box = geo.getLength();
        CHECK(box.x() == Approx(2 * radius));
        CHECK(box.y() == Approx(2 * radius));
        CHECK(box.z() == Approx(height));

        // check random position
        Point a(2. * radius, 0, 0); // out of the box
        bool container_overlap = false;
        for (int i = 0; i < 1e4; i++) {
            geo.randompos(a, slump);
            if (geo.collision(a))
                container_overlap = true;
        }
        CHECK(container_overlap == false);

        // volume scaling
        geo.setVolume(9.0, VolumeMethod::XY);
        CHECK(geo.getVolume() == Approx(9.0));
        box = geo.getLength();
        CHECK(box.x() == Approx(3 * 2 * radius));
        CHECK(box.y() == Approx(3 * 2 * radius));
        CHECK(box.z() == Approx(height));

        // check json
        json j = {{"type", "cylinder"}, {"radius", 2.0}, {"length", 2 / pc::pi}};
        geo.from_json(j);
        CHECK(geo.getVolume() == doctest::Approx(8.0));
    }

    SUBCASE("hexagonal prism") {
        double side = 1., volume = 1.;
        double outer_radius = side, inner_radius = side * std::sqrt(3.0) / 2.;
        double height = volume / (3. * outer_radius * inner_radius);
        Point box;
        HexagonalPrism geo(side, height);

        CHECK(geo.getVolume() == Approx(volume));
        CHECK(geo.collision({-1.01 * inner_radius, 0, 0}) == true);
        CHECK(geo.collision({0.99 * inner_radius, 0, 0}) == false);
        CHECK(geo.collision({0.0, -1.01 * outer_radius, 0}) == true);
        CHECK(geo.collision({0.0, 0.99 * outer_radius, 0}) == false);
        CHECK(geo.collision({0.99 * std::cos(pc::pi / 3.) * inner_radius, 0.99 * std::sin(pc::pi / 3.) * inner_radius,
                             0}) == false);
        CHECK(geo.collision({1.01 * std::cos(pc::pi / 3.) * inner_radius, 1.01 * std::sin(pc::pi / 3.) * inner_radius,
                             0}) == true);
        CHECK(geo.collision({0, 0, -0.51 * height}) == true);
        CHECK(geo.collision({0, 0, 0.49 * height}) == false);

        // check that geometry is properly inscribed in a cuboid
        box = geo.getLength();
        CHECK(box.x() == Approx(2 * inner_radius));
        CHECK(box.y() == Approx(2 * outer_radius));
        CHECK(box.z() == Approx(height));

        // check random position
        Point a;
        bool container_overlap = false;
        for (int i = 0; i < 1e4; i++) {
            geo.randompos(a, slump);
            if (geo.collision(a))
                container_overlap = true;
        }
        CHECK(container_overlap == false);

        json j = {{"type", "hexagonal"}, {"radius", 3 * inner_radius}, {"length", 5 * height}};
        geo.from_json(j);
        CHECK(geo.getVolume() == Approx(9 * 5 * volume));
    }
}

TEST_CASE("[Faunus] Chameleon") {

    using doctest::Approx;
    Random slump;

    //! function compares if Chamelon's and Geometry's boundary methods produce the same result
    //! using n random points
    auto compare_boundary = [&slump](Chameleon &chameleon, GeometryImplementation &geo, Cuboid &box, int n = 100) {
        Point a, b;
        for (int i = 0; i < n; i++) {
            box.randompos(a, slump);
            b = a;
            chameleon.boundary(a);
            geo.boundary(b);
            CHECK(a.x() == Approx(b.x()));
            CHECK(a.y() == Approx(b.y()));
            CHECK(a.z() == Approx(b.z()));
        }
    };

    //! function compares if Chamelon's and Geometry's vdist methods produce the same result
    //! using n random points
    auto compare_vdist = [&slump](Chameleon &chameleon, GeometryImplementation &geo, Cuboid &box, int n = 100) {
        Point a, b, d_cham, d_geo;
        for (int i = 0; i < n; i++) {
            box.randompos(a, slump);
            box.randompos(b, slump);
            d_cham = chameleon.vdist(a, b);
            d_geo = geo.vdist(a, b);
            CHECK(d_cham.x() == Approx(d_geo.x()));
            CHECK(d_cham.y() == Approx(d_geo.y()));
            CHECK(d_cham.z() == Approx(d_geo.z()));
            CHECK(chameleon.sqdist(a, b) == Approx(d_cham.squaredNorm()));
        }
    };

    SUBCASE("cuboid") {
        double x = 2.0, y = 3.0, z = 4.0;
        Point box_size = std::cbrt(2.0) * Point(x, y, z);
        Cuboid box(box_size);
        Cuboid geo({x, y, z});
        Chameleon chameleon(geo, Variant::CUBOID);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
        CHECK(geo.boundary_conditions.isPeriodic()[0] == true);
        CHECK(geo.boundary_conditions.isPeriodic()[1] == true);
        CHECK(geo.boundary_conditions.isPeriodic()[2] == true);
    }

    SUBCASE("slit") {
        double x = 2.0, y = 3.0, z = 4.0;
        Point box_size = std::cbrt(2.0) * Point(x, y, z);
        Cuboid box(box_size);
        Slit geo(x, y, z);
        Chameleon chameleon(geo, Variant::SLIT);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
        CHECK(geo.boundary_conditions.isPeriodic()[0] == true);
        CHECK(geo.boundary_conditions.isPeriodic()[1] == true);
        CHECK(geo.boundary_conditions.isPeriodic()[2] == false);
    }

    SUBCASE("sphere") {
        double radius = 10.0;
        Point box_size;
        box_size.setConstant(std::cbrt(2.0) * 2 * radius);
        Cuboid box(box_size);
        Sphere geo(radius);
        Chameleon chameleon(geo, Variant::SPHERE);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
        CHECK(geo.boundary_conditions.isPeriodic()[0] == false);
        CHECK(geo.boundary_conditions.isPeriodic()[1] == false);
        CHECK(geo.boundary_conditions.isPeriodic()[2] == false);
    }

    SUBCASE("cylinder") {
        double radius = 2.0, height = 10.0;
        Point box_size = std::cbrt(2.0) * Point(2 * radius, 2 * radius, height);
        Cuboid box(box_size);
        Cylinder geo(radius, height);
        Chameleon chameleon(geo, Variant::CYLINDER);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
        CHECK(geo.boundary_conditions.isPeriodic()[0] == false);
        CHECK(geo.boundary_conditions.isPeriodic()[1] == false);
        CHECK(geo.boundary_conditions.isPeriodic()[2] == true);
    }

    SUBCASE("hexagonal prism") {
        double edge = 5.0, height = 20.0;
        Point box_size = std::cbrt(2.0) * Point(2 * edge, 2 * edge, height); // a bit larger in x-direction
        Cuboid box(box_size);
        HexagonalPrism geo(edge, height);
        Chameleon chameleon(geo, Variant::HEXAGONAL);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
    }

    SUBCASE("truncated octahedron") {
        double edge = 5.0;
        Point box_size;
        box_size.setConstant(std::cbrt(2.0) * edge * std::sqrt(5.0 / 2.0)); // enlarged circumradius
        Cuboid box(box_size);
        TruncatedOctahedron geo(edge);
        Chameleon chameleon(geo, Variant::OCTAHEDRON);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
    }

    SUBCASE("Cereal serialisation") {
        double x = 2.0, y = 3.0, z = 4.0;
        std::ostringstream os(std::stringstream::binary);
        { // write
            Cuboid geo({x, y, z});
            cereal::BinaryOutputArchive archive(os);
            archive(geo);
        }

        { // read
            Cuboid geo({10, 20, 30});
            std::istringstream in(os.str());
            cereal::BinaryInputArchive archive(in);
            archive(geo);
            CHECK(geo.getLength().x() == Approx(x));
            CHECK(geo.getLength().y() == Approx(y));
            CHECK(geo.getLength().z() == Approx(z));
        }
    }
}

TEST_CASE("[Faunus] anyCenter") {
    Chameleon cyl = json({{"type", "cuboid"}, {"length", 100}, {"radius", 20}});
    std::vector<Particle> p;

    CHECK(!atoms.empty()); // set in a previous test
    p.push_back(atoms[0]);
    p.push_back(atoms[0]);

    p.front().pos = {10, 10, -10};
    p.back().pos = {15, -10, 10};

    Point cm = Geometry::massCenter(p.begin(), p.end(), cyl.getBoundaryFunc());
    CHECK(cm.x() == doctest::Approx(12.5));
    CHECK(cm.y() == doctest::Approx(0));
    CHECK(cm.z() == doctest::Approx(0));
}

TEST_CASE("[Faunus] rootMeanSquareDeviation") {
    std::vector<double> v1 = {1.3, 4.4, -1.1};
    std::vector<double> v2 = {1.1, 4.6, -1.0};
    auto f = [](double a, double b) { return std::pow(a - b, 2); };
    double rmsd = Geometry::rootMeanSquareDeviation(v1.begin(), v1.end(), v2.begin(), f);
    CHECK(rmsd == doctest::Approx(0.17320508075688745));
}

} // namespace Faunus::Geometry
