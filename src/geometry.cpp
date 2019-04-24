#include "geometry.h"

namespace Faunus {
    Point ranunit(Random &rand, const Point &dir) {
        double r2;
        Point p;
        do {
            for (size_t i=0; i<3; i++)
                p[i] = ( rand()-0.5 ) * dir[i];
            r2 = p.squaredNorm();
        } while ( r2 > 0.25 );
        return p / std::sqrt(r2);
    }

    Point ranunit_polar(Random &rand) {
        return rtp2xyz( {1, 2*pc::pi*rand(), std::acos(2*rand()-1)} );
    }

    Point xyz2rth(const Point &p, const Point &origin, const Point &dir, const Point &dir2) {
        assert( fabs(dir.norm()-1.0)<1e-6 );
        assert( fabs(dir2.norm()-1.0)<1e-6 );
        assert( fabs(dir.dot(dir2))<1e-6 ); // check if unit-vectors are perpendicular
        Point xyz = p - origin;
        double h = xyz.dot(dir);
        Point xy = xyz - dir*h;
        double x = xy.dot(dir2);
        Point y = xy - dir2*x;
        double theta = std::atan2(y.norm(),x);
        double radius = xy.norm();
        return {radius,theta,h};
    }

    Point xyz2rtp(const Point &p, const Point &origin) {
        Point xyz = p - origin;
        double radius = xyz.norm();
        return {
            radius,
                std::atan2( xyz.y(), xyz.x() ),
                std::acos( xyz.z()/radius) };
    }

    Point rtp2xyz(const Point &rtp, const Point &origin) {
        return origin + rtp.x() * Point(
                std::cos(rtp.y()) * std::sin(rtp.z()),
                std::sin(rtp.y()) * std::sin(rtp.z()),
                std::cos(rtp.z()) );
    }

    namespace Geometry {

        // =============== GeometryBase ===============

        GeometryBase::~GeometryBase() {}


        // =============== Cuboid ===============

        Cuboid::Cuboid(const Point &p) {
            boundary_conditions = BoundaryCondition(ORTHOGONAL, {PERIODIC, PERIODIC, PERIODIC});
            setLength(p);
        }

        Cuboid::Cuboid(double x, double y, double z) : Cuboid(Point(x, y, z)) {
        }

        Cuboid::Cuboid(double x) : Cuboid(x, x, x) {
        }

        inline Point Cuboid::getLength() const {
            return box;
        }

        void Cuboid::setLength(const Point &len) {
            box = len;
            box_half = 0.5 * box;
            box_inv = box.cwiseInverse();
        }

        inline double Cuboid::getVolume(int) const {
            return box.x() * box.y() * box.z();
        }

        Point Cuboid::setVolume(double volume, const VolumeMethod method) {
            const double old_volume = getVolume();
            double alpha;
            Point new_box, box_scaling;
            switch (method) {
                case ISOTROPIC:
                    alpha = std::cbrt(volume / old_volume);
                    box_scaling = {alpha, alpha, alpha};
                    break;
                case XY:
                    alpha = std::sqrt(volume / old_volume);
                    box_scaling = {alpha, alpha, 1.0};
                    break;
                case ISOCHORIC:
                    // z is scaled by 1/alpha/alpha, x and y are scaled by alpha
                    alpha = std::cbrt(volume / old_volume);
                    box_scaling = {alpha, alpha, 1 / (alpha * alpha)};
                    break;
                default:
                    throw std::invalid_argument("unsupported volume scaling method for the cuboid geometry");
            }
            new_box = box.cwiseProduct(box_scaling);
            setLength(new_box);
            assert(fabs(getVolume() - volume) < 1e-6);
            return box_scaling; // this will scale any point to new volume
        }

        inline void Cuboid::boundary(Point &a) const {
            if (boundary_conditions.direction.x() == PERIODIC) {
                if (std::fabs(a.x()) > box_half.x())
                    a.x() -= box.x() * anint(a.x() * box_inv.x());
            }
            if (boundary_conditions.direction.y() == PERIODIC) {
                if (std::fabs(a.y()) > box_half.y())
                    a.y() -= box.y() * anint(a.y() * box_inv.y());
            }
            if (boundary_conditions.direction.z() == PERIODIC) {
                if (std::fabs(a.z()) > box_half.z())
                    a.z() -= box.z() * anint(a.z() * box_inv.z());
            }
        }

        inline Point Cuboid::vdist(const Point &a, const Point &b) const {
            Point distance(a - b);
            if (boundary_conditions.direction.x() == PERIODIC) {
                if (distance.x() > box_half.x())
                    distance.x() -= box.x();
                else if (distance.x() < -box_half.x())
                    distance.x() += box.x();
            }
            if (boundary_conditions.direction.y() == PERIODIC) {
                if (distance.y() > box_half.y())
                    distance.y() -= box.y();
                else if (distance.y() < -box_half.y())
                    distance.y() += box.y();
            }
            if (boundary_conditions.direction.z() == PERIODIC) {
                if (distance.z() > box_half.z())
                    distance.z() -= box.z();
                else if (distance.z() < -box_half.z())
                    distance.z() += box.z();
            }
            return distance;
        }

        inline void Cuboid::randompos(Point &m, Random &rand) const {
            m.x() = (rand() - 0.5) * box.x();
            m.y() = (rand() - 0.5) * box.y();
            m.z() = (rand() - 0.5) * box.z();
        }

        inline bool Cuboid::collision(const Point &a) const {
            bool collision =
                    std::fabs(a.x()) > box_half.x() ||
                    std::fabs(a.y()) > box_half.y() ||
                    std::fabs(a.z()) > box_half.z();
            return collision;
        }

        void Cuboid::from_json(const json &j) {
            box.setZero();

            auto m = j.at("length");
            if (m.is_number()) {
                double l = m.get<double>();
                setLength({l, l, l});
            } else if (m.is_array()) {
                if (m.size() == 3) {
                    setLength(m.get<Point>());
                } else {
                    // TODO warning
                }
            } else {
                // TODO warning
            }
        }

        void Cuboid::to_json(json &j) const {
            j = {{"length", box}};
        }


        // =============== Slit ===============

        Slit::Slit(const Point &p) : Tbase(p) {
            boundary_conditions = BoundaryCondition(ORTHOGONAL, {PERIODIC, PERIODIC, FIXED});
        }

        Slit::Slit(double x, double y, double z) : Slit(Point(x, y, z)) {
        }

        Slit::Slit(double x) : Slit(x, x, x) {
        }


        // =============== Sphere ===============

        Sphere::Sphere(double radius) : radius(radius) {
            boundary_conditions = BoundaryCondition(ORTHOGONAL, {FIXED, FIXED, FIXED});
        }

        Point Sphere::getLength() const {
            return {2 * radius, 2 * radius, 2 * radius};
        }

        inline double Sphere::getVolume(int dim) const {
            double result;
            switch (dim) {
                case 3:
                    result = 4.0 / 3.0 * pc::pi * radius * radius * radius; // volume
                    break;
                case 2:
                    result = 4.0 * pc::pi * radius * radius; // surface area
                    break;
                case 1:
                    result = 2 * radius; // diameter
                    break;
                default:
                    throw std::invalid_argument("unsupported volume dimension for the sphere: " +
                                             std::to_string(dim));
            }
            return result;
        }

        Point Sphere::setVolume(double volume, const VolumeMethod method) {
            const double old_radius = radius;
            Point box_scaling;
            if (method == ISOTROPIC) {
                radius = std::cbrt(volume / (4.0 / 3.0 * pc::pi));
                assert(std::fabs(getVolume() - volume) < 1e-6 && "error setting sphere volume");
            } else {
                throw std::invalid_argument("unsupported volume scaling method for the spherical geometry");
            }
            box_scaling.setConstant(radius / old_radius);
            return box_scaling;
        }

        inline void Sphere::boundary(Point &) const {
            // no pbc
        }

        inline bool Sphere::collision(const Point &a) const {
            bool collision = a.squaredNorm() > radius * radius;
            return collision;
        }

        inline Point Sphere::vdist(const Point &a, const Point &b) const {
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

        void Sphere::from_json(const json &j) {
            radius = j.at("radius").get<double>();
        }

        void Sphere::to_json(json &j) const {
            j = {{"radius", radius}};
        }


        // =============== Hypersphere 2D ===============

        Hypersphere2d::Hypersphere2d(double radius) : Sphere(radius) {
            boundary_conditions = BoundaryCondition(NON3D);
        }

        inline Point Hypersphere2d::vdist(const Point &a, const Point &b) const {
            // ugly but works, needs fixing though...
            Point distance3d(a - b);
            double angle = std::acos(a.dot(b) / radius / radius);
            double distance = radius * angle;
            return distance3d / distance3d.norm() * distance;
        }

        void Hypersphere2d::randompos(Point &m, Random &rand) const {
            Sphere::randompos(m, rand);
            m = m / m.norm() * radius;
        }

        inline bool Hypersphere2d::collision(const Point &a) const {
            bool collision = std::fabs(a.norm() - radius) > 1e-6;
            return collision;
        }


        // =============== Hexagonal Prism ===============

        const Eigen::Matrix3d HexagonalPrism::rhombic2cartesian = (Eigen::Matrix3d() <<
                std::cos(-pc::pi/6), 0.0, 0.0,
                std::sin(-pc::pi/6), 1.0, 0.0,
                0.0,                 0.0, 1.0).finished();

        const  Eigen::Matrix3d HexagonalPrism::cartesian2rhombic = rhombic2cartesian.inverse();

        HexagonalPrism::HexagonalPrism(double side, double height) {
            // the current implementation is hardcoded as bellow and ignores other periodic_directions settings
            boundary_conditions = BoundaryCondition(ORTHOHEXAGONAL, {PERIODIC, PERIODIC, FIXED});
            set_box(side, height);
        }

        inline Point HexagonalPrism::getLength() const {
            return box;
        }

        double HexagonalPrism::getVolume(int) const {
            return 3. / 4. * box.x() * box.y() * box.z(); // 3 * inner_radius * outer_radius * height
        }

        void HexagonalPrism::set_box(double side, double height) {
            box = {std::sqrt(3.) * side, 2 * side, height};
        }

        Point HexagonalPrism::setVolume(double volume, const VolumeMethod method) {
            const double old_volume = getVolume();
            double alpha;
            Point box_scaling;

            switch (method) {
                case ISOTROPIC:
                    alpha = std::cbrt(volume / old_volume);
                    box_scaling = {alpha, alpha, alpha};
                    break;
                case XY:
                    alpha = std::sqrt(volume / old_volume);
                    box_scaling = {alpha, alpha, 1.0};
                    break;
                case ISOCHORIC:
                    // radius is scaled by alpha, z is scaled by 1/alpha/alpha
                    alpha = std::cbrt(volume / old_volume);
                    box_scaling = {alpha, alpha, 1 / (alpha * alpha)};
                    break;
                default:
                    throw std::invalid_argument("unsupported volume scaling method for the hexagonal-prism geometry");
            }
            box = box.cwiseProduct(box_scaling);
            assert(fabs(getVolume() - volume) < 1e-6);
            return box_scaling;
        }

        inline Point HexagonalPrism::vdist(const Point &a, const Point &b) const {
            Point distance(a - b);
            boundary(distance);
            return distance;
        }

        inline bool HexagonalPrism::collision(const Point &a) const {
            const double height = box.z();
            const double outer_radius = 0.5 * box.y();

            // Hexagon can be divided into three rhombuses. Using natural rhombic coordinates, it can be
            // straightforwardly determined if the particular point is inside the rhombus. If the point is mirrored to
            // the first quadrant taking the absolute values of coordinates, only one rhombus needs to be evaluated.
            Point b = cartesian2rhombic * a.cwiseAbs();
            bool collision = b.z() > 0.5 * height || b.x() > outer_radius || b.y() > outer_radius;
            return collision;
        }

        inline void HexagonalPrism::boundary(Point &a) const {
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

        void HexagonalPrism::to_json(json &j) const {
            j = {{"radius", box.y()},
                 {"length", box.z()}};
        }


        // =============== Cylinder ===============

        Cylinder::Cylinder(double radius, double height) : radius(radius), height(height) {
            boundary_conditions = BoundaryCondition(ORTHOGONAL, {FIXED, FIXED, PERIODIC});
        }

        inline Point Cylinder::getLength() const {
            return {2 * radius, 2 * radius, height};
        }

        inline double Cylinder::getVolume(int) const {
            return pc::pi * radius * radius * height;
        }

        Point Cylinder::setVolume(double volume, const VolumeMethod method) {
            const double old_volume = getVolume();
            double alpha;
            Point box_scaling;

            switch (method) {
                case ISOTROPIC:
                    alpha = std::cbrt(volume / old_volume);
                    radius *= alpha;
                    height *= alpha;
                    box_scaling = {alpha, alpha, alpha};
                    break;
                case XY: // earlier wrongly named as ISOTROPIC!
                    alpha = std::sqrt(volume / old_volume);
                    radius *= alpha;;
                    box_scaling = {alpha, alpha, 1.0};
                    break;
                case ISOCHORIC:
                    // height is scaled by 1/alpha/alpha, radius is scaled by alpha
                    alpha = std::cbrt(volume / old_volume);
                    radius *= alpha;
                    height /= (alpha * alpha);
                    box_scaling = {alpha, alpha, 1.0 / (alpha * alpha)};
                    break;
                default:
                    throw std::invalid_argument("unsupported volume scaling method for the cylindrical geometry");
            }
            assert(std::fabs(getVolume() - volume) < 1e-6 && "error setting sphere volume");
            return box_scaling;
        }

        inline void Cylinder::boundary(Point &a) const {
            // z-pbc
            if (std::fabs(a.z()) > 0.5 * height)
                a.z() -= height * anint(a.z() / height);
        }

        bool Cylinder::collision(const Point &a) const {
            bool collision = std::fabs(a.z()) > 0.5 * height || a.x() * a.x() + a.y() * a.y() > radius * radius;
            return collision;
        }

        inline Point Cylinder::vdist(const Point &a, const Point &b) const {
            Point distance(a - b);
            if (distance.z() > 0.5 * height)
                distance.z() -= height;
            else if (distance.z() < -0.5 * height)
                distance.z() += height;
            return distance;
        }

        void Cylinder::randompos(Point &m, Random &rand) const {
            double r2 = radius * radius, d = 2 * radius;
            m.z() = (rand() - 0.5) * height;
            do {
                // x,y shall be always generated as a pair; 78.5% chance to hit
                m.x() = (rand() - 0.5) * d;
                m.y() = (rand() - 0.5) * d;
            } while (m.x() * m.x() + m.y() * m.y() > r2);
        }

        void Cylinder::from_json(const json &j) {
            radius = j.at("radius").get<double>();
            height = j.at("length").get<double>();
        }

        void Cylinder::to_json(json &j) const {
            j = {{"radius", radius},
                 {"length", height}};
        }


        // =============== Truncated Octahedron ===============

        TruncatedOctahedron::TruncatedOctahedron(double side) : side(side) {
            // the current implementation is hardcoded as bellow and ignores other periodic_directions settings
            boundary_conditions = BoundaryCondition(TRUNC_OCTAHEDRAL, {PERIODIC, PERIODIC, PERIODIC});
        }

        inline Point TruncatedOctahedron::getLength() const {
            // todo check orientation in xyz
            return Point::Constant(2. * std::sqrt(2.) * side); // distance between opposite square faces
        }

        double TruncatedOctahedron::getVolume(int) const {
            return std::sqrt(128.) * side * side * side;
        }

        Point TruncatedOctahedron::setVolume(double volume, const VolumeMethod method) {
            const double old_side = side;
            Point box_scaling;

            if (method == ISOTROPIC) {
                side = std::cbrt(volume / std::sqrt(128.));
                assert(std::fabs(getVolume() - volume) < 1e-6 && "error setting sphere volume");
            } else {
                throw std::invalid_argument("unsupported volume scaling method for the truncated-octahedral geometry");
            }
            box_scaling.setConstant(side / old_side);
            assert(fabs(getVolume() - volume) < 1e-6);
            return box_scaling;
        }

        inline Point TruncatedOctahedron::vdist(const Point &a, const Point &b) const {
            Point distance(a - b);
            boundary(distance);
            return distance;
        }

        inline bool TruncatedOctahedron::collision(const Point &a) const {
            const double sqrtThreeI = 1.0 / std::sqrt(3.0);
            const double origin_to_square_face = std::sqrt(2.) * side;
            const double origin_to_hexagonal_face = std::sqrt(1.5) * side;

            // ugly
            if (std::fabs(a.dot(Point(1.0, 0.0, 0.0))) > origin_to_square_face) return true;
            if (std::fabs(a.dot(Point(0.0, 1.0, 0.0))) > origin_to_square_face) return true;
            if (std::fabs(a.dot(Point(0.0, 0.0, 1.0))) > origin_to_square_face) return true;
            if (std::fabs(a.dot(Point(1.0, 1.0, 1.0) * sqrtThreeI)) > origin_to_hexagonal_face) return true;
            if (std::fabs(a.dot(Point(1.0, 1.0, -1.0) * sqrtThreeI)) > origin_to_hexagonal_face) return true;
            if (std::fabs(a.dot(Point(1.0, -1.0, -1.0) * sqrtThreeI)) > origin_to_hexagonal_face) return true;
            if (std::fabs(a.dot(Point(1.0, -1.0, 1.0) * sqrtThreeI)) > origin_to_hexagonal_face) return true;
            return false;
        }

        inline void TruncatedOctahedron::boundary(Point &a) const {
            const double sqrtThreeI = 1.0 / std::sqrt(3.0);
            const double square_face_distance = std::sqrt(8.) * side;
            const double hexagonal_face_distance = std::sqrt(6.) * side;
            const Point unitvXYZ = Point(1, 1, 1) * sqrtThreeI;
            const Point unitvXiYZ = Point(1, 1, -1) * sqrtThreeI;
            const Point unitvXYiZ = Point(1, -1, -1) * sqrtThreeI;
            const Point unitvXYZi = Point(1, -1, 1) * sqrtThreeI;

            // todo improve
            bool outside = false;
            do {
                outside = false;
                double tmp = a.dot(unitvXYZ);
                if (std::fabs(tmp) > hexagonal_face_distance / 2) {
                    a -= hexagonal_face_distance * anint(tmp / hexagonal_face_distance) * unitvXYZ;
                    outside = true;
                }
                tmp = a.dot(unitvXiYZ);
                if (std::fabs(tmp) > hexagonal_face_distance / 2) {
                    a -= hexagonal_face_distance * anint(tmp / hexagonal_face_distance) * unitvXiYZ;
                    outside = true;
                }
                tmp = a.dot(unitvXYiZ);
                if (std::fabs(tmp) > hexagonal_face_distance / 2) {
                    a -= hexagonal_face_distance * anint(tmp / hexagonal_face_distance) * unitvXYiZ;
                    outside = true;
                }
                tmp = a.dot(unitvXYZi);
                if (std::fabs(tmp) > hexagonal_face_distance / 2) {
                    a -= hexagonal_face_distance * anint(tmp / hexagonal_face_distance) * unitvXYZi;
                    outside = true;
                }
            } while (outside);

            if (std::fabs(a.x()) > square_face_distance / 2)
                a.x() -= square_face_distance * anint(a.x() / square_face_distance);

            if (std::fabs(a.y()) > square_face_distance / 2)
                a.y() -= square_face_distance * anint(a.y() / square_face_distance);

            if (std::fabs(a.z()) > square_face_distance / 2)
                a.z() -= square_face_distance * anint(a.z() / square_face_distance);
        }

        void TruncatedOctahedron::randompos(Point &m, Random &rand) const {
            const double d = std::sqrt(10) * side;  // use circumdiameter
            const double r2 = d * d / 4.;
            do {
                do {
                    m.x() = (rand() - 0.5) * d;
                    m.y() = (rand() - 0.5) * d;
                    m.z() = (rand() - 0.5) * d;
                } while (m.squaredNorm() > r2);
            } while (collision(m));
        }

        void TruncatedOctahedron::from_json(const json &j) {
            side = j.at("radius").get<double>();
        }

        void TruncatedOctahedron::to_json(json &j) const {
            j = {{"radius", side}};
        }


        // =============== Chameleon==============

        const std::map<std::string, Variant> Chameleon::names = {
                {
                        {"cuboid", CUBOID},
                        {"cylinder", CYLINDER},
                        {"slit", SLIT},
                        {"sphere", SPHERE},
                        {"hexagonal", HEXAGONAL},
                        {"octahedron", OCTAHEDRON},
                        {"hypersphere2d", HYPERSPHERE2D}
                }
        };

        void from_json(const json &j, Chameleon &g) {
            try { 
                g.from_json(j);
            } catch(std::exception& e) {
                throw std::runtime_error("geometry construction error: "s + e.what() + usageTip["geometry"]);
            }
        }

        void to_json(json &j, const Chameleon &g) {
            g.to_json(j);
        }

        inline Point Chameleon::getLength() const {
            assert(geometry);
            return geometry->getLength();
        }

        inline void Chameleon::_setLength(const Point &l) {
            len = l;
            len_half = l * 0.5;
            len_inv = l.cwiseInverse();
        }

        inline void Chameleon::setLength(const Point &l) {
            assert(geometry);
            _setLength(l);
            // ugly
            if(type == CUBOID) {
                Cuboid &cuboid = dynamic_cast<Cuboid &>(*geometry);
                cuboid.setLength(l);
            } else {
                throw std::runtime_error("setLength allowed only for the Cuboid geometry");
            }
        }

        inline double Chameleon::getVolume(int dim) const {
            assert(geometry);
            return geometry->getVolume(dim);
        }

        inline Point Chameleon::setVolume(double V, VolumeMethod method) {
            auto scale = geometry->setVolume(V, method);
            _setLength(geometry->getLength());
            return scale;
        }

        inline void Chameleon::randompos(Point &m, Random &rand) const {
            assert(geometry);
            geometry->randompos(m, rand);
        }

        inline bool Chameleon::collision(const Point &a) const {
            assert(geometry);
            return geometry->collision(a);
        }

        void Chameleon::makeGeometry(const Variant type) {
            switch (type) {
                case CUBOID:
                    geometry = std::make_unique<Cuboid>();
                    break;
                case SLIT:
                    geometry = std::make_unique<Slit>();
                    break;
                case SPHERE:
                    geometry = std::make_unique<Sphere>();
                    break;
                case CYLINDER:
                    geometry = std::make_unique<Cylinder>();
                    break;
                case HEXAGONAL:
                    geometry = std::make_unique<HexagonalPrism>();
                    break;
                case OCTAHEDRON:
                    geometry = std::make_unique<TruncatedOctahedron>();
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

    } // end of Geometry namespace

} // end of Faunus namespace
