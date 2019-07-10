#pragma once

#include "units.h"
#include "core.h"
#include "particle.h"
#include "tensor.h"
#include <Eigen/Geometry>
#include "spdlog/spdlog.h"

/** @brief Faunus main namespace */
namespace Faunus {

struct Random;

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] spherical coordinates") {
        using doctest::Approx;

        Point sph1 = {2, 0.5, -0.3};
        auto pnt1 = rtp2xyz(sph1); // sph --> cart
        auto sph2 = xyz2rtp(pnt1); // cart --> sph

        CHECK( pnt1.norm() == Approx(2));
        CHECK( sph1.x() == Approx(sph2.x()));
        //CHECK( sph1.y() == Approx(sph2.y()));
        //CHECK( sph1.z() == Approx(sph2.z()));
    }
#endif

    /**
     * @brief Simulation geometries and related operations.
     *
     * Other parts of Faunus use directly only Chameleon geometry which serves as an interface. Based on the provided
     * configuration, Chameleon initializes an appropriate concrete implementation, which it encapsulates.
     *
     * To add a new geometry implementation, a class derived from GeometryImplementation is created. Geometry::Variant
     * enum type is extended and an initialization within Chameleon::makeGeometry() is provided. In order to make
     * geometry constructable from a json configuration, the map Chameleon::names is extended. When performance is
     * an issue, inlineable implementation of vdist and boundary can be added into respective methods of Chameleon.
     *
     * All geometry implementation shall be covered by unit tests.
     *
     */
    namespace Geometry {
        typedef std::function<void(Point&)> BoundaryFunction;
        typedef std::function<Point(const Point&, const Point&)> DistanceFunction;

        //! Geometry variant used for Chameleon.
        enum Variant {CUBOID = 0, SPHERE, CYLINDER, SLIT, HEXAGONAL, OCTAHEDRON, HYPERSPHERE2D};

        //! Various methods of volume scaling, @see GeometryBase::setVolume.
        enum VolumeMethod {ISOTROPIC, ISOCHORIC, XY, Z};

        enum Coordinates {ORTHOGONAL, ORTHOHEXAGONAL, TRUNC_OCTAHEDRAL, NON3D};
        enum Boundary {FIXED, PERIODIC};


        /**
         * @brief A structure containing a type of boundary condition in each direction.
         *
         * A stub. It can be extended to fully json-configurable boundary conditions.
         */
        struct BoundaryCondition {
            typedef Eigen::Matrix<Boundary, 3, 1> BoundaryXYZ;
            //typedef std::pair<std::string, BoundaryXYZ> BoundaryName;
            //static const std::map<std::string, BoundaryXYZ> names; //!< boundary names

            Coordinates coordinates;
            BoundaryXYZ direction;

            BoundaryCondition(Coordinates coordinates = ORTHOGONAL, BoundaryXYZ boundary = {FIXED, FIXED, FIXED}) :
                coordinates(coordinates), direction(boundary) {};
        };


        /**
         * @brief An interface for all geometries.
         */
        struct GeometryBase {
            virtual Point setVolume(double, VolumeMethod = ISOTROPIC) = 0; //!< Set volume
            virtual double getVolume(int = 3) const = 0; //!< Get volume
            virtual void boundary(Point &) const = 0; //!< Apply boundary conditions
            virtual bool collision(const Point &) const = 0; //!< Overlap with boundaries
            virtual void randompos(Point &, Random &) const = 0; //!< Generate random position
            virtual Point vdist(const Point &, const Point &) const = 0; //!< (Minimum) distance between two points
            virtual Point getLength() const = 0; //!< Side lengths

            inline double sqdist( const Point &a, const Point &b ) const {
                return vdist(a,b).squaredNorm();
            } //!< Squared (minimum) distance between two points

            virtual ~GeometryBase();
            virtual void to_json(json &j) const = 0;
            virtual void from_json(const json &j) = 0;

            inline BoundaryFunction getBoundaryFunc() const {
                return [this](Point &i){boundary(i);};
            } //!< returns lambda to boundary()

            inline DistanceFunction getDistanceFunc() const {
                return [this](const Point &i, const Point &j){return vdist(i,j);};
            } //!< returns lambda to vdist()

          protected:
            template<typename T=double>
            inline int anint( T x ) const {
                return int(x > 0.0 ? x + 0.5 : x - 0.5);
            } //!< Round to int

        }; //!< Base class for all geometries


        /**
         * @brief A base class for various geometries implementations.
         */
        class GeometryImplementation : public GeometryBase {
          protected:
            std::shared_ptr<spdlog::logger> faunus_logger, mcloop_logger;

            GeometryImplementation() :
                    faunus_logger(spdlog::get("faunus")),
                    mcloop_logger(spdlog::get("mcloop")) {};

          public:
            BoundaryCondition boundary_conditions;

            virtual ~GeometryImplementation();

            //! A unique pointer to a copy of self. To be used in copy constructors.
            virtual std::unique_ptr<GeometryImplementation> clone() const = 0;
        };


        /**
         * @brief The cuboid geometry with periodic boundary conditions possibly applied in all three directions.
         */
        class Cuboid : public GeometryImplementation {
          protected:
            Point box, box_half, box_inv;

          public:
            Point getLength() const override;
            double getVolume(int dim = 3) const final; // finalized to help the compiler with inlining
            void setLength(const Point &len); // todo shall be protected
            Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            Cuboid(const Point &p);
            Cuboid(double x, double y, double z);
            Cuboid(double x = 0.0);

            std::unique_ptr<GeometryImplementation> clone() const override {
                return std::make_unique<Cuboid>(*this);
            }; //!< A unique pointer to a copy of self.
        };


        /**
         * @brief A legacy class for the cuboid geometry with periodic boundary conditions only in xy directions.
         *
         * @deprecated Shall be replaced by Cuboid with a proper periodic boundary set on initialization.
         */
        class Slit : public Cuboid {
            using Tbase = Cuboid;
          public:
            Slit(const Point &p);
            Slit(double x, double y, double z);
            Slit(double x = 0.0);

            std::unique_ptr<GeometryImplementation> clone() const override {
                return std::make_unique<Slit>(*this);
            }; //!< A unique pointer to a copy of itself.
        };


        /**
         * @brief The spherical geometry where no periodic boundary condition could be applied.
         */
        class Sphere : public GeometryImplementation {
          protected:
            double radius;

          public:
            Point getLength() const override;
            double getVolume(int dim=3) const override;
            Point setVolume(double volume, VolumeMethod method=ISOTROPIC) override;
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            Sphere(double radius = 0.0);

            std::unique_ptr<GeometryImplementation> clone() const override {
                return std::make_unique<Sphere>(*this);
            }; //!< A unique pointer to a copy of self.
        };


        class Hypersphere2d: public Sphere {

          public:
            Point vdist(const Point &a, const Point &b) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            Hypersphere2d(double radius = 0.0);

            std::unique_ptr<GeometryImplementation> clone() const override {
                return std::make_unique<Hypersphere2d>(*this);
            }; //!< A unique pointer to a copy of self.
        };


        /**
         * @brief The cylindrical geometry with periodic boundary conditions in z-axis (the height of the cylinder).
         */
        class Cylinder : public GeometryImplementation {
          protected:
            double radius, height;

          public:
            Point getLength() const override;
            double getVolume(int dim = 3) const override;
            Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            Cylinder(double radius = 0.0, double height = 0.0);

            std::unique_ptr<GeometryImplementation> clone() const override {
                return std::make_unique<Cylinder>(*this);
            }; //!< A unique pointer to a copy of self.
        };


        /**
         * @brief The hexagonal prism geometry with periodic boundary conditions.
         *
         * The prism is oriented in the coordination system as follows: z height, xy base ⬢ with a shorter length
         * (the diameter of an inscribed circle d = 2r) in x direction, and a longer length (the diameter of a
         * circumscribed circle D = 2R) in y direction.
         */
        class HexagonalPrism : public GeometryImplementation {
            //! Change matrices from rhombic to cartesian coordinates and back.
            //! The Y (and Z) axis is identical in both coordination systems, while the X-axis is tilted by -30deg (i.e.,
            //! clockwise), forming a 120deg angle with the Y axis in the rhombic coordinates.
            static const Eigen::Matrix3d rhombic2cartesian;
            static const Eigen::Matrix3d cartesian2rhombic;

            Point box; //!< x = inscribed circle diameter, y = circumscribed circle diameter, z = height
            void set_box(double side, double height);

        public:
            Point getLength() const override;
            double getVolume(int dim=3) const override;
            Point setVolume(double volume, VolumeMethod method=ISOTROPIC) override;
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            HexagonalPrism(double side = 0.0, double height = 0.0);

            std::unique_ptr<GeometryImplementation> clone() const override {
                return std::make_unique<HexagonalPrism>(*this);
            }; //!< A unique pointer to a copy of self.
        };


        /**
         * @brief The truncated octahedron geoemtry with periodic boundary conditions in all directions.
         */
        class TruncatedOctahedron : public GeometryImplementation {
            double side;

          public:
            Point getLength() const override;
            double getVolume(int dim = 3) const override;
            Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            TruncatedOctahedron(double side = 0.0);

            std::unique_ptr<GeometryImplementation> clone() const override {
                return std::make_unique<TruncatedOctahedron>(*this);
            }; //!< A unique pointer to a copy of self.
        };

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Geometry") {
        using doctest::Approx;
        Random slump;

        SUBCASE("cuboid") {
            double x = 2, y = 3, z = 4;
            Cuboid geo(x, y, z);
            CHECK( geo.getVolume() == doctest::Approx(x*y*z) );

            // check boundaries and pbc
            Point a(1.1, 1.5, -2.001);
            CHECK( geo.collision(a) == true );
            geo.getBoundaryFunc()(a);
            CHECK( geo.collision(a) == false );
            CHECK( a.x() == Approx(-0.9) );
            CHECK( a.y() == Approx(1.5) );
            CHECK( a.z() == Approx(1.999) );
            Point b = a;
            geo.boundary(b);
            CHECK( a == b );

            // check distances
            Point distance = geo.vdist({0.1,0.5,-1.001}, a);
            CHECK( distance.x() == Approx(1.0) );
            CHECK( distance.y() == Approx(-1.0) );
            CHECK( distance.z() == Approx(1.0) );
            CHECK( geo.vdist({1,2,3}, a) == geo.getDistanceFunc()({1,2,3},a) );

            // check that geometry is properly enscribed in a cuboid
            Point box = geo.getLength();
            CHECK( box.x() == Approx(x) );
            CHECK( box.y() == Approx(y) );
            CHECK( box.z() == Approx(z) );

            // check random position
            Point c(x+1, y+1, z+1); // out of the box
            bool container_overlap = false;
            for (int i=0; i<1e4; i++) {
                geo.randompos(c, slump);
                if (geo.collision(c))
                    container_overlap = true;
            }
            CHECK( container_overlap == false );

            // volume scaling
            double sf = 2.;
            auto scaling = geo.setVolume(sf*sf*sf*x*y*z);
            CHECK( geo.getVolume() == doctest::Approx(sf*sf*sf*x*y*z) );
            CHECK( geo.getLength().x() == Approx(sf*x) );
            CHECK( geo.getLength().y() == Approx(sf*y) );
            CHECK( geo.getLength().z() == Approx(sf*z) );
            CHECK( scaling.x() == Approx(sf) );
            CHECK( scaling.y() == Approx(sf) );
            CHECK( scaling.z() == Approx(sf) );

            // check json
            geo.from_json(R"( {"type": "cuboid", "length": [2.5,3.5,4.5]} )"_json);
            CHECK( geo.getVolume() == doctest::Approx(2.5*3.5*4.5) );
        }


        SUBCASE("slit") {
            double x = 2, y = 4, z = 3;
            Slit geo(x, y, z);
            CHECK( geo.getVolume() == doctest::Approx(x*y*z) );

            // check boundaries and pbc
            Point a(1.1, -2.001, 1.499);
            CHECK( geo.collision(a) == true);
            geo.getBoundaryFunc()(a);
            CHECK( geo.collision(a) == false);
            CHECK( a.x() == Approx(-0.9) );
            CHECK( a.y() == Approx(1.999) );
            CHECK( a.z() == Approx(1.499) );
            Point b = a;
            geo.boundary(b);
            CHECK( a == b );
            Point c(0, 0, -0.51*z);

            // check distances
            Point distance = geo.vdist({0.1,-1.001,-0.501}, a);
            CHECK( distance.x() == Approx(1.0));
            CHECK( distance.y() == Approx(1.0));
            CHECK( distance.z() == Approx(-2.0));
            CHECK( geo.vdist({1,2,3}, a) == geo.getDistanceFunc()({1,2,3},a) );

            // check that geometry is properly enscribed in a cuboid
            Point box = geo.getLength();
            CHECK( box.x() == Approx(x) );
            CHECK( box.y() == Approx(y) );
            CHECK( box.z() == Approx(z) );

            // check random position
            Point d(x+1, y+1, z+1); // out of the box
            bool containerOverlap = false;
            for (int i=0; i<1e4; i++) {
                geo.randompos(d, slump);
                if (geo.collision(d))
                    containerOverlap = true;
            }
            CHECK( containerOverlap == false );

            // volume scaling
            double sf = 2.;
            auto scaling = geo.setVolume(sf*sf*x*y*z, XY);
            CHECK( geo.getVolume() == doctest::Approx(sf*sf*x*y*z) );
            CHECK( geo.getLength().x() == Approx(sf*x) );
            CHECK( geo.getLength().y() == Approx(sf*y) );
            CHECK( geo.getLength().z() == Approx(z) );
            CHECK( scaling.x() == Approx(sf) );
            CHECK( scaling.y() == Approx(sf) );
            CHECK( scaling.z() == Approx(1.0) );

            // check json
            geo.from_json(R"( {"type": "cuboid", "length": [2.5,3.5,4.5]} )"_json);
            CHECK( geo.getVolume() == doctest::Approx(2.5*3.5*4.5) );
        }


        SUBCASE("sphere") {
            double radius = 5.;
            Sphere geo(radius);
            CHECK( geo.getVolume() == doctest::Approx(4./3.*pc::pi*radius*radius*radius) );

            // check boundaries
            CHECK( geo.collision( { 5.01, 0, 0 } ) == true );
            CHECK( geo.collision( { 4.99, 0, 0 } ) == false );
            Point a(radius - 1, 0, -0.5 * radius);
            Point b = a;
            geo.boundary(a);
            CHECK (a == b);

            // check distances
            Point distance = geo.vdist({3.0,1.0,-2.0}, {-3.0,-1.0,2.0});
            CHECK( distance.x() == Approx(6.0));
            CHECK( distance.y() == Approx(2.0));
            CHECK( distance.z() == Approx(-4.0));

            // check that geometry is properly enscribed in a cuboid
            Point box = geo.getLength();
            CHECK( box.x() == Approx(10) );
            CHECK( box.y() == Approx(10) );
            CHECK( box.z() == Approx(10) );

            // check random position
            Point c(radius+1, radius+1, radius+1); // out of the box
            bool container_overlap = false;
            for (int i=0; i<1e4; i++) {
                geo.randompos(c, slump);
                if (geo.collision(c))
                    container_overlap = true;
            }
            CHECK( container_overlap == false );

            // volume scaling
            geo.setVolume(123.4);
            CHECK( geo.getVolume() == Approx(123.4) );
            CHECK_THROWS_AS( geo.setVolume(100., ISOCHORIC), std::invalid_argument );
            CHECK_THROWS_AS( geo.setVolume(100., XY), std::invalid_argument );

            // check json
            geo.from_json(R"( { "type": "sphere", "radius": 2.0 } )"_json);
            CHECK( geo.getVolume() == doctest::Approx(4./3.*pc::pi*2.0*2.0*2.0) );
        }


        SUBCASE("cylinder") {
            double radius = 1., volume = 1.;
            double height = volume / (pc::pi * radius * radius);
            Point box;
            Cylinder geo(radius, height);

            // check boundaries
            CHECK( geo.getVolume() == Approx( volume ) );
            CHECK( geo.collision({-1.01*radius, 0, 0}) == true );
            CHECK( geo.collision({0.99*radius, 0, 0}) == false );
            CHECK( geo.collision({-0.99*radius, 0.15*radius, 0}) == true );
            CHECK( geo.collision({0, 0, -0.51*height}) == true );
            CHECK( geo.collision({0, 0, 0.49*height}) == false );

            // check that geometry is properly enscribed in a cuboid
            box = geo.getLength();
            CHECK( box.x() == Approx(2*radius) );
            CHECK( box.y() == Approx(2*radius) );
            CHECK( box.z() == Approx(height) );

            // check random position
            Point a(2.*radius, 0, 0);  // out of the box
            bool container_overlap = false;
            for (int i=0; i<1e4; i++) {
                geo.randompos(a, slump);
                if (geo.collision(a))
                    container_overlap = true;
            }
            CHECK( container_overlap == false );

            // volume scaling
            geo.setVolume(9.0, XY);
            CHECK( geo.getVolume() == Approx(9.0) );
            box = geo.getLength();
            CHECK( box.x() == Approx(3*2*radius) );
            CHECK( box.y() == Approx(3*2*radius) );
            CHECK( box.z() == Approx(height) );

            // check json
            json j = {  {"type","cylinder"}, {"radius",2.0}, {"length", 2/pc::pi} };
            geo.from_json(j);
            CHECK( geo.getVolume() == doctest::Approx(8.0) );
        }


        SUBCASE("hexagonal prism") {
            double side = 1., volume = 1.;
            double outer_radius = side, inner_radius = side * std::sqrt(3.0) / 2.;
            double height = volume / (3. * outer_radius * inner_radius);
            Point box;
            HexagonalPrism geo(side, height);

            CHECK( geo.getVolume() == Approx(volume) );
            CHECK( geo.collision({-1.01 * inner_radius, 0, 0}) == true );
            CHECK( geo.collision({0.99 * inner_radius, 0 ,0}) == false );
            CHECK( geo.collision({0.0, -1.01 * outer_radius, 0}) == true );
            CHECK( geo.collision({0.0, 0.99 * outer_radius, 0}) == false );
            CHECK( geo.collision({0.99 * std::cos(pc::pi/3.) * inner_radius, 0.99 * std::sin(pc::pi/3.) * inner_radius, 0}) == false );
            CHECK( geo.collision({1.01 * std::cos(pc::pi/3.) * inner_radius, 1.01 * std::sin(pc::pi/3.) * inner_radius, 0}) == true );
            CHECK( geo.collision({0, 0, -0.51 * height}) == true );
            CHECK( geo.collision({0, 0, 0.49 * height}) == false );

            // check that geometry is properly inscribed in a cuboid
            box = geo.getLength();
            CHECK( box.x() == Approx(2*inner_radius) );
            CHECK( box.y() == Approx(2*outer_radius) );
            CHECK( box.z() == Approx(height) );

            // check random position
            Point a;
            bool container_overlap=false;
            for (int i=0; i<1e4; i++) {
                geo.randompos(a, slump);
                if (geo.collision(a))
                    container_overlap = true;
            }
            CHECK( container_overlap == false );

            json j = {  {"type","hexagonal"}, {"radius", 3*inner_radius}, {"length", 5*height} };
            geo.from_json(j);
            CHECK( geo.getVolume() == Approx(9*5*volume) );
        }

    }
#endif
        /**
         * @brief Geometry class for spheres, cylinders, cuboids, hexagonal prism, truncated octahedron, slits. It is
         * a wrapper of a concrete geometry implementation.
         *
         * The class re-implements the time-critical functions vdist and boundary for the orthogonal periodic boundary
         * conditions. Hence the call can be inlined by the compiler. That would not be possible otherwise due to the
         * polymorphism of the concrete implementations. Other functions calls are delegated directly to the concrete
         * implementation.
         *
         * Note that the class implements a copy constructor and overloads the assignment operator.
         *
         * @note
         * - [Efficient Coding of the Minimum Image Convention](http://doi.org/kvs)
         * - [Fast Coding of the Minimum Image Convention](http://doi.org/ck2nrd)
         *
         * @todo Implement unit tests
         */
        class Chameleon : public GeometryBase {
          private:
            Point len, len_half, len_inv; //!< Cached box dimensions, their half-values, and reciprocal values.
            std::unique_ptr<GeometryImplementation> geometry = nullptr; //!< A concrete geometry implementation.
            Variant _type; //!< Type of concrete geometry.
            std::string _name; //!< Name of concrete geometry, e.g., for json.
            void makeGeometry(const Variant type = CUBOID); //!< Creates and assigns a concrete geometry implementation.
            void _setLength(const Point &l);

          public:
            const Variant& type = _type; //!< Type of concrete geometry, read-only.
            const std::string& name = _name; //!< Name of concrete geometry, e.g., for json, read-only.
            double getVolume(int dim = 3) const override;
            Point setVolume(double V, VolumeMethod method = ISOTROPIC) override;
            Point getLength() const override; //!< A minimal containing cubic box.
            // setLength() needed only for IO::FormatXTC::loadnextframe().
            void setLength(const Point &l); //!< Sets the box dimensions.
            void boundary(Point &a) const override; //!< Apply boundary conditions
            Point vdist(const Point &a, const Point &b) const override; //!< (Minimum) distance between two points
            void randompos(Point &m, Random &rand) const override;
            bool collision(const Point &a) const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;

            static const std::map<std::string, Variant> names; //!< Geometry names.
            typedef std::pair<std::string, Variant> VariantName;

            static VariantName variantName(const std::string &name);

            static VariantName variantName(const json &j);

            Chameleon(const Variant type = CUBOID);
            Chameleon(const GeometryImplementation &geo, const Variant type);

            //! Copy everything, but clone the geometry.
            Chameleon(const Chameleon &geo) : GeometryBase(geo),
                                              len(geo.len), len_half(geo.len_half), len_inv(geo.len_inv),
                                              _type(geo._type), _name(geo._name) {
                geometry = geo.geometry != nullptr ? geo.geometry->clone() : nullptr;
            }

            //! During the assignment copy everything, but clone the geometry.
            Chameleon &operator=(const Chameleon &geo);
        };

        inline void Chameleon::_setLength(const Point &l) {
            len = l;
            len_half = l * 0.5;
            len_inv = l.cwiseInverse();
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

        inline void Chameleon::boundary(Point &a) const {
            const auto &boundary_conditions = geometry->boundary_conditions;
            if (boundary_conditions.coordinates == ORTHOGONAL) {
                if (boundary_conditions.direction.x() == PERIODIC) {
                    if (std::fabs(a.x()) > len_half.x())
                        a.x() -= len.x() * anint(a.x() * len_inv.x());
                }
                if (boundary_conditions.direction.y() == PERIODIC) {
                    if (std::fabs(a.y()) > len_half.y())
                        a.y() -= len.y() * anint(a.y() * len_inv.y());
                }
                if (boundary_conditions.direction.z() == PERIODIC) {
                    if (std::fabs(a.z()) > len_half.y())
                        a.z() -= len.z() * anint(a.z() * len_inv.z());
                }
            } else {
                geometry->boundary(a);
            }
        }

        inline Point Chameleon::vdist(const Point &a, const Point &b) const {
            Point distance;
            const auto &boundary_conditions = geometry->boundary_conditions;
            if (boundary_conditions.coordinates == ORTHOGONAL) {
                distance = a - b;
                if (boundary_conditions.direction.x() == PERIODIC) {
                    if (distance.x() > len_half.x())
                        distance.x() -= len.x();
                    else if (distance.x() < -len_half.x())
                        distance.x() += len.x();
                }
                if (boundary_conditions.direction.y() == PERIODIC) {
                    if (distance.y() > len_half.y())
                        distance.y() -= len.y();
                    else if (distance.y() < -len_half.y())
                        distance.y() += len.y();
                }
                if (boundary_conditions.direction.z() == PERIODIC) {
                    if (distance.z() > len_half.z())
                        distance.z() -= len.z();
                    else if (distance.z() < -len_half.z())
                        distance.z() += len.z();
                }
            } else if (boundary_conditions.coordinates == ORTHOHEXAGONAL ||
                       boundary_conditions.coordinates == TRUNC_OCTAHEDRAL) {
                distance = a - b;
                boundary(distance);
            } else {
                distance = geometry->vdist(a, b);
            }
            return distance;
        }


        void to_json(json&, const Chameleon&);
        void from_json(const json&, Chameleon&);

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Chameleon") {

        using doctest::Approx;
        Random slump;

        //! function compares if Chamelon's and Geometry's boundary methods produce the same result
        //! using n random points
        auto compare_boundary = [&slump](Chameleon &chameleon, GeometryImplementation &geo, Cuboid &box,
                int n = 100) {
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
        auto compare_vdist = [&slump](Chameleon &chameleon, GeometryImplementation &geo, Cuboid &box,
                int n = 100) {
            Point a, b, d_cham, d_geo;
            for (int i = 0; i < n; i++) {
                box.randompos(a, slump);
                box.randompos(b, slump);
                d_cham = chameleon.vdist(a, b);
                d_geo = geo.vdist(a, b);
                CHECK(d_cham.x() == Approx(d_geo.x()));
                CHECK(d_cham.y() == Approx(d_geo.y()));
                CHECK(d_cham.z() == Approx(d_geo.z()));
            }
        };

        SUBCASE("cuboid") {
            double x = 2.0, y = 3.0, z = 4.0;
            Point box_size = std::cbrt(2.0) * Point(x, y, z);
            Cuboid box(box_size);
            Cuboid geo(x, y, z);
            Chameleon chameleon(geo, CUBOID);
            compare_boundary(chameleon, geo, box);
            compare_vdist(chameleon, geo, box);
        }

        SUBCASE("slit") {
            double x = 2.0, y = 3.0, z = 4.0;
            Point box_size = std::cbrt(2.0) * Point(x, y, z);
            Cuboid box(box_size);
            Slit geo(x, y, z);
            Chameleon chameleon(geo, SLIT);
            compare_boundary(chameleon, geo, box);
            compare_vdist(chameleon, geo, box);
        }

        SUBCASE("sphere") {
            double radius = 10.0;
            Point box_size;
            box_size.setConstant(std::cbrt(2.0) * 2*radius);
            Cuboid box(box_size);
            Sphere geo(radius);
            Chameleon chameleon(geo, SPHERE);
            compare_boundary(chameleon, geo, box);
            compare_vdist(chameleon, geo, box);
        }

        SUBCASE("cylinder") {
            double radius = 2.0, height = 10.0;
            Point box_size = std::cbrt(2.0) * Point(2*radius, 2*radius, height);
            Cuboid box(box_size);
            Cylinder geo(radius, height);
            Chameleon chameleon(geo, CYLINDER);
            compare_boundary(chameleon, geo, box);
            compare_vdist(chameleon, geo, box);
        }

        SUBCASE("hexagonal prism") {
            double edge = 5.0, height = 20.0;
            Point box_size = std::cbrt(2.0) * Point(2*edge, 2*edge, height); // a bit larger in x-direction
            Cuboid box(box_size);
            HexagonalPrism geo(edge);
            Chameleon chameleon(geo, HEXAGONAL);
            compare_boundary(chameleon, geo, box);
            compare_vdist(chameleon, geo, box);
        }

        SUBCASE("truncated octahedron") {
            double edge = 5.0;
            Point box_size;
            box_size.setConstant(std::cbrt(2.0) * edge * std::sqrt(5.0/2.0)); // enlarged circumradius
            Cuboid box(box_size);
            TruncatedOctahedron geo(edge);
            Chameleon chameleon(geo, OCTAHEDRON);
            compare_boundary(chameleon, geo, box);
            compare_vdist(chameleon, geo, box);
        }
    }
#endif

        /*
           void unwrap( Point &a, const Point &ref ) const {
           a = vdist(a, ref) + ref;
           } //!< Remove PBC with respect to a reference point
         */
        enum class weight { MASS, CHARGE, GEOMETRIC };

        template<typename Titer, typename Tparticle=typename Titer::value_type, typename weightFunc>
            Point anyCenter( Titer begin, Titer end, BoundaryFunction boundary, const weightFunc &weight,
                    const Point &shift={0,0,0})
            {
                double sum=0;
                Point c(0,0,0), t;
                for (auto &i=begin; i!=end; ++i) {
                    t = i->pos + shift;       // translate to origo
                    boundary(t);
                    double w = weight(*i);
                    c += w * t;
                    sum += w;
                }
                if (sum<=pc::epsilon_dbl)
                    return {0,0,0};
                else
                    c = c/sum - shift;
                boundary(c);
                return c;
            } //!< Mass, charge, or geometric center of a collection of particles

        template<typename Titer, typename Tparticle=typename Titer::value_type>
            Point massCenter(Titer begin, Titer end, BoundaryFunction boundary=[](Point&){}, const Point &shift={0,0,0}) {
                return anyCenter(begin, end, boundary,
                        []( const Tparticle &p ){ return atoms.at(p.id).mw; }, shift );
            } // Mass center

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] anyCenter") {
            Chameleon cyl = json( {{"type","cuboid"}, {"length",100}, {"radius",20}} );
            std::vector<Particle> p;

            CHECK( !atoms.empty() ); // set in a previous test
            p.push_back( atoms[0] );
            p.push_back( atoms[0] );

            p.front().pos = {10, 10, -10};
            p.back().pos  = {15, -10, 10};

            Point cm = Geometry::massCenter(p.begin(), p.end(), cyl.getBoundaryFunc());
            CHECK( cm.x() == doctest::Approx(12.5) );
            CHECK( cm.y() == doctest::Approx(0) );
            CHECK( cm.z() == doctest::Approx(0) );
        }
#endif

        template<class Titer=typename std::vector<T>::iterator>
            void translate( Titer begin, Titer end, const Point &d,
                    BoundaryFunction boundary=[](Point&){}  )
            {
                for ( auto i=begin; i!=end; ++i )
                {
                    i->pos += d;
                    boundary(i->pos);
                }
            } //!< Vector displacement of a range of particles

        template<typename Titer>
            void cm2origo( Titer begin, Titer end, BoundaryFunction boundary=[](Point&){} )
            {
                Point cm = massCenter(begin, end, boundary);
                translate(begin, end, -cm, boundary);
            } //!< Translate to that mass center is in (0,0,0)

        template<typename Titer>
            void rotate(
                    Titer begin,
                    Titer end,
                    const Eigen::Quaterniond &q,
                    BoundaryFunction boundary=[](Point&){},
                    const Point &shift=Point(0,0,0) )
            {
                auto m = q.toRotationMatrix(); // rotation matrix
                for (auto i=begin; i!=end; ++i) {
                    i->rotate(q,m); // rotate internal coordinates
                    i->pos += shift;
                    boundary(i->pos);
                    i->pos = q*i->pos;
                    boundary(i->pos);
                    i->pos -= shift;
                    boundary(i->pos);
                }
            } //!< Rotate particle pos and internal coordinates

        /* 
         * @brief Calculate mass center of cluster of particles in unbounded environment 
         *
         * [More info](http://dx.doi.org/10.1080/2151237X.2008.10129266)
         */
        template<class Tspace, class GroupIndex>
            Point trigoCom( const Tspace &spc, const GroupIndex &groups, const std::vector<int> &dir = {0, 1, 2} )
            {
                assert(!dir.empty() && dir.size() <= 3);
                Point xhi(0, 0, 0), zeta(0, 0, 0), theta(0, 0, 0), com(0, 0, 0);
                for ( auto k : dir ) {
                    double q = 2 * pc::pi / spc.geo.getLength()[k];
                    size_t N=0;
                    for (auto i : groups)
                        for (auto &particle : spc.groups[i]) {
                            theta[k] = particle.pos[k] * q;
                            zeta[k] += std::sin(theta[k]);
                            xhi[k] += std::cos(theta[k]);
                            N++;
                        }
                    theta[k] = std::atan2(-zeta[k] / N, -xhi[k] / N) + pc::pi;
                    com[k] = spc.geo.getLength()[k] * theta[k] / (2 * pc::pi);
                }
                spc.geo.boundary(com); // is this really needed?
                return com;
            }

        /**
         * @brief Calculates the gyration tensor of a molecular group
         * The gyration tensor is computed from the dyadic product of the position
         * vectors in the c-o-m reference system, \f$ t_{i} = r_{i} - shift \f$:
         * \f$ S = \sum_{i=0}^{N} t_{i} t_{i}^{T} \f$
         */
        template<typename iter>
            Tensor gyration(iter begin, iter end, BoundaryFunction boundary=[](const Point&){}, const Point shift=Point(0,0,0) ) {
                Tensor S;
                size_t n = std::distance(begin,end);
                if (n>0) {
                    for (auto it=begin; it!=end; ++it) {
                        Point t = it->pos - shift;
                        boundary(t);
                        // S += t * t.transpose();
                        S += t.squaredNorm() * Eigen::Matrix<double, 3, 3>::Identity() - t * t.transpose();
                    }
                    return S*(1.0/n);
                }
                return S;
            }

    } //geometry namespace
}//end of faunus namespace
