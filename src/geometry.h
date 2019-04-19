#pragma once

#include "core.h"
#include "atomdata.h"
#include "particle.h"

/** @brief Faunus main namespace */
namespace Faunus {

    /**
     * @brief Convert cartesian- to cylindrical-coordinates
     * @note Input (x,y,z), output \f$ (r,\theta, h) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [-\pi,\pi) \f$, and \f$ h\in (-\infty,\infty] \f$.
     */
    Point xyz2rth(const Point&, const Point &origin={0,0,0}, const Point &dir={0,0,1}, const Point &dir2={1,0,0});

    /**
     * @brief Convert cartesian- to spherical-coordinates
     * @note Input (x,y,z), output \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [-\pi,\pi) \f$, and \f$ \phi\in [0,\pi] \f$.
     */
    Point xyz2rtp(const Point&, const Point &origin={0,0,0});

    /**
     * @brief Convert spherical- to cartesian-coordinates
     * @param origin The origin to be added (optional)
     * @note Input \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [0,2\pi) \f$, and \f$ \phi\in [0,\pi] \f$, and output (x,y,z).
     */
    Point rtp2xyz(const Point &rtp, const Point &origin = {0,0,0});

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

    Point ranunit(Random&, const Point &dir={1,1,1}); //!< Random unit vector using Neuman's method ("sphere picking")
    Point ranunit_polar(Random&); //!< Random unit vector using polar coordinates ("sphere picking")

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] ranunit") {
        Random r;
        int n=2e5;
        Point rtp(0,0,0);
        for (int i=0; i<n; i++)
            rtp += xyz2rtp( ranunit(r) );
        rtp = rtp / n;
        CHECK( rtp.x() == doctest::Approx(1) );
        CHECK( rtp.y() == doctest::Approx(0).epsilon(0.005) ); // theta [-pi:pi] --> <theta>=0
        CHECK( rtp.z() == doctest::Approx(pc::pi/2).epsilon(0.005) );// phi [0:pi] --> <phi>=pi/2
    }

    TEST_CASE("[Faunus] ranunit_polar") {
        Random r;
        int n=2e5;
        Point rtp(0,0,0);
        for (int i=0; i<n; i++)
            rtp += xyz2rtp( ranunit_polar(r) );
        rtp = rtp / n;
        CHECK( rtp.x() == doctest::Approx(1) );
        CHECK( rtp.y() == doctest::Approx(0).epsilon(0.005) ); // theta [-pi:pi] --> <theta>=0
        CHECK( rtp.z() == doctest::Approx(pc::pi/2).epsilon(0.005) );// phi [0:pi] --> <phi>=pi/2
    }
#endif

    /** @brief Simulation geometries and related operations */
    namespace Geometry {
        typedef std::function<void(Point&)> BoundaryFunction;
        typedef std::function<Point(const Point&, const Point&)> DistanceFunction;

        enum Variant {CUBOID = 0, SPHERE, CYLINDER, SLIT, HEXAGONAL, OCTAHEDRON};
        enum VolumeMethod {ISOTROPIC, ISOCHORIC, XY, Z};

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

            virtual std::unique_ptr<GeometryBase> clone() const = 0; //!< To be used in copy constructors
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


        class Cuboid : public GeometryBase {
          protected:
            Point box, box_half, box_inv;

          public:
            Point getLength() const override;
            double getVolume(int dim = 3) const override final; // finalized to help the compiler with inlining
            void setLength(const Point &len); // todo shall be protected
            Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            void from_json(const json &j);
            void to_json(json &j) const;
            Cuboid(const Point &p);
            Cuboid(double x, double y, double z);
            Cuboid(double x = 0.0);

            std::unique_ptr<GeometryBase> clone() const override {
                return std::make_unique<Cuboid>(*this);
            };
        };


        class Slit : public Cuboid {
            using Tbase = Cuboid;
          public:
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            Slit(const Point &p);
            Slit(double x, double y, double z);
            Slit(double x = 0.0);

            std::unique_ptr<GeometryBase> clone() const override {
                return std::make_unique<Slit>(*this);
            };
        };


        class Sphere : public GeometryBase {
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
            void from_json(const json &j);
            void to_json(json &j) const;
            Sphere(double radius = 0.0);

            std::unique_ptr<GeometryBase> clone() const override {
                return std::make_unique<Sphere>(*this);
            };
        };


        class Cylinder : public GeometryBase {
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
            void from_json(const json &j);
            void to_json(json &j) const;
            Cylinder(double radius = 0.0, double height = 0.0);

            std::unique_ptr<GeometryBase> clone() const override {
                return std::make_unique<Cylinder>(*this);
            };
        };


        /**
         * The orientation in the coordinate system ⬢ (x→, y↑).
         */
        class HexagonalPrism : public GeometryBase {
            // Change matrices from rhombic to cartesian coordinates and back.
            // Y (and Z) axis is identical, X-axis is tilted by -30 deg (i.e., clockwise)
            static const Eigen::Matrix3d rhombic2cartesian;
            static const Eigen::Matrix3d cartesian2rhombic;
            Point box; //< x = inscribed circle diameter, y = circumscribed circle diameter, z = height
            double volume; //< volume of the prism
            void set_box(double side, double height);

        public:
            Point getLength() const override;
            double getVolume(int dim=3) const override;
            Point setVolume(double volume, VolumeMethod method=ISOTROPIC) override;
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            void from_json(const json &j);
            void to_json(json &j) const;
            HexagonalPrism(double side = 0.0, double height = 0.0);

            std::unique_ptr<GeometryBase> clone() const override {
                return std::make_unique<HexagonalPrism>(*this);
            };
        };

        class TruncatedOctahedron : public GeometryBase {
            double side;

          public:
            Point getLength() const override;
            double getVolume(int dim = 3) const override;
            Point setVolume(double volume, VolumeMethod method = ISOTROPIC) override;
            Point vdist(const Point &a, const Point &b) const override;
            void boundary(Point &a) const override;
            bool collision(const Point &a) const override;
            void randompos(Point &m, Random &rand) const override;
            void from_json(const json &j);
            void to_json(json &j) const;
            TruncatedOctahedron(double side = 0.0);

            std::unique_ptr<GeometryBase> clone() const override {
                return std::make_unique<TruncatedOctahedron>(*this);
            };
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
         * @brief Geometry class for spheres, cylinders, cuboids, hexagonal prism, truncated octahedron, slits
         *
         * All geometries shares the same distance calculation where
         * PBC is disabled by artificially setting long sidelengths.
         *
         * @note
         * - [Efficient Coding of the Minimum Image Convention](http://doi.org/kvs)
         * - [Fast Coding of the Minimum Image Convention](http://doi.org/ck2nrd)
         *
         * @todo Implement unit tests
         */
        class Chameleon : public GeometryBase {
          private:
            Point len, len_half, len_inv;
            std::unique_ptr<GeometryBase> geometry = nullptr;
            Variant _type; //!< type of contained geometry
            std::string _name; //!< name of contained geometry, e.g., for json
            void makeGeometry(const Variant type = CUBOID);
            void _setLength(const Point &l);

          public:
            const Variant& type = _type; //!< type of contained geometry, read-only
            const std::string& name = _name; //!< name of contained geometry, e.g., for json, read-only
            double getVolume(int dim = 3) const override;
            Point setVolume(double V, VolumeMethod method = ISOTROPIC) override;
            Point getLength() const override; //!< Enscribed box
            void setLength(const Point &l);
            void randompos(Point &m, Random &rand) const override;
            bool collision(const Point &a) const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;

            static const std::map<std::string, Variant> names; //!< geometry names
            typedef std::pair<std::string, Variant> VariantName;

            inline static VariantName variantName(const std::string &name) {
                auto it = names.find(name);
                if (it == names.end()) {
                    throw std::runtime_error("unknown geometry: " + name);
                }
                return *it;
            }

            inline static VariantName variantName(const json &j) {
                return variantName(j.at("type").get<std::string>());
            }

            std::unique_ptr<GeometryBase> clone() const override {
                return std::make_unique<Chameleon>(*this);
            };

            explicit Chameleon() {}

            Chameleon(const Chameleon &geo) : GeometryBase(geo),
                                              len(geo.len), len_half(geo.len_half), len_inv(geo.len_inv),
                                              _type(geo._type), _name(geo._name) {
                geometry = geo.geometry != nullptr ? geo.geometry->clone() : nullptr;
            }

            Chameleon &operator=(const Chameleon &geo) {
                if (&geo != this) {
                    GeometryBase::operator=(geo);
                    len = geo.len;
                    len_half = geo.len_half;
                    len_inv = geo.len_inv;
                    _type = geo._type;
                    _name = geo._name;
                    geometry = geo.geometry != nullptr ? geo.geometry->clone() : nullptr;
                }
                return *this;
            }

                inline void boundary( Point &a ) const override {
                    // TODO refactor
                    if ( type!=HEXAGONAL and type!=OCTAHEDRON and type!=SPHERE) {
                        if ( std::fabs(a.x()) > len_half.x())
                            a.x() -= len.x() * anint(a.x() * len_inv.x());
                        if ( std::fabs(a.y()) > len_half.y())
                            a.y() -= len.y() * anint(a.y() * len_inv.y());
                        if ( type != SLIT ) {
                            if ( std::fabs(a.z()) > len_half.z())
                                a.z() -= len.z() * anint(a.z() * len_inv.z());
                        }
                    } else if ( type == HEXAGONAL ) {
                        const double sqrtThreeByTwo = 0.5*sqrt(3.0);
                        const Point unitvX = {1,0,0};
                        const Point unitvY = {0.5,sqrtThreeByTwo,0};
                        const Point unitvZ = {-0.5,sqrtThreeByTwo,0.0};

                        double tmp = a.dot(unitvX);
                        if (std::fabs(tmp) > len_half.x())
                            a -= len.x() * anint(tmp * len_inv.x())*unitvX;

                        if (a.dot(unitvY) > len_half.x()) {
                            a = a - len.x()*unitvY;
                            if (a.dot(unitvX) < -len_half.x()) // Check that point did not get past x-limit
                                a = a + len.x()*unitvX;
                        }
                        if (a.dot(unitvY) < -len_half.x()) {
                            a = a + len.x()*unitvY;
                            if (a.dot(unitvX) > len_half.x()) // Check that point did not get past x-limit
                                a = a - len.x()*unitvX;
                        }

                        tmp = a.dot(unitvZ);
                        if (std::fabs(tmp) > len_half.x())
                            a -= len.x() * anint(tmp * len_inv.x())*unitvZ;
                        if (std::fabs(a.z()) > len_half.z())
                            a.z() -= len.z() * anint(a.z() * len_inv.z());

                    } else if ( type == OCTAHEDRON ) {
                        const double sqrtThreeI = 1.0/std::sqrt(3.0);
                        const Point unitvXYZ  = Point(1, 1, 1)*sqrtThreeI;
                        const Point unitvXiYZ = Point(1, 1,-1)*sqrtThreeI;
                        const Point unitvXYiZ = Point(1,-1,-1)*sqrtThreeI;
                        const Point unitvXYZi = Point(1,-1, 1)*sqrtThreeI;

                        bool outside = false;
                        do {
                            outside = false;
                            double tmp = a.dot(unitvXYZ);
                            if ( std::fabs(tmp) > len_half.x()) {
                                a -= len.x() * anint(tmp * len_inv.x()) * unitvXYZ;
                                outside = true;
                            }
                            tmp = a.dot(unitvXiYZ);
                            if ( std::fabs(tmp) > len_half.x()) {
                                a -= len.x() * anint(tmp * len_inv.x()) * unitvXiYZ;
                                outside = true;
                            }
                            tmp = a.dot(unitvXYiZ);
                            if ( std::fabs(tmp) > len_half.x()) {
                                a -= len.x() * anint(tmp * len_inv.x()) * unitvXYiZ;
                                outside = true;
                            }
                            tmp = a.dot(unitvXYZi);
                            if ( std::fabs(tmp) > len_half.x()) {
                                a -= len.x() * anint(tmp * len_inv.x()) * unitvXYZi;
                                outside = true;
                            }
                        } while (outside);

                        if ( std::fabs(a.x()) > len_half.y())
                            a.x() -= len.y() * anint(a.x() * len_inv.y());

                        if ( std::fabs(a.y()) > len_half.y())
                            a.y() -= len.y() * anint(a.y() * len_inv.y());

                        if ( std::fabs(a.z()) > len_half.y())
                            a.z() -= len.y() * anint(a.z() * len_inv.y());
                    }
                } //!< Apply boundary conditions

                inline Point vdist(const Point &a, const Point &b) const override {
                    Point r(a - b);
                    if (type==CUBOID or type==SLIT) {
                        if (r.x() > len_half.x())
                            r.x() -= len.x();
                        else if (r.x() < -len_half.x())
                            r.x() += len.x();

                        if (r.y() > len_half.y())
                            r.y() -= len.y();
                        else if (r.y() < -len_half.y())
                            r.y() += len.y();

                        if (type == CUBOID) {
                            if (r.z() > len_half.z())
                                r.z() -= len.z();
                            else if (r.z() < -len_half.z())
                                r.z() += len.z();
                        }
                    } else if (type==SPHERE) {
                        // do nothing
                    } else
                        boundary(r);
                    return r;
                }
        };

        void to_json(json&, const Chameleon&);
        void from_json(const json&, Chameleon&);

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] Chameleon") {

            using doctest::Approx;
            Chameleon geo;
            Random slump;

            SUBCASE("cuboid") {
                geo = R"( {"type": "cuboid", "length": [2,3,4]} )"_json;
                CHECK( geo.getVolume() == doctest::Approx(2*3*4) );
                Point a(1.1, 1.5, -2.001);
                CHECK( geo.collision(a) == true);
                geo.getBoundaryFunc()(a);
                CHECK( geo.collision(a) == false);
                CHECK( a.x() == Approx(-0.9) );
                CHECK( a.y() == Approx(1.5) );
                CHECK( a.z() == Approx(1.999) );
                CHECK( geo.vdist({1,2,3}, a) == geo.getDistanceFunc()({1,2,3},a) );
                Point b = a;
                geo.boundary(b);
                CHECK( a == b );

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(2) );
                CHECK( L.y() == Approx(3) );
                CHECK( L.z() == Approx(4) );

                // check random position
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }

            SUBCASE("slit") {
                geo = R"( {"type": "slit", "length": [2,3,4]} )"_json;
                CHECK( geo.getVolume() == doctest::Approx(2*3*4) );
                Point a(1.1, 1.5, -2.001);
                CHECK( geo.collision(a) == true);
                geo.getBoundaryFunc()(a);
                CHECK( geo.collision(a) == true);
                CHECK( a.x() == doctest::Approx(-0.9) );
                CHECK( a.y() == doctest::Approx(1.5) );
                CHECK( a.z() == doctest::Approx(-2.001) );

                Point b = a;
                b.z() = -b.z();
                geo.boundary(b);
                CHECK( b.z() == doctest::Approx(2.001) );
                CHECK( geo.vdist(a,b).x()==0);
                CHECK( geo.vdist(a,b).y()==0);
                CHECK( geo.vdist(a,b).z()==-4.002);

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(2) );
                CHECK( L.y() == Approx(3) );
                CHECK( L.z() == Approx(4) );

                // check random position
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }

            SUBCASE("cylinder") {
                json j = {  {"type","cylinder"}, {"radius",1.0}, {"length", 1/pc::pi} };
                geo = j;
                CHECK( geo.getVolume() == doctest::Approx( 1.0 ) );
                CHECK( geo.collision({1.01,0,0}) == true);
                CHECK( geo.collision({0.99,0,0}) == false);
                CHECK( geo.collision({0.99,1/pc::pi+0.01,0}) == true);

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(2) );
                CHECK( L.y() == Approx(2) );
                CHECK( L.z() == Approx(1/pc::pi) );

                // check random position
                Point a;
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }

            SUBCASE("hexagonal") {
                json j = {  {"type","hexagonal"}, {"radius",1.0}, {"length", 1.0/2.0/std::sqrt(3.0)} };
                geo = j;
                CHECK( geo.getVolume() == doctest::Approx( 1.0 ) );
                CHECK( geo.collision({1.01,0,0}) == true);
                CHECK( geo.collision({0.99,0,0}) == false);
                CHECK( geo.collision({0.0,2.0/std::sqrt(3.0)+0.01,0}) == true);

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                //CHECK( L.x() == Approx(1.0) );
                //CHECK( L.y() == Approx(1.0) );
                CHECK( L.z() == Approx(1.0/2.0/std::sqrt(3.0)) );

                // check random position
                Point a;
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }

            SUBCASE("sphere") {
                geo = R"( { "type": "sphere", "radius": 2 } )"_json;
                CHECK( geo.getVolume() == doctest::Approx( 4*pc::pi*8/3.0 ) );
                CHECK( geo.collision( { 2.01, 0, 0 } ) == true );
                CHECK( geo.collision( { 1.99, 0, 0 } ) == false );
                CHECK( geo.getLength().squaredNorm() == doctest::Approx(48) );

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(4) );
                CHECK( L.y() == Approx(4) );
                CHECK( L.z() == Approx(4) );

                geo.setVolume(123.4);
                CHECK( geo.getVolume() == doctest::Approx(123.4) );

                // check random position
                Point a;
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
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
            typedef Particle<Radius, Charge, Dipole, Cigar> T;

            Chameleon cyl = json( {{"type","cuboid"}, {"length",100}, {"radius",20}} );
            std::vector<T> p;

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

        template<class Titer>
            double monopoleMoment( Titer begin, Titer end ) {
                double z=0;
                for (auto it=begin; it!=end; ++it)
                    z += it->charge;
                return z;
            } //!< Calculates dipole moment vector

        template<class Titer>
            Point dipoleMoment( Titer begin, Titer end, BoundaryFunction boundary=[](const Point&){}, double cutoff=pc::infty) {
                Point mu(0,0,0);
                for (auto it=begin; it!=end; ++it) {
                    Point t = it->pos - begin->pos;
                    boundary(t);
                    if (t.squaredNorm() < cutoff*cutoff)
                        mu += t * it->charge;
                }
                return mu;
            } //!< Calculates dipole moment vector

        template<class Titer>
            Tensor quadrupoleMoment( Titer begin, Titer end, BoundaryFunction boundary=[](const Point&){}, double cutoff=pc::infty) {
                Tensor theta;
                theta.setZero();
                for (auto it=begin; it!=end; ++it) {
                    Point t = it->pos - begin->pos;
                    boundary(t);
                    if (t.squaredNorm() < cutoff*cutoff)
                        theta += t * t.transpose() * it->charge;
                }
                return 0.5 * theta;
            } //!< Calculates quadrupole moment tensor (with trace)

        template<class Tgroup>
            auto toMultipole(const Tgroup &g, BoundaryFunction boundary=[](const Point&){}, double cutoff=pc::infty) {
                Particle<Charge,Dipole,Quadrupole> m;
                m.pos = g.cm;
                m.charge = Geometry::monopoleMoment(g.begin(), g.end());                   // monopole
                m.mu = Geometry::dipoleMoment(g.begin(), g.end(), boundary, cutoff);   // dipole
                m.Q = Geometry::quadrupoleMoment(g.begin(), g.end(), boundary, cutoff);// quadrupole
                m.mulen = m.mu.norm();
                if (m.mulen>1e-9)
                    m.mu.normalize();
                return m;
            } //<! Group --> Multipole

    } //geometry namespace
}//end of faunus namespace
