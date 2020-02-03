#include "geometry.h"
#include "units.h"
#include <Eigen/Geometry>
#include <range/v3/view/sample.hpp>
#include <range/v3/view/bounded.hpp>
#include <cereal/archives/binary.hpp>

namespace Faunus {
namespace Geometry {

using doctest::Approx;

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
        Cuboid geo(x, y, z);
        CHECK(geo.getVolume() == doctest::Approx(x * y * z));

        // check boundaries and pbc
        Point a(1.1, 1.5, -2.001);
        CHECK(geo.collision(a) == true);
        geo.getBoundaryFunc()(a);
        CHECK(geo.collision(a) == false);
        CHECK(a.x() == Approx(-0.9));
        CHECK(a.y() == Approx(1.5));
        CHECK(a.z() == Approx(1.999));
        Point b = a;
        geo.boundary(b);
        CHECK(a == b);

        // check distances
        Point distance = geo.vdist({0.1, 0.5, -1.001}, a);
        CHECK(distance.x() == Approx(1.0));
        CHECK(distance.y() == Approx(-1.0));
        CHECK(distance.z() == Approx(1.0));
        CHECK(geo.vdist({1, 2, 3}, a) == geo.getDistanceFunc()({1, 2, 3}, a));

        // check that geometry is properly enscribed in a cuboid
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
        auto scaling = geo.setVolume(sf * sf * x * y * z, XY);
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
        CHECK_THROWS_AS(geo.setVolume(100., ISOCHORIC), std::invalid_argument);
        CHECK_THROWS_AS(geo.setVolume(100., XY), std::invalid_argument);

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
        geo.setVolume(9.0, XY);
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
        box_size.setConstant(std::cbrt(2.0) * 2 * radius);
        Cuboid box(box_size);
        Sphere geo(radius);
        Chameleon chameleon(geo, SPHERE);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
    }

    SUBCASE("cylinder") {
        double radius = 2.0, height = 10.0;
        Point box_size = std::cbrt(2.0) * Point(2 * radius, 2 * radius, height);
        Cuboid box(box_size);
        Cylinder geo(radius, height);
        Chameleon chameleon(geo, CYLINDER);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
    }

    SUBCASE("hexagonal prism") {
        double edge = 5.0, height = 20.0;
        Point box_size = std::cbrt(2.0) * Point(2 * edge, 2 * edge, height); // a bit larger in x-direction
        Cuboid box(box_size);
        HexagonalPrism geo(edge);
        Chameleon chameleon(geo, HEXAGONAL);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
    }

    SUBCASE("truncated octahedron") {
        double edge = 5.0;
        Point box_size;
        box_size.setConstant(std::cbrt(2.0) * edge * std::sqrt(5.0 / 2.0)); // enlarged circumradius
        Cuboid box(box_size);
        TruncatedOctahedron geo(edge);
        Chameleon chameleon(geo, OCTAHEDRON);
        compare_boundary(chameleon, geo, box);
        compare_vdist(chameleon, geo, box);
    }

    SUBCASE("Cereal serialisation") {
        double x = 2.0, y = 3.0, z = 4.0;
        std::ostringstream os(std::stringstream::binary);
        { // write
            Cuboid geo(x, y, z);
            cereal::BinaryOutputArchive archive(os);
            archive(geo);
        }

        { // read
            Cuboid geo(10, 20, 30);
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

} // namespace Geometry
} // namespace Faunus

