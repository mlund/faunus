#include "particle.h"

namespace Faunus {

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
