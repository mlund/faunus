#include "core.h"
#include "units.h"
#include "particle.h"
#include <range/v3/view/filter.hpp>

namespace Faunus {

using doctest::Approx;

TEST_SUITE_BEGIN("Core");

TEST_CASE("[Faunus] distance") {
    std::vector<long long int> v = {10, 20, 30, 40, 30};
    auto rng = v | ranges::view::filter([](int i) { return i == 30; });
    CHECK(Faunus::distance(v.begin(), rng.begin()) == 2);
    auto it = rng.begin();
    CHECK(Faunus::distance(v.begin(), ++it) == 4);
}

/*
TEST_CASE("[Faunus] asEigenMatrix") {
    using doctest::Approx;
    typedef Particle<Radius, Charge, Dipole, Cigar> T;
    std::vector<T> v(4);
    v[0].pos.x()=5;
    v[1].pos.y()=10;
    v[2].pos.z()=2;
    auto m = asEigenMatrix(v.begin(), v.end(), &T::pos);

    CHECK( m.cols()==3 );
    CHECK( m.rows()==4 );
    CHECK( m.row(0).x() == 5 );
    CHECK( m.row(1).y() == 10 );
    CHECK( m.row(2).z() == 2 );
    CHECK( m.sum() == 17);
    m.row(0).z()+=0.5;
    CHECK( v[0].pos.z() == Approx(0.5) );

    v[2].charge = 2;
    v[3].charge = -12;
    auto m2 = asEigenVector(v.begin()+1, v.end(), &T::charge);
    CHECK( m2.cols()==1 );
    CHECK( m2.rows()==3 );
    CHECK( m2.col(0).sum() == Approx(-10) );
}*/

TEST_CASE("[Faunus] ranunit") {
    Random r;
    int n = 2e5;
    Point rtp(0, 0, 0);
    for (int i = 0; i < n; i++)
        rtp += xyz2rtp(ranunit(r));
    rtp = rtp / n;
    CHECK(rtp.x() == doctest::Approx(1));
    CHECK(rtp.y() == doctest::Approx(0).epsilon(0.005));          // theta [-pi:pi] --> <theta>=0
    CHECK(rtp.z() == doctest::Approx(pc::pi / 2).epsilon(0.005)); // phi [0:pi] --> <phi>=pi/2
}

TEST_CASE("[Faunus] ranunit_polar") {
    Random r;
    int n = 2e5;
    Point rtp(0, 0, 0);
    for (int i = 0; i < n; i++)
        rtp += xyz2rtp(ranunit_polar(r));
    rtp = rtp / n;
    CHECK(rtp.x() == doctest::Approx(1));
    CHECK(rtp.y() == doctest::Approx(0).epsilon(0.005));          // theta [-pi:pi] --> <theta>=0
    CHECK(rtp.z() == doctest::Approx(pc::pi / 2).epsilon(0.005)); // phi [0:pi] --> <phi>=pi/2
}
TEST_SUITE_END();
} // namespace Faunus
