#include <doctest/doctest.h>
#include "io.h"
#include "space.h"

namespace Faunus {

#ifdef DOCTEST_LIBRARY_INCLUDED

using doctest::Approx;

TEST_SUITE_BEGIN("IO");

TEST_CASE("[Faunus] FormatPQR") {
    using doctest::Approx;

    // Space object with two salt pairs, i.e. four particles
    Space spc;
    SpaceFactory::makeNaCl(spc, 2, R"( {"type": "cuboid", "length": [20,30,40]} )"_json);
    double d = 0;

    // fix positions
    for (int i = 0; i < 4; i++) {
        spc.p[i].pos = {d, d + 0.1, d + 0.2};
        spc.p[i].charge = double(i);
        d += 0.5;
    }

    // write PQR stream
    std::ostringstream out;
    FormatPQR::save(out, spc.p, spc.geo.getLength());

    // read from stream
    std::stringstream in(out.str());
    spc.p.clear();
    Point box_length = FormatPQR::load(in, spc.p, true);
    CHECK(box_length.x() == 20);
    CHECK(box_length.y() == 30);
    CHECK(box_length.z() == 40);
    CHECK(spc.p.size() == 4);

    // check if positions are restored
    d = 0;
    for (int i = 0; i < 4; i++) {
        spc.p[i].pos = spc.p[i].pos - box_length / 2;
        CHECK(spc.p[i].pos.x() == Approx(d));
        CHECK(spc.p[i].pos.y() == Approx(d + 0.1));
        CHECK(spc.p[i].pos.z() == Approx(d + 0.2));
        CHECK(spc.p[i].charge == double(i));
        d += 0.5;
    }
};

#endif
} // namespace Faunus
