#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>
#include "atomdata.h"
#include "molecule.h"
#include "random.h"
#include "space.h"
#include "analysis.h"
#include "energy.h"

namespace Faunus {

TEST_CASE("sofq") {
    using namespace Faunus::Scatter;
    Point box = {80.0, 80.0, 80.0};
    std::vector<Point> pos(1000);
    for (auto& p : pos)
        p = Eigen::Vector3d::Random() * box.x();
    ankerl::nanobench::Bench bench;
    bench.minEpochIterations(100);
    bench.run("S(Q)_SIMD", [&] { StructureFactorPBC<double, SIMD>(10).sample(pos, box); });
    bench.run("S(Q)_EIGEN", [&] { StructureFactorPBC<double, EIGEN>(10).sample(pos, box); });
    bench.run("S(Q)_GENERIC", [&] { StructureFactorPBC<double, GENERIC>(10).sample(pos, box); });
}

TEST_CASE("ewald_ionion") {
    using namespace Energy;
    Space spc;
    Faunus::atoms.resize(1);
    Faunus::molecules.resize(1);
    Faunus::molecules.front().id() = 0;
    spc.geometry = R"( {"type": "cuboid", "length": 80} )"_json;
    spc.particles.resize(200);
    for (auto& p : spc.particles) {
        p.id = 0;
        p.charge = 1.0;
        p.pos = (Faunus::random() - 0.5) * spc.geometry.getLength();
    }
    Group g(0, spc.particles.begin(), spc.particles.end());
    spc.groups.push_back(g);

    EwaldData data(R"({
                "epsr": 1.0, "alpha": 0.894427190999916, "epss": 1.0,
                "ncutoff": 11.0, "spherical_sum": true, "cutoff": 9.0})"_json);
    Change c;
    c.everything = true;
    data.policy = EwaldData::PBC;

    {
        PolicyIonIon pbc;
        PolicyIonIonEigen pbc_eigen;
        pbc.updateBox(data, spc.geometry.getLength());
        pbc_eigen.updateBox(data, spc.geometry.getLength());

        ankerl::nanobench::Bench bench;
        bench.minEpochIterations(20);
        bench.run("EwaldIonIon_PBC", [&] { pbc.updateComplex(data, spc.groups); });
        bench.run("EwaldIonIon_PBCEigen", [&] { pbc_eigen.updateComplex(data, spc.groups); });
    }
}
} // namespace Faunus
