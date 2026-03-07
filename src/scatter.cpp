#include <doctest/doctest.h>
#include "analysis.h"
#include "io.h"

// #define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

namespace Faunus::Scatter {

using doctest::Approx;

TEST_CASE_TEMPLATE("[Faunus] StructureFactorPBC", T,
                   StructureFactorPBC<FormFactorUnity<float>, float, SIMD>,
                   StructureFactorPBC<FormFactorUnity<float>, float, EIGEN>,
                   StructureFactorPBC<FormFactorUnity<float>, float, GENERIC>)
{
    size_t cnt = 0;
    Point box = {80.0, 80.0, 80.0};
    const std::vector<Point> positions = {{10, 20, 30},  {-32, 19, 1},  {34, -2, 23},  {0, 0, 1},
                                          {25, 0, -12},  {-6, -4, -29}, {-12, 23, -3}, {3, 1, -4},
                                          {-31, 29, -20}}; // random position vector

    std::vector<float> result = {0.0785, 1.48621,  0.1111, 0.567279, 0.136,  1.39515,
                                 0.1571, 0.730579, 0.2221, 0.701547, 0.2721, 0.692064};
    T scatter(2);
    scatter.sample(positions, box);
    for (auto [q, S] : scatter.getSampling()) {
        CHECK_EQ(q, Approx(result[cnt++]));
        CHECK_EQ(S, Approx(result[cnt++]));
    }
    CHECK_EQ(cnt, result.size());
}

#ifdef ANKERL_NANOBENCH_H_INCLUDED
TEST_CASE("Benchmark")
{
    Point box = {80.0, 80.0, 80.0};
    std::vector<Point> positions(1000);
    for (auto& position : positions) {
        position = Eigen::Vector3d::Random() * box.x();
    }
    ankerl::nanobench::Bench bench;
    bench.minEpochIterations(100);
    bench.run("SIMD", [&] { StructureFactorPBC<FormFactorUnity<double>, double, SIMD>(10).sample(positions, box); });
    bench.run("EIGEN", [&] { StructureFactorPBC<FormFactorUnity<double>, double, EIGEN>(10).sample(positions, box); });
    bench.run("GENERIC", [&] { StructureFactorPBC<FormFactorUnity<double>, double, GENERIC>(10).sample(positions, box); });
}
#endif

TEST_CASE("[Faunus] StructureFactorIPBC")
{
    size_t cnt = 0;
    Point box = {80.0, 80.0, 80.0};
    const std::vector<Point> positions = {{10, 20, 30},  {-32, 19, 1},  {34, -2, 23},  {0, 0, 1},
                                          {25, 0, -12},  {-6, -4, -29}, {-12, 23, -3}, {3, 1, -4},
                                          {-31, 29, -20}}; // random position vector
    std::vector<double> result = {0.0785, 0.384363, 0.1111, 1.51652, 0.136,  1.18027,
                                  0.1571, 1.40662,  0.2221, 2.06042, 0.2721, 1.53482};
    StructureFactorIPBC scatter(2);
    scatter.sample(positions, box);
    for (auto [q, S] : scatter.getSampling()) {
        CHECK_EQ(q, Approx(result[cnt++]));
        CHECK_EQ(S, Approx(result[cnt++]));
    }
    CHECK_EQ(cnt, result.size());
}

TEST_CASE("[Faunus] FormFactorAtomicConstant")
{
    // Setup atom types with different scattering_f0 values
    Faunus::atoms = R"([
        { "A": { "sigma": 2.0, "scattering_f0": 6.0 } },
        { "B": { "sigma": 3.0, "scattering_f0": 7.5 } }
    ])"_json.get<decltype(Faunus::atoms)>();

    Scatterer s1{{0, 0, 0}, 0}; // type A
    Scatterer s2{{1, 1, 1}, 1}; // type B
    Scatterer s_unity{{2, 2, 2}, -1}; // special id for unity form factor

    FormFactorAtomicConstant<double> ff;
    CHECK_EQ(ff(0.1, s1), Approx(6.0));  // q-independent
    CHECK_EQ(ff(0.5, s1), Approx(6.0));  // same for different q
    CHECK_EQ(ff(0.1, s2), Approx(7.5));  // different atom type
    CHECK_EQ(ff(0.5, s2), Approx(7.5));  // q-independent for type B
    CHECK_EQ(ff(0.1, s_unity), Approx(1.0));  // id=-1 returns unity

    // Verify FormFactorUnity still returns 1
    FormFactorUnity<double> ff_unity;
    CHECK_EQ(ff_unity(0.1, s1), Approx(1.0));
    CHECK_EQ(ff_unity(0.1, s2), Approx(1.0));
}

} // namespace Faunus::Scatter
