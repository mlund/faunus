#include <doctest/doctest.h>
#include "analysis.h"
#include "io.h"

//#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

namespace Faunus::Scatter {

using doctest::Approx;

TEST_CASE_TEMPLATE("[Faunus] StructureFactorPBC", T, StructureFactorPBC<float, SIMD>, StructureFactorPBC<float, EIGEN>,
                   StructureFactorPBC<float, GENERIC>) {
    size_t cnt = 0;
    Point box = {80.0, 80.0, 80.0};
    const std::vector<Point> positions = {
        {10, 20, 30},  {-32, 19, 1},  {34, -2, 23}, {0, 0, 1},     {25, 0, -12},
        {-6, -4, -29}, {-12, 23, -3}, {3, 1, -4},   {-31, 29, -20}}; // random position vector

    std::vector<float> result = {0.0785, 1.48621,  0.1111, 0.567279, 0.136,  1.39515,
                                 0.1571, 0.730579, 0.2221, 0.701547, 0.2721, 0.692064};
    T scatter(2);
    scatter.sample(positions, box);
    for (auto [q, S] : scatter.getSampling()) {
        CHECK(q == Approx(result[cnt++]));
        CHECK(S == Approx(result[cnt++]));
    }
    CHECK(cnt == result.size());
}

#ifdef ANKERL_NANOBENCH_H_INCLUDED
TEST_CASE("Benchmark") {
    Point box = {80.0, 80.0, 80.0};
    std::vector<Point> pos(1000);
    for (auto &p : pos)
        p = Eigen::Vector3d::Random() * box.x();
    ankerl::nanobench::Config bench;
    bench.minEpochIterations(100);
    bench.run("SIMD", [&] { StructureFactorPBC<double, SIMD>(10).sample(pos, box); }).doNotOptimizeAway();
    bench.run("EIGEN", [&] { StructureFactorPBC<double, EIGEN>(10).sample(pos, box); }).doNotOptimizeAway();
    bench.run("GENERIC", [&] { StructureFactorPBC<double, GENERIC>(10).sample(pos, box); }).doNotOptimizeAway();
}
#endif

TEST_CASE("[Faunus] StructureFactorIPBC") {
    size_t cnt = 0;
    Point box = {80.0, 80.0, 80.0};
    const std::vector<Point> positions = {
        {10, 20, 30},  {-32, 19, 1},  {34, -2, 23}, {0, 0, 1},     {25, 0, -12},
        {-6, -4, -29}, {-12, 23, -3}, {3, 1, -4},   {-31, 29, -20}}; // random position vector
    std::vector<double> result = {0.0785, 0.384363, 0.1111, 1.51652, 0.136,  1.18027,
                                  0.1571, 1.40662,  0.2221, 2.06042, 0.2721, 1.53482};
    StructureFactorIPBC scatter(2);
    scatter.sample(positions, box);
    for (auto [q, S] : scatter.getSampling()) {
        CHECK(q == Approx(result[cnt++]));
        CHECK(S == Approx(result[cnt++]));
    }
    CHECK(cnt == result.size());
}

} // namespace