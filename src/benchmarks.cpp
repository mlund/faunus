#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>
#include "atomdata.h"
#include "molecule.h"
#include "random.h"
#include "space.h"
#include "analysis.h"
#include "energy.h"
#include "tabulate.h"
#include "aux/arange.h"

int main() {
    using namespace Faunus;
    { // Andrea
        Tabulate::Andrea<double> spline;
        spline.setTolerance(1e-5, 1e-4); // ftol carries no meaning
        double radius = 50;
        auto f = [=](double r) { return std::pow((2.0 / r), 12) - 7.0 / r * std::exp(-r / radius); };
        auto knots = spline.generate(f, 2.0, radius);

        // Distance distribution weighted by volume element, p(r) ∝ 4πr²
        WeightedDistribution<double> distance_distribution;
        for (auto r : arange(2.0 + 0.01, radius - 0.01, 0.5)) {
            distance_distribution.push_back(r, r * r); // register distance and it's relative weight
        }
        // Generate set of weighted distances
        std::vector<double> rvec(10000);
        std::generate(rvec.begin(), rvec.end(), [&]() { return distance_distribution.sample(Faunus::random.engine); });

        ankerl::nanobench::Bench bench;
        bench.minEpochIterations(300).warmup(100).relative(true).performanceCounters(true);
        bench.run("Andrea", [&] {
            double sum = 0;
            for (auto r : rvec) {
                sum += spline.eval(knots, r * r);
            }
            ankerl::nanobench::doNotOptimizeAway(sum);
        });
    }

    { // sofq
        using namespace Faunus::Scatter;
        Point box = {80.0, 80.0, 80.0};
        std::vector<Point> pos(1000);
        for (auto& p : pos)
            p = Eigen::Vector3d::Random() * box.x();
        ankerl::nanobench::Bench bench;
        bench.minEpochIterations(300).warmup(100).relative(true).performanceCounters(true);
        bench.run("S(Q)_SIMD", [&] { StructureFactorPBC<double, SIMD>(10).sample(pos, box); });
        bench.run("S(Q)_EIGEN", [&] { StructureFactorPBC<double, EIGEN>(10).sample(pos, box); });
        bench.run("S(Q)_GENERIC", [&] { StructureFactorPBC<double, GENERIC>(10).sample(pos, box); });
    }

    { // ewald ion-ion
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
            bench.minEpochIterations(20).warmup(100).relative(true).performanceCounters(true);
            bench.run("EwaldIonIon_PBC", [&] { pbc.updateComplex(data, spc.groups); });
            bench.run("EwaldIonIon_PBCEigen", [&] { pbc_eigen.updateComplex(data, spc.groups); });
        }
    }
}
