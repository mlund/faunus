#include "core.h"
#include "move.h"

using namespace Faunus;
using namespace std;

typedef Geometry::Cuboid Tgeometry;
typedef Particle<Radius, Charge> Tparticle;

int main() {

    // run with:
    // ./yason.py < example.yml | ./example
    //json j;
    //std::cin >> j;

    json j = openjson("example.json");

    MCSimulation<Tgeometry,Tparticle> sim(j);
    Analysis::CombinedAnalysis analysis(j, sim.space(), sim.pot());

    for (int i=0; i<10000; i++) {
        sim.move();
        analysis.sample();
    }

    std::ofstream f("out.json");
    if (f) {
        json j = sim;
        j["analysis"] = analysis;
        f << std::setw(4) << j << endl;
    }

    FormatPQR::save("confout.pqr", sim.particles(), sim.geometry().getLength());
}
