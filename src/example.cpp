#include "core.h"
#include "move.h"
#include "multipole.h"
#include <cstdlib>

using namespace Faunus;
using namespace std;

typedef Geometry::Cuboid Tgeometry;
typedef Particle<Radius, Charge> Tparticle;

int main() {

    try {
        json j;
        std::cin >> j;

        MCSimulation<Tgeometry,Tparticle> sim(j);
        Analysis::CombinedAnalysis analysis(j, sim.space(), sim.pot());

        auto& loop = j.at("mcloop");
        int macro = loop.at("macro");
        int micro = loop.at("micro");

        for (int i=0; i<macro; i++) {
            for (int j=0; j<micro; j++) {
                sim.move();
                analysis.sample();
            }
            cout << "absolute drift (kT) = " << sim.drift() << endl;
        }

        std::ofstream f("out.json");
        if (f) {
            json j = sim;
            j["analysis"] = analysis;
            f << std::setw(4) << j << endl;
        }

    } catch (std::exception &e) {
        std::cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
