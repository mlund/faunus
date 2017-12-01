#include "core.h"
#include "move.h"

using namespace Faunus;
using namespace std;

typedef Particle<Radius, Charge> Tparticle;
typedef Space<Geometry::Cuboid, Tparticle> Tspace;

int main() {

    // run with:
    // ./yason.py < example.yml | ./example
    //json j;
    //std::cin >> j;

    json j = openjson("example.json");


    MCSimulation<Tspace> sim(j);

    for (int i=0; i<10000; i++)
        sim.move();

    FormatPQR::save("confout.pqr", sim.p(), sim.geo().getLength());
}
