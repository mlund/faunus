#include "core.h"
#include "move.h"

using namespace Faunus;
using namespace std;

typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
typedef Space<Geometry::Cuboid, Tparticle> Tspace;

int main() {

    Tspace spc;

    json j = openjson("example.json");

    MCSimulation<Tspace> sim(j);
    sim.move();
}
