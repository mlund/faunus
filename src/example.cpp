#include "core.h"
#include "move.h"
#include "system.h"

using namespace Faunus;
using namespace std;

typedef Particle<Radius, Dipole, Cigar, Charge> Tparticle;
typedef Space<Geometry::Cuboid, Tparticle> Tspace;

int main() {

    //Tspace spc;

    json j = openjson("example.json");

    MCSimulation<Tspace> sim(j);
    sim.move();
}
