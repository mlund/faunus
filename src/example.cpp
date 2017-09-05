#include "core.h"
#include "move.h"

using namespace Faunus;
using namespace std;

typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
typedef Space<Geometry::Cuboid, Tparticle> Tspace;

int main() {

    Tspace spc;

    MCSimulation<Tspace> sim;

    sim.move();





}
