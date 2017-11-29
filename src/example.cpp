#include "core.h"
#include "move.h"

using namespace Faunus;
using namespace std;

typedef Particle<Radius, Dipole, Cigar, Charge> Tparticle;
typedef Space<Geometry::Cuboid, Tparticle> Tspace;

int main() {

    // run with:
    // ./yason.py < example.yml | ./example
    json j;
    std::cin >> j;

    MCSimulation<Tspace> sim(j);
    sim.move();
}
