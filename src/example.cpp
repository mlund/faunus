#include "core.h"
#include "move.h"
#include "system.h"

using namespace Faunus;
using namespace std;

typedef Particle<Radius, Dipole, Cigar, Charge> Tparticle;
typedef Space<Geometry::Cuboid, Tparticle> Tspace;

int main() {

    // run with:
    // ./yason.py < example.yml | ./example 
    json j;
    std::cin >> j;

    auto _j = j["moves"];
    for (auto &m : j["moves"])
        for (auto it=m.begin(); it!=m.end(); ++it)
            cout << it.key() << endl; 

    MCSimulation<Tspace> sim(j);
    sim.move();
}
