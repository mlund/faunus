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

    BasePointerVector<Move::Movebase> vec;

    MCSimulation<Tspace> sim(j);
    //sim.move();

    vec.push_back<Move::TranslateRotate<Tspace>>(sim.old.spc, sim.trial.spc);

    FormatPQR::save("confout.pqr", sim.p(), sim.geo().getLength());
}
