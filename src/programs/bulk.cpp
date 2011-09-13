#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/energy.h>
#include <faunus/potentials.h>
#include <faunus/average.h>
#include <faunus/space.h>

#include <typeinfo>

using namespace Faunus;

typedef Geometry::cuboid Tgeometry;                // select simulation geometry
typedef Potential::coulomb_lj<Tgeometry> Tpairpot; // select particle-particle pairpotential

int main() {
  atom.includefile("atomlist.inp"); // load atom properties
  inputfile in("bulk.inp");         // read user input

  Tgeometry geo(in); // simulation geometry
  space spc(geo);    // generate simulation space

  Energy::nonbonded<Tpairpot> pot(in);

  cout << atom.info() << spc.info() << pot.info() << endl;
}
