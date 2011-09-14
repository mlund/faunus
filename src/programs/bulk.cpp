#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/energy.h>
#include <faunus/potentials.h>
#include <faunus/average.h>
#include <faunus/space.h>
#include <faunus/move.h>
#include <faunus/mcloop.h>
#include <faunus/group.h>

#include <typeinfo>

using namespace Faunus;

typedef Geometry::cuboid Tgeometry;                // select simulation geometry
typedef Potential::coulomb_lj<Tgeometry> Tpairpot; // select particle-particle pairpotential

int main() {
  atom.includefile("atomlist.inp"); // load atom properties
  inputfile in("bulk.inp");         // read user input
  mcloop loop(in);

  Tgeometry geo(in); // simulation geometry
  space spc(geo);    // generate simulation space
  Energy::nonbonded<Tpairpot> pot(in);

  atom["NA"].dp=40.;
  atom["CL"].dp=40.;

  spc.insert("NA",10);
  spc.insert("CL",10);
  group salt(0,19);
  salt.name="Salt particles";
  Move::translate_particle mv(in, pot, spc);
  mv.igroup = spc.add(salt);

  particle a,b;
  cout << pot.pair.info();


  cout << atom.info() << spc.info() << pot.info() << mv.info() << in.info() << endl;

  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      mv.move();
    }
    cout << loop.timing();
  }

  cout << mv.info();
}
