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
#include <faunus/io.h>

#include <typeinfo>

using namespace Faunus;

typedef Geometry::cuboid Tgeometry;                // select simulation geometry
typedef Potential::coulomb_lj<Tgeometry> Tpairpot; // select particle-particle pairpotential

int main() {
  atom.includefile("atomlist.inp"); // load atom properties
  inputfile in("bulk.inp");         // read user input
  mcloop loop(in);
  iopqr pqr;

  Tgeometry geo(in); // simulation geometry
  space spc(geo);    // generate simulation space
  Energy::nonbonded<Tpairpot> pot(in);

  atom["NA"].dp=40.;
  atom["CL"].dp=800.;

  spc.insert("NA",100);
  spc.insert("CL",100);
  group salt(0,199);
  salt.name="Salt particles";
  Move::translate_particle mv(in, pot, spc);
  mv.igroup = spc.add(salt);

  pqr.save("init.pqr", spc.p);

  cout << pot.pair.info();

  cout << atom.info() << spc.info() << pot.info() << mv.info() << in.info() << endl;

  double du=0;
  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      du+=mv.move();
    }
    cout << loop.timing();
  }

  for (int i=0; i<spc.p.size(); i++)
    if ( geo.collision(spc.trial[i])==true )
      cout << "!" << endl;

  pqr.save("confout.pqr", spc.p);

  cout << mv.info();
}
