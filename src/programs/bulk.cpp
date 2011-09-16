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

typedef Geometry::sphere Tgeometry;                // select simulation geometry
typedef Potential::coulomb_lj<Tgeometry> Tpairpot; // select particle-particle pairpotential

int main() {
  atom.includefile("atomlist.inp");    // load atom properties
  inputfile in("bulk.inp");            // read user input
  mcloop loop(in);                     // class for handling mc loops
  iopqr pqr;                           // PQR structure file I/O
  energydrift sys;                     // class for tracking system energy drifts

  Tgeometry geo(in);                   // simulation geometry
  space spc(geo);                      // generate simulation space
  Energy::nonbonded<Tpairpot> pot(in); // energy calculation class

  // Handle particles
  atom["NA"].dp=80.;                   // Displacement parameter
  atom["CL"].dp=80.;
  spc.insert("NA",100);                // Insert particles into space
  spc.insert("CL",100);
  group salt(0,199);                   // Define salt range
  salt.name="Salt particles";
  Move::translate_particle mv(in, pot, spc);  // Particle move class
  mv.igroup = spc.enroll(salt);               // Enroll salt in space and select it for particle moves

  spc.load("confout.spc");
  sys.init( pot.all2all(spc.p) );

  cout << atom.info() << spc.info() << pot.info() << mv.info() << pot.pair.info() << endl;

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move();
    }
    sys.checkdrift( pot.all2all(spc.p) );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("confout.spc");
  cout << mv.info() << sys.info();
}
