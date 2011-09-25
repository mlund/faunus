#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/energy.h>
#include <faunus/potentials.h>
#include <faunus/average.h>
#include <faunus/space.h>
#include <faunus/faunus.h>

#include <typeinfo>

using namespace Faunus;
using namespace Faunus::Geometry;

int main() {
  cout << textio::splash();
  atom.includefile("atomlist.inp");
  inputfile in("cigar.inp");

  particle p,q;
  p.patchangle=1.0;
  p=atom["NA"];

  cuboid geo(in); // simulation geometry
  space spc(geo); // generate space (geometry+particle pool)

  p.sqdist<cuboid>(q);

  //hamiltonian
  Energy::nonbonded< Potential::coulomb_lj<cuboid> > pot_nb(in);
  Energy::hamiltonian pot( &geo ) ;
  pot+=pot_nb;

  cout << atom.info() << spc.info() << pot_nb.pair.info() << endl;
  cout << p << endl;
}
