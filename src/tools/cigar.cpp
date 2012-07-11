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

typedef Geometry::Cuboid Tgeometry;
typedef Potential::DebyeHuckel Tpairpot;

int main() {
  cout << textio::splash();
  atom.includefile("atomlist.inp");
  InputMap in("cigar.inp");

  particle p,q;
  p.patchangle=1.0;
  p=atom["Na"];


  //hamiltonian
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(in) );

  Space spc( pot.getGeometry() ); // generate space (geometry+particle pool)

  cout << pot.info();

  //cout << atom.info() << spc.info() << pot_nb.pair.info() << endl;
  cout << p << endl;
}
