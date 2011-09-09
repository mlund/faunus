#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/container.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/energy_base.h>
#include <faunus/pot_base.h>

using namespace Faunus;

int main() {

  atom.includefile("atomlist.inp");
  inputfile in("cigar.inp");

  particle p;
  p.patchangle=1.0;
  p=atom["NA"];

  //hamiltonian
  Energy::bonded pot_b;
  Energy::nonbonded<Potential::coulomb_lj_cuboid> pot_nb(in);
  Energy::hamiltonian pot( &pot_nb.pair.geo ) ;
  pot+=pot_b;
  pot+=pot_nb;

  space spc;
  spc.geo=pot.geo;

  cout << atom.info() << spc.info() << pot_nb.pair.info() << endl;
  cout << p << endl;

}
