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
  inputfile in;
  cuboid con(in);

  particle p;
  p.patchangle=1.0;
  p=atom["NA"];

  //hamiltonian
  Energy::hamiltonian pot;
  Energy::bonded pot_b;
  Energy::nonbonded<Potential::coulomb_lj_mi> pot_nb(in);
  pot+=pot_b;
  pot+=pot_nb;

  cout << atom.info() << con.info() << pot_nb.pair.info() << endl;
  cout << p << endl;

}
