#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/container.h>
#include <faunus/species.h>
#include <faunus/inputfile.h>
#include <faunus/energy_base.h>
#include <faunus/pot_base.h>
#include <faunus/average.h>
#include <faunus/space.h>

#include <typeinfo>

using namespace Faunus;
using namespace Faunus::Geometry;

class systemenergy {
  private:
    double delta;
    double initial;
    double drift;
    average<double> total;
  public:
    systemenergy& operator+=(double u) {
      delta+=u;
      total+=u;
      return *this;
    }
    double calc_drift(double current) {
      drift=current-(initial+delta);
      return drift;
    }
};

int main() {
  atom.includefile("atomlist.inp");
  inputfile in("cigar.inp");

  particle p,q;
  p.patchangle=1.0;
  p=atom["NA"];

  cuboid geo(in); // simulation geometry
  space spc(geo); // generate space (geometry+particle pool)

  p.sqdist<cuboid>(q);

  //hamiltonian
  Energy::bonded pot_b;
  Energy::nonbonded< Potential::coulomb_lj<cuboid> > pot_nb(in);
  Energy::hamiltonian pot( &geo ) ;
  pot+=pot_b;
  pot+=pot_nb;

  cout << atom.info() << spc.info() << pot_nb.pair.info() << endl;
  cout << p << endl;

  //typedef typeid(p).name() hej;
  //cout << typeid(p).name() << endl;

}
