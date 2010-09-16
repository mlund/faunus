#include "faunus/faunus.h"
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/potentials/pot_debyehuckelP3.h"
#include "faunus/bottles/npt_molecular.h"

/*!
 * \brief Replica exchange program - proof of concept.
 * \author Mikael Lund
 * \date Malmo 2010
 */

using namespace std;
using namespace Faunus;

typedef interaction<pot_debyehuckelP3> Tpot;

int main() {
  cout << faunus_splash();
  replicaexchange temper;
  npt_molecular<box,Tpot> a("a"), b("b"), c("c");

  a.prepare();
  b.prepare();
  c.prepare();

  for (int i=0; i<10; i++) {
    cout << "m" << i+1 << " " << flush;
    for (int j=0; j<1e3; j++) {
      a.microloop();
      b.microloop();
      c.microloop();
      if (slp.random_one()>0.95) {
        temper.swap(a,b);
        temper.swap(b,c);
      }
    }
    a.macroloop();
    b.macroloop();
    c.macroloop();
  }
  a.save();
  b.save();
  c.save();

  cout << a.postinfo() << b.postinfo() << temper.info();
};
