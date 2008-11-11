/*
 * Program to print a pair-potentials
 */

#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main() {
  species spc;
  inputfile in("printpotential.conf");
  particle p1,p2,p3;
  pot_netz netz(in);
  pot_coulomb coulomb(in);
  pot_tableMI pmf(in);

  p1.charge=-0;
  p1.radius=2.;
  p1.id=particle::CL;

  p2.radius=3.5;
  p2.hydrophobic=true;
  p2.id=particle::NA;

  p3.charge=-0;
  p3.radius=0.;
  p3.id=particle::NA;

  pmf.loadpmf(spc, p2.id, p1.id);
  pmf.loadpmf(spc, p1.id, p1.id);
  pmf.loadpmf(spc, p2.id, p2.id);
 
  //cout << coulomb.info() << endl << netz.info() << endl << pmf.info();
  cout << pmf.info(spc);
 
  cout << "# r/AA  U/kT" << endl;
  for (float r=2.0; r<40; r+=.1) {
    p1.z=r;
    p3.z=r;
    cout << p1.dist(p2) << " "
      //<< netz.f*netz.hypairpot(p1,p2, p1.dist(p2) ) << " "
      //<< coulomb.f*coulomb.pairpot(p1,p2) << endl << " "
      << pmf.f*pmf.pairpot(p1,p2) << endl;
  }
}
