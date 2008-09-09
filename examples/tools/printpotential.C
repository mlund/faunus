/*
 * Program to print a pair-potentials
 */

#include "faunus/pot_hydrophobic.h"
#include "faunus/pot_netz.h"

using namespace Faunus;
using namespace std;

int main() {
  particle p1,p2,p3;
  pot_setup cfg;
  cfg.eps=0.3;
  pot_netz netz(cfg);
  pot_coulomb coulomb(cfg);
  //coulomb.eps=coulomb.eps*3.;

  p1.charge=-0;
  p1.radius=2.;
  p1.id=particle::CL;

  p2.radius=3.5;
  p2.hydrophobic=true;

  p3.charge=-0;
  p3.radius=2.;
  p3.id=particle::NA;
 
  cout << coulomb.info() << endl << netz.info();
 
  cout << "# r/AA  U/kT" << endl;
  for (float r=p1.radius+p2.radius-2.; r<15; r+=.1) {
    p1.z=r;
    p3.z=r;
    cout << p1.dist(p2) << " "
      << netz.f*netz.hypairpot(p1,p2, p1.dist(p2) ) << " "
      << coulomb.f*coulomb.pairpot(p1,p2) << endl;
  }
}
