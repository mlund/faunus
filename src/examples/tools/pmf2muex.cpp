/*
 * Program to calculate activity coefficients from
 * effective PMF's at finite concentrations.
 */

#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main() {
  inputfile in("pmf2muex.conf");

  double conc=in.getflt("pmfconc",0);
  double cellvol=1/(conc*1e-27*6.022e23),
         cellr=pow(3*cellvol/(4*acos(-1.)), 1./3);
  in.add("ionicstr",conc);
  in.add("cellradius", cellr );

  pot_debyehuckel dh(in);
  cell cell(in);
  //pot_table pmf(in);

  particle p1,p2;
  p1=atom("NA");
  p1.charge=+1;
  p1.radius=2.0;
  p2=atom("CL");
  p2.charge=-1;
  p2.radius=2.0;

  //pmf.loadpmf(atom, p2.id, p1.id);
  //pmf.loadpmf(atom, p1.id, p1.id);
  //pmf.loadpmf(atom, p2.id, p2.id);
 
  //cout << coulomb.info() << endl << netz.info() << endl << pmf.info();
  //cout << pmf.info(atom);
  cout << in.info() << dh.info() << cell.info();

  double u, a2=pow(p1.radius+p2.radius,2);
  average<double> expmu;
  for (int i=0; i<1e6; i++) {
    cell.randompos(p2);

    if (p1.sqdist(p2)>a2)
      u=dh.f*dh.pairpot(p1,p2);
    else
      u=1.0e5;
    expmu+=exp(-u);
  }
  cout << "mu     = " << -log(expmu.avg()) << endl
       << "gamma  = " << exp(-log(expmu.avg())) << endl
       << "gamma_dh = "
       << exp(-dh.f*dh.k/2. / (1+dh.k*(p1.radius+p2.radius))     )
       << endl;

  return 0;
 
  cout << "# r/AA  U/kT" << endl;
  for (float r=2.0; r<40; r+=.1) {
    p1.z=r/2;
    p2.z=-p1.z;
    cout << sqrt(p1.sqdist(p2)) << " "
      << dh.f*dh.pairpot(p1,p2) << endl;
  }
}
