/*
 * Program to create a gromacs .gro file for
 * a "star" molecule (a central atom to which
 * nchains is attached, each containing nmon
 * monomers). The monomers are generated to
 * be as close to the previous monomer as possible.
 * 
 * M. Lund, 2007
 */

#include <iostream>
#include <vector>
#include "point.h"
#include "slump.h"

using namespace std;

int main() {
  double a;
  int nmon,nchain,cnt=0;
  bool coll;
  slump s;
  vector<particle> p(1);
  
  cout << "Number of chains attached to center ? "; cin >> nchain;
  cout << "Number of monomers per chain ? "; cin >> nmon;
  cout << "Monomer radius ? "; cin >> a;

  p[cnt++].radius=a;      // center

  for (int n=0; n<nchain; n++) {
    for (int i=0; i<nmon; i++) {
      coll=true;
      p.resize( cnt+1 );
      if (i==0)
        p[cnt]=p[0]; // use central atom for first monomer
      else
        p[cnt]=p[cnt-1];
      while (coll==true) {
        p[cnt].x+=2*a*s.random_half();
        p[cnt].y+=2*a*s.random_half();
        p[cnt].z+=2*a*s.random_half();
        coll=false;
        for (int j=0; j<cnt; j++)
          if (p[j].overlap(p[cnt])==true) {
            coll=true;
            break;
          } 
      }
      cnt++;
    }
  }

  cout << "STAR WITH " << nchain << " CHAIN(S) EACH CONTAINING " << nmon << " MONOMERS.\n";
  printf("%5d\n",p.size());
  for (int i=0; i<p.size(); i++)
    printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
        1,"STR","MON",i+1,p[i].x/10,p[i].y/10,p[i].z/10);
}
