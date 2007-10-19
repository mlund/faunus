/*
 * Calculate volume of a protein
 * M. Lund, October 2005
 */

#include <iostream>
#include <string>
#include "interact.h"
#include "group.h"
#include "space.h"
#include "slump.h"
#include "histogram.h"
#include "cell.h"

int main(int argc, char* argv[] ) {
  string protein = "../struct/barnase.aam";
  collision col;
  space s;
  group g;
  peptide pep;
  long int hit=0, miss=0;

  g = s.append( pep.loadAAModel(protein)  ); // Load protein

  s.addmasscenter(g);
  s.move( g, -s.center_of_mass(g), space::AUTOACCEPT);
  g.radius=s.radius( g, s.p[g.cm] );
  double cell_r=g.radius+10.;
  cell cell(cell_r);

  particle p;
  p.radius=0.1;
  double r,rmin=1e3;
  for (int i=0; i<4000000; i++) {
    cell.randomPos(p);
    if (col.overlap(s.p, p)==true)
      hit++;
    else {
      miss++;
      r=p.len();
      if (r<rmin) rmin=r;
    };
  };

  double vol=double(hit)/double(miss)*4/3*3.14*pow(cell_r,3);
   
  cout << "Volume   = " << vol << endl;
  cout << "R_min    = " << rmin << endl;
  cout << "R_sphere = " << pow(3*vol/2/3.14,0.333333) << endl;
};
