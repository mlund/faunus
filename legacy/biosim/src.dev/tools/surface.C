/*
 * Prope the surface of a protein and output the surface
 * coordinates either in cartesian- or spherical coordinates.
 *
 * Plot for example using GnuPlot:
 *   set polar
 *   plot "data" index 1 u 2:1
 *
 * M. Lund, November 2005
 *
 */

#include <iostream>
#include <string>
#include <vector>
#include "point.h"
#include "interact.h"
#include "group.h"
#include "space.h"
#include "slump.h"
#include "histogram.h"
#include "cell.h"

int main(int argc, char* argv[] ) {

  vector<point> surf(0);
  point cm,origo;     //here we manually specify CM
  cm.x= 34.2938;      //(needs to be center-of-mass
  cm.y= 1.35419;       //  in this case...)
  cm.z= 10.2892;

  string protein = "../struct/ribonuclease-allatom.aam";
  //string protein = "../struct/ribonuclease.aam";
  collision col;
  space s,chg;
  group g,g2;
  peptide pep;

  g = s.append( pep.loadAAModel(protein)  );
  s.addmasscenter(g);
  s.p[g.cm]=cm;
  s.move( g, -cm, space::AUTOACCEPT);
  g.radius=s.radius( g, s.p[g.cm] );
  double cell_r=g.radius+3.;
  cell cell(cell_r);

  particle p;
  spherical pol;
  p.radius=0.3;
 
  double pi=acos(-1.);
  for (pol.theta=0.05; pol.theta<pi; pol.theta+=0.05) {
    cout << "#" << pol.theta << endl;
    for (pol.phi=0; pol.phi<2*pi; pol.phi+=0.05)
      for (pol.r=cell_r; pol.r>0; pol.r-=0.1) {
        p=pol.cartesian();
        if (col.overlap(s.p, p)==true) {
          cout << pol.phi << " " <<pol.r << endl;
          surf.push_back(p);
          break;
        };
      };
    cout << endl << endl;
  };
  return 0;

  /*
  spherical org;
  point b;
  g2=chg.append(pep.loadAAModel("../struct/lysozyme-scm.aam")); 
  chg.move( g2, -cm, space::AUTOACCEPT);
  for (int i=0; i<chg.p.size(); i++) {
    pol=chg.p[i];
    org=pol;
    for (pol.r=cell_r; pol.r>0; pol.r-=0.1) {
      p=pol.cartesian();
      if (col.overlap(s.p, p)==true) {
	cout << org.r << " " << pol.r << endl; 
	break;
      };
    };
  };
  return 0;
  */

  double d2,d3;
  s.p.clear();
  g = s.append( pep.loadAAModel("../struct/ribonuclease-scm.aam")); 
  s.move( g, -cm, space::AUTOACCEPT);        
  for (int j=0; j<s.p.size(); j++) {
    double d=10000;
    for (int i=0; i<surf.size(); i++)
      if (s.p[j].dist(surf[i])<d) {
        d=s.p[j].dist(surf[i]);
        d2=surf[i].len();
	d3=s.p[j].len();
      };
    cout << d3 << " " << d2 << endl; 
  };
        
  
};
