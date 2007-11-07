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

  g = s.append( pep.loadAAModel(protein)  );

  s.addmasscenter(g);
  s.move( g, -s.center_of_mass(g), space::AUTOACCEPT);
  g.radius=s.radius( g, s.p[g.cm] );
  double cell_r=g.radius+10.;
  cell cell(cell_r);

  vector<histogram::data> h(1);
  histogram hist( 0.5, cell_r);
  hist.init(h[0]);
  h[0].name = "P(r)";
  h[0].type = histogram::AVERAGE;

  particle p, origo;
  p.radius=2.0;
  double r;
  for (int i=0; i<200000; i++) {
    cell.randomPos(p);
    r=p.dist(origo);
    if (col.overlap(s.p, p)==true)
      hist.add(h[0], r, 1.0);
    else
      hist.add(h[0], r, 0.0);
  };

  hist.show(h, 0, cell_r);
};
