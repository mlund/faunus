/* This program will calculate the 'born moment'
 * of a molecule, defined as
 *
 * $\sum q_i^2*r_i/a_i*(1-1/eps)
 */

#include <iostream>
#include "point.h"
#include "group.h"
#include "io.h"
#include "container.h"

using namespace std;

int main() {
  cell cell(100.);
  ioaam aam(cell);
  macromolecule mol;
  mol.add(cell, aam.load("test.aam"));
  mol.move(cell, -mol.cm);
  mol.accept(cell);

  point b;
  double f;
  for (int i=mol.beg; i<=mol.end; i++) {
    f=pow(cell.p[i].charge,2)/cell.p[i].radius;
    b.x+=f*cell.p[i].x;
    b.y+=f*cell.p[i].y;
    b.z+=f*cell.p[i].z;
  }
  cout << b << " " << b.len() << endl;
}
