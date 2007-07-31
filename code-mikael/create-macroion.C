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
#include <string>
#include <vector>
#include "point.h"
#include "slump.h"

using namespace std;

int main() {
  int cnt=0;
  bool coll;
  string atom,countername;
  slump s;

  double box=30;
  int macrocharge=10;
  int macroradius=10;
  int counterradius=2;

  vector<particle> p(2+2*macrocharge);
  p[cnt].radius=macroradius;
  p[cnt++].z=0;
  p[cnt].radius=macroradius;
  p[cnt++].z=2*macroradius;
  
  if (macrocharge>0)
    countername="CL";
  else
    countername="NA";
  
  while (cnt<p.size()) {
    coll=true;
    while (coll==true) {
      p[cnt].x+=box*s.random_one();
      p[cnt].y+=box*s.random_one();
      p[cnt].z+=box*s.random_one();
      coll=false;
      for (int j=0; j<cnt; j++)
        if (p[j].overlap(p[cnt])==true) {
          coll=true;
          break;
        } 
    }
    cnt++;
  }

  cout << "Two macroions with radius " << macroradius << " and " << 2*macrocharge << " counterions.\n";
  printf("%5d\n",p.size());
  for (int i=0; i<p.size(); i++) {
    if (i<2)
      atom="MAC";
    else
      atom=countername;
    printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
        1,atom.c_str(),atom.c_str(),i+1,p[i].x/10,p[i].y/10,p[i].z/10);
  }
}
