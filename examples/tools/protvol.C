/*! \file Calculate the volume of a protein
 *  
 *
 */
#include <iostream>
#include "faunus/io.h"
#include "faunus/analysis.h"
#include "faunus/container.h"
#include "faunus/potentials.h"
namespace Faunus {
  typedef pot_minimage T_pairpot;
}
#include "faunus/markovmove.h"
#include "faunus/analysis.h"
#include "faunus/hardsphere.h"

using namespace Faunus;
using namespace std;

int main() {
  inputfile in("protvol.conf");
  cell con(in.getflt("radius"));         // We want a spherical cell
 
  macromolecule protein;
  ioaam aam(con);
  protein.add(con, aam.load(in.getstr("protein")));
  protein.center(con);
  protein.accept(con);
  hardsphere hs;
  particle probe;
  probe.radius=0;
  int hit=0;
 
  for (int macro=0; macro<in.getint("macro"); macro++) {        // Markov chain
    for (int micro=0; micro<in.getint("micro"); micro++) {
      con.randompos(probe);
      if(hs.overlap(con.p, protein, probe)==true)
        hit++;
      
    }
  }
  cout << con.info()<<endl <<protein.info()<<endl
       << "#  Protein radius      = "<< protein.radius(con.p) <<endl
       << "#  Number of shots     = "<< in.getint("macro")*in.getint("micro") <<endl
       << "#  Number of hist      = "<< hit << endl
       << "#  Hit ratio           = "<< hit/double(in.getint("macro")*in.getint("micro")) <<endl
       << "#  _____________________________________" << endl
       << "#  Protein volume is " << con.getvolume()*hit/(in.getint("macro")*in.getint("micro"))<<" A^3"<<endl;
   
}

