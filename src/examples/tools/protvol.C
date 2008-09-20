/*! \file Calculate the volume of a protein
 *  
 *
 */
#include "faunus/faunus.h"
#include "faunus/energy.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/moves/markovmove.h"

using namespace Faunus;
using namespace std;

int main() {
  inputfile in("protvol.conf");
  cell con(in.getflt("radius"));         // We want a spherical cell
 
  macromolecule protein;
  ioaam aam(con);
  protein.add(con, aam.load(in.getflt("protein")));
  protein.center(con);
  protein.accept(con);
  hardsphere hs;
  particle probe;
  probe.radius=0;
  int hit=0;
 
  for (int macro=0; macro<in.getflt("macro"); macro++) {        // Markov chain
    for (int micro=0; micro<in.getflt("micro"); micro++) {
      con.randompos(probe);
      if(hs.overlap(con.p, protein, probe)==true)
        hit++;
      
    }
  }
  cout << con.info()<<endl <<protein.info()<<endl
       << "#  Protein radius      = "<< protein.radius(con.p) <<endl
       << "#  Number of shots     = "<< in.getflt("macro")*in.getflt("micro") <<endl
       << "#  Number of hist      = "<< hit << endl
       << "#  Hit ratio           = "<< hit/double(in.getflt("macro")*in.getflt("micro")) <<endl
       << "#  _____________________________________" << endl
       << "#  Protein volume is " << con.getvolume()*hit/(in.getflt("macro")*in.getflt("micro"))<<" A^3"<<endl;
   
}

