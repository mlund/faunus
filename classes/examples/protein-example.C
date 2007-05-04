/*! \file Example how to simulate a NaCl solution
 *  \example protein-example.C
 */
#include <iostream>
#include "../io.h"
#include "../analysis.h"
#include "../container.h"
#include "../potentials.h"
typedef pot_coulomb T_pairpot;
#include "../markovmove.C"

using namespace std;

int main() {
  cell cell(100.);                      // We want a spherical cell
  iopov povray(cell);                   // We want a POVRAY snapshot
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential
  interaction<T_pairpot> pot(cfg);      // Functions for interactions

  macromolecule protein;                // Group for the protein
  ioaam aam(cell);                      // Protein file format is AAM
  protein=cell.append(
      aam.load("protein-example.aam")); // Load protein from disk
  cell.move(protein, -protein.cm);      // ..and move it to origo
  cell.accept(protein);                 // (accept move)

  group salt;                           // Group for mobile ions
  salt+=cell.insert( particle::NA, 34 );// Insert sodium ions
  salt+=cell.insert( particle::CL, 15 );// Insert chloride ions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  chargereg tit(nvt,cell,pot,salt,7);   // Prepare titration. pH 7

  cout << cell.info() << salt << protein << tit.info();

  double u=0,u0=pot.energy(cell.p);
  for (int macro=0; macro<10; macro++) {        // Markov chain
    for (int micro=0; micro<1e3; micro++) {
      sm.move(salt);                            // Displace salt particles
      tit.titrateall();                         // Titrate groups
    }
    cout << pot.energy(cell.p) - (u0 + sm.utot + tit.utot) << endl;
  }
  cout << sm.info() << tit.info();
  povray.save("protein-example.pov", cell.p);    // Save POVRAY file
}

