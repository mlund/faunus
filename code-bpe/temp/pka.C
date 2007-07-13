/*! \file Example how to simulate a NaCl solution
 *  \example protein-example.C
 */
#include <iostream>
#include "classes/io.h"
#include "classes/analysis.h"
#include "classes/container.h"
#include "classes/potentials.h"
#include "classes/countdown.h"
typedef pot_coulomb T_pairpot;
#include "classes/markovmove.C"

using namespace std;

int main() {

  /****************************
  Declarations
  ****************************/

  double cell_r;
  double peeage;
  double Zobs;

  int micro;
  int macro;
  int nion1;
  int nion2;
  

  cell cell(92.);                       // We want a spherical cell
  iopov povray(cell);                   // We want a POVRAY snapshot
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  countdown<int> clock(10);             // Estimate simulation time
  widom<T_pairpot> wid(cell, pot, particle::NA
        , particle::CL);
  slump slp;

  macromolecule protein;                // Group for the protein
  ioaam aam(cell);                      // Protein file format is AAM
  protein.add( cell, aam.load(
        "protein-example.aam" ));    // Load protein from disk
  protein.move(cell, -protein.cm);      // ..and move it to origo
  protein.accept(cell);                 // (accept move)

  group salt;                           // Group for mobile ions
  salt.add( cell, particle::NA, 11+19); // Insert sodium ions
  salt.add( cell, particle::CL, 11 );   // Insert chloride ions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  chargereg tit(nvt,cell,pot,salt,7);   // Prepare titration. pH 7
  systemenergy sys(pot.energy(cell.p)); // System energy analysis

  cout << cell.info() << tit.info();    // Some information

  for (int macro=1; macro<=10; macro++) {       // Markov chain
    for (int micro=1; micro<=1e3; micro++) {
      sm.move(salt);                            // Displace salt particles
      if (tit.titrateall()) {                   // Titrate groups
        protein.charge(cell.p);                 // Re-calc. protein charge
        protein.dipole(cell.p);                 // Re-calc. dipole moment
      }
      sys+=sm.du+tit.du;                        // Keep system energy updated
      if (slp.random_one()<0.01) wid.insert();
    }
    cout << "Macro step " << macro << " completed. ETA: " << clock.eta(macro);
    sys.update(pot.energy(cell.p));
  }
  cout << sys.info() << sm.info() << tit.info() // More information...
    << salt.info() << protein.info() <<wid.info();
  povray.save("protein-example.pov", cell.p);   // Save POVRAY file
}

