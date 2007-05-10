/*! \file Example how to simulate a NaCl solution
 *  \example protein-example.C
 */
#include <iostream>
#include "../io.h"
#include "../analysis.h"
#include "../container.h"
#include "../potentials.h"
#include "../countdown.h"
typedef pot_test T_pairpot;
#include "../markovmove.C"

using namespace std;

int main() {
  clutch cell(100.,-15,15);             // We want a spherical cell
  iopov povray(cell);                   // We want a POVRAY snapshot
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  countdown<int> clock(10);             // Estimate simulation time

  macromolecule protein;                // Group for the protein
  ioaam aam(cell);                      // Protein file format is AAM
  protein=cell.append( aam.load(
   "examples/clutch-example-small.aam" ));    // Load protein from disk

  particle ghost;
  ghost.z=90.;
  ghost.charge=+1;
  ghost.radius=2.0;
  widompath<T_pairpot> widom(ghost,cell.p[3]);

  group salt;                           // Group for mobile ions
  salt+=cell.insert( particle::NA, 32+9);// Insert sodium ions
  salt+=cell.insert( particle::CL, 32);// Insert chloride ions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  chargereg tit(nvt,cell,pot,salt,7);   // Prepare titration. pH 7
  systemenergy sys(pot.energy(cell.p)); // System energy analysis

  cout << cell.info() << tit.info() << protein.info();

  for (int macro=1; macro<=10; macro++) {       // Markov chain
    for (int micro=1; micro<=3e4; micro++) {
      sm.move(salt);                            // Displace salt particles
      if (tit.titrateall()) {                   // Titrate groups
        protein.charge(cell.p);                 // Re-calc. protein charge
        protein.dipole(cell.p);                 // Re-calc. dipole moment
      }
      widom.update(cell,pot,povray);
      sys+=sm.du+tit.du;                        // Keep system energy updated
    }
    cout << "Macro step " << macro << " completed. ETA: " << clock.eta(macro);
    sys.update(pot.energy(cell.p));
    cout << widom.info();
  }
  cout << sys.info() << sm.info() << tit.info() << salt.info() << protein.info();
  widom.povpath(povray);
  povray.save("clutch-example.pov", cell.p);    // Save POVRAY file
}

