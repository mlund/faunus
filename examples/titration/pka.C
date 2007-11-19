#include <iostream>
#include "io.h"
#include "analysis.h"
#include "container.h"
#include "potentials.h"
#include "countdown.h"
typedef pot_coulomb T_pairpot;          // Specific pair interaction function
#include "markovmove.C"

using namespace std;

int main() {
  slump slump;
  cell cell(100.);                      // We want a spherical cell
  iopov povray(cell);                   // We want a POVRAY snapshot
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential (default values)
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  countdown<int> clock(10);             // Estimate simulation time
  macromolecule protein;                // Group for the protein
  ioaam aam(cell);                      // Protein input file format is AAM
  protein.add( cell, aam.load(
        "calbindin.aam" ) );            // Load protein from disk
  protein.move(cell, -protein.cm);      // ..and translate it to origo (0,0,0)
  protein.accept(cell);                 // ..accept translation

  group salt;                           // Group for mobile ions
  salt.add( cell, particle::NA, 43);    // Insert sodium ions
  salt.add( cell, particle::CL, 24 );   // Insert chloride ions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  chargereg tit(nvt,cell,pot,salt,7.6); // Prepare titration. pH 7.6
  systemenergy sys(pot.energy(cell.p)); // System energy analysis

  cout << cell.info() << tit.info();    // Some information

  for (int macro=1; macro<=10; macro++) {       // Markov chain
    for (int micro=1; micro<=1e3; micro++) {
      sm.move(salt);                            // Displace salt particles
      sys+=sm.du;
      if (tit.titrateall()==true) {             // Titrate groups
        protein.charge(cell.p);                 // Re-calc. protein charge
        protein.dipole(cell.p);                 // Re-calc. dipole moment
        sys+=tit.du;
      }
      sys+=sm.du;                               // Keep system energy updated
    }

    cout << "Macro step " << macro << " completed. ETA: " << clock.eta(macro);
    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();
  }
  cout << sys.info() << sm.info() << tit.info() // More information...
       << salt.info() << protein.info();
  povray.save("protein-example.pov", cell.p);   // Save POVRAY file
}

