/*
 * This program will simulate an arbitrary number of
 * protein molecules in a dielectric solvent with
 * explicit mobile ions.
 *
 * \author Mikael Lund
 * \date Prague 2007
 */

#include <iostream>
#include "io.h"
#include "analysis.h"
#include "container.h"
#include "potentials.h"
#include "countdown.h"
#include "histogram.h"
typedef pot_minimage T_pairpot;          // Specific pair interaction function
#include "markovmove.C"

using namespace std;

int main() {
  slump slump;                          // A random number generator
  box cell(100.);                       // We want a spherical cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential (default values)
  cfg.box = cell.len;
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  countdown<int> clock(10);             // Estimate simulation time
  rdf rdf;                              // Protein radial distribution, g(r)
  ioxyz xyz(cell);

  ioaam aam(cell);                      // Protein input file format is AAM
  vector<macromolecule> g(12);          // Vector of proteins
  for (short i=0; i<g.size(); i++) {    // Loop over all proteins...
    g[i].add( cell, aam.load(
          "mrh4a.aam" ) );              // Load structure from disk
    g[i].move(cell, -g[i].cm);          // ..and translate it to origo (0,0,0)
    g[i].accept(cell);                  // ..accept translation
    while (g[i].overlap(cell)==true) {  // Place proteins
      point a;                          // ..at random positions
      cell.randompos(a);                // ..within the cell
      g[i].move(cell, a);
      g[i].accept(cell);                // ..accept translation
    }
  }

  group salt;                           // Group for mobile ions
  salt.add( cell, particle::NA, 0 ); // Insert sodium ions
  salt.add( cell, particle::CL, 36 ); // Insert chloride ions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  translate mt(nvt, cell, pot);         // Class for macromolecular translation
  systemenergy sys(pot.energy(cell.p)); // System energy analysis
  cout << cell.info();                  // Some information

  #ifdef GROMACS
  ioxtc xtc(cell);                      // Gromacs xtc output (if installed)
  #endif

  for (int macro=1; macro<=10; macro++) {       // Markov chain -- outer loop
    for (int micro=1; micro<=1e2; micro++) {    //  - // -      -- inner loop

      sm.move(salt);                            // Displace salt particles
      sm.adjust_dp(40,50);
      sys+=sm.du;
      for (int i=0; i<g.size(); i++) {          // Loop over proteins
        if (slump.random_one()<0.5)
          sys+=mr.move(g[i]);                   // Rotate...
        else {
          sys+=mt.move(g[i]);                   // ...or translate
          for (int j=0; j<g.size(); j++)
            if (j!=i)
              rdf.update(cell,g[i].cm,g[j].cm); // Update mass center g(r)
        }
      }

      #ifdef GROMACS
      if (slump.random_one()>0.8)
        xtc.save("ignored-name.xtc", cell.p);
      #endif
    }
    cout << "Macro step " << macro << " completed. ETA: " << clock.eta(macro);
    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector
  }

  xyz.save("coord.xyz", cell.p);
  rdf.write("rdf.dat");

  cout << sys.info() << salt.info()               // Final information...
       << sm.info() << mr.info() << mt.info();

  #ifdef GROMACS
  xtc.close();
  #endif
}

