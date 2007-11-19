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
#include "inputfile.h"
typedef pot_minimage T_pairpot;         // Specific pair interaction function
#include "markovmove.C"

using namespace std;

int main() {
  inputfile in("twobody.conf");         // Read input file
  slump slump;                          // A random number generator
  box cell(90.);                        // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential (default values)
  cfg.box = cell.len;
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  countdown<int> clock(10);             // Estimate simulation time
  ioxyz xyz(cell);                      // xyz output for VMD etc.
  ioaam aam(cell);                      // Protein input file format is AAM
  rdf protrdf(0,0,.5,cell.len/2);       // Protein and salt radial distributions, g(r)
  rdf saltrdf(particle::NA,particle::SO4, .5, cell.len/2);

  vector<macromolecule> g(12);          // Vector of proteins
  for (short i=0; i<g.size(); i++)      // Insert proteins...
    g[i].add( cell,                     // ...at random, non-overlapping
        aam.load("mrh4a.aam"), true );  // ...positions.

  group salt;                           // Group for mobile ions
  salt.add( cell, particle::NA, 0+880); // Insert sodium ions
  salt.add( cell, particle::SO4,30+440);// Insert chloride ions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  translate mt(nvt, cell, pot);         // Class for macromolecular translation
  systemenergy sys(pot.energy(cell.p)); // System energy analysis
  cout << cell.info();                  // Some information

  #ifdef GROMACS
  ioxtc xtc(cell);                      // Gromacs xtc output (if installed)
  #endif

  for (int macro=1; macro<=10; macro++) {       // Markov chain 
    for (int micro=1; micro<=6e4; micro++) {

      sys+=sm.move(salt);                       // Displace salt particles
      for (int i=0; i<g.size(); i++) {          // Loop over proteins
        if (slump.random_one()<0.5)
          sys+=mr.move(g[i]);                   // Rotate...
        else {
          sys+=mt.move(g[i]);                   // ...or translate
          for (int j=0; j<g.size(); j++)        // Analyse protein g(r)
            if (j!=i)
              protrdf.update(cell,g[i].cm,g[j].cm);
        }
      }

      if (slump.random_one()>0.9) {
        saltrdf.update(cell);                   // Analyse salt g(r)
        sm.adjust_dp(60,70);                    // Tune displacement
        mr.adjust_dp(60,70);                    // ...parameters
        mt.adjust_dp(60,70);
        #ifdef GROMACS
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
        #endif
      }
    } // End of inner loop

    cout << "Macro step " << macro << " completed. ETA: " << clock.eta(macro);
    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector

    xyz.save("coord.xyz", cell.p);              // Write .xyz coordinate file
    protrdf.write("rdfprot.dat");               // Write g(r)'s
    saltrdf.write("rdfsalt.dat");               //   -//-

  } // End of outer loop

  cout << sys.info() << salt.info()             // Final information...
       << sm.info() << mr.info() << mt.info();

  #ifdef GROMACS
  xtc.close();
  #endif
}

