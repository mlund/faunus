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
#include "potentials.h"
#include "container.h"
#include "countdown.h"
#include "histogram.h"
#include "inputfile.h"
typedef pot_minimage T_pairpot;         // Specific pair interaction function
#include "markovmove.h"

using namespace std;

int main() {
  inputfile in("manybody.conf");        // Read input file
  slump slump;                          // A random number generator
  box cell(90.);                        // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg(in);                    // Setup pair potential (default values)
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  countdown<int> clock(10);             // Estimate simulation time
  ioxyz xyz(cell);                      // xyz output for VMD etc.
  rdf protrdf(0,0,.5,cell.len/2);       // Protein and salt radial distributions, g(r)
  rdf saltrdf(particle::NA,particle::SO4, .5, cell.len/2);

  vector<macromolecule> g;              // Protein groups
  ioaam aam(cell);                      // Protein input file format is AAM
  aam.load(cell, in, g);                // Load and insert proteins

  group salt;                           // Group for mobile ions
  salt.add(cell, in);                   // Add salt particles
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  translate mt(nvt, cell, pot);         // Class for macromolecular translation
  systemenergy sys(pot.energy(cell.p)); // System energy analysis
  cout << cell.info() << pot.info();    // Print information to screen

  #ifdef GROMACS
  ioxtc xtc(cell, cell.len/2.);         // Gromacs xtc output (if installed)
  #endif

  for (int macro=1; macro<=10; macro++) {       // Markov chain 
    for (int micro=1; micro<=1e3; micro++) {
      short i,j,n;
      switch (rand() % 3) {                     // Pick a random MC move
        case 0:                                 // Displace salt
          sys+=sm.move(salt);                   //   Do the move.
          break;
        case 1:                                 // Rotate proteins
          for (n=0; n<g.size(); n++) {          //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            sys+=mr.move(g[i]);                 //   Do the move.
          }
          break;
        case 2:                                 // Translate proteins
          for (n=0; n<g.size(); n++) {          //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            sys+=mt.move(g[i]);                 //   Do the move.
            for (j=0; j<g.size(); j++)          //   Analyse g(r)...
              if (j!=i && macro>1)
                protrdf.update(cell,g[i].cm,g[j].cm);
          }
          break;
      }
      if (macro==1 && micro<1e3) {
        mt.adjust_dp(30,40);                    // Adjust displacement
        mr.adjust_dp(40,50);                    // parameters. Use ONLY
        sm.adjust_dp(20,30);                    // during equillibration!
      }
      if (slump.random_one()>.8 && macro>1)
        saltrdf.update(cell);                   // Analyse salt g(r)

      #ifdef GROMACS
      if (slump.random_one()>.95 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif

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

