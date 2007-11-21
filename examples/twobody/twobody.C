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
typedef pot_coulomb T_pairpot;         // Specific pair interaction function
#include "markovmove.h"

using namespace std;

int main() {
  inputfile in("twobody.conf");         // Read input file
  slump slump;                          // A random number generator
  cell cell(90.);                        // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential (default values)
  //cfg.box = cell.len;                   // Pass box len to pair potential
  cfg.eps = 1.0;
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  countdown<int> clock(10);             // Estimate simulation time
  ioxyz xyz(cell);                      // xyz output for VMD etc.
  ioaam aam(cell);                      // Protein input file format is AAM
  rdf protrdf(0,0,.5, 45.);       // Protein and salt radial distributions, g(r)
  rdf saltrdf(particle::NA,particle::SO4, .5, 45.);

  vector<macromolecule> g(0);          // Vector of proteins
  for (short i=0; i<g.size(); i++)      // Insert proteins...
    g[i].add( cell,
        aam.load("mrh4a.aam"), true );

  group salt;                           // Group for mobile ions
  salt.add( cell, particle::NA, 0+80); // Insert sodium ions
  salt.add( cell, particle::SO4,0+40);// Insert chloride ions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  translate mt(nvt, cell, pot);         // Class for macromolecular translation
  systemenergy sys(pot.energy(cell.p)); // System energy analysis
  cout << cell.info() << pot.info();    // Some information

  #ifdef GROMACS
    ioxtc xtc(cell, 45.);                      // Gromacs xtc output (if installed)
  #endif

  for (int macro=1; macro<=10; macro++) {       // Markov chain 
    for (int micro=1; micro<=3e2; micro++) {

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
              if (j!=i)
                protrdf.update(cell,g[i].cm,g[j].cm);
          }
          break;
      }

      if (slump.random_one()>0.8)
        saltrdf.update(cell);                   // Analyse salt g(r)

      #ifdef GROMACS
      if (slump.random_one()>0.3)
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

