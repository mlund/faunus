/*
 * This program will simulate two protein molecules
 * in a dielectric solvent with explicit mobile ions.
 * The container is SPHERICAL w. hard walls.
 *
 * \author Mikael Lund
 * \date Prague 2007
 * \todo Maybe use a cylindrical cell?
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
  enum dstkeys { Q1=0,Q2,MU1,MU2, LAST };
  slump slump;                          // A random number generator
  inputfile in("twobody.conf");         // Read input file
  cell cell(in);                        // We want a spherical, hard cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg(in);                    // Setup pair potential (default values)
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  countdown<int> clock(10);             // Estimate simulation time
  ioxyz xyz(cell);                      // xyz output for VMD etc.
  distributions dst;
  rdf saltrdf(particle::NA,particle::CL, .5, cell.r);

  vector<macromolecule> g;              // Group for proteins
  dualmove dm(nvt, cell, pot);          //   Class for 1D macromolecular translation
  dm.load( in, g, 40.);                 //   Load proteins and separate them 
  salt salt;                            // Group for mobile ions
  salt.add(cell, in);                   //   Add salt particles
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation

  ioaam aam(cell);                      // Protein input file format is AAM
  if (aam.load(cell,"confout.aam")) {
    g[0].masscenter(cell.p);            // Load old config (if present)
    g[1].masscenter(cell.p);            // ...and recalc mass centers
  }
  chargereg tit(nvt,cell,pot,salt,4.0); // Prepare titration.
  systemenergy sys(pot.energy(cell.p)); // System energy analysis
  cout << cell.info() << pot.info();    // Print information to screen

  #ifdef GROMACS
  ioxtc xtc(cell, cell.r);              // Gromacs xtc output (if installed)
  #endif

  for (int macro=1; macro<=10; macro++) {       // Markov chain 
    for (int micro=1; micro<=2e3; micro++) {
      short i,n;
      switch (rand() % 4) {                     // Pick a random MC move
        case 0:                                 // Displace salt
          sys+=sm.move(salt);                   //   Do the move.
          break;
        case 1:                                 // Rotate proteins
          for (n=0; n<2; n++) {                 //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            sys+=mr.move(g[i]);                 //   Do the move.
          }
          break;
        case 2:                                 // Translate proteins
          sys+=dm.move(g[0], g[1]);             //   Do the move.
          break;
        case 3:                                 // Fluctuate charges
          sys+=tit.titrateall();                // Titrate sites on the protein
          if (tit.du!=0) {                      // Average charge and dipole
            dst.add("Q1",dm.r, g[0].charge(cell.p));
            dst.add("Q2",dm.r, g[1].charge(cell.p));
            dst.add("MU1",dm.r,g[0].dipole(cell.p));
            dst.add("MU2",dm.r,g[1].dipole(cell.p));
          }
          break;
      }
      if (slump.random_one()>.8 && macro>1)
        saltrdf.update(cell);                   // Analyse salt g(r)

      #ifdef GROMACS
      if (slump.random_one()>.95 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif

    } // End of inner loop

    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector

    xyz.save("coord.xyz", cell.p);              // Write .xyz coordinate file
    aam.save("confout.aam", cell.p);            // Save config. for next run
    dm.gofr.write("rdfprot.dat");
    saltrdf.write("rdfsalt.dat");
    dst.write("distributions.dat");
    cout << "Macro step " << macro
         << " completed. ETA: " << clock.eta(macro);
 
  } // End of outer loop

  cout << salt.info(cell)                       // Final information...
       << sm.info() << mr.info() << dm.info()
       << tit.info() << sys.info()
       << g[0].info() << g[1].info() << cell.info();

  #ifdef GROMACS
  xtc.close();
  #endif
}

