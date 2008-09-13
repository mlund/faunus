/*!\page test_twobody Twobody
 * This program will simulate two protein molecules
 * in a dielectric solvent with explicit mobile ions.
 * The container is SPHERICAL w. hard walls.
 *
 * \author Mikael Lund
 * \date Prague 2007
 * \todo Maybe use a cylindrical cell?
 * \include twobody.C
 */

#include "faunus/faunus.h"
//#include "faunus/analysis.h"
//#include "faunus/histogram.h"
//#include "faunus/mcloop.h"
namespace Faunus {
  typedef pot_coulomb T_pairpot; // Specific pair interaction function
}
#include "faunus/markovmove.h"
using namespace Faunus;
using namespace std;

int main() {
  slump slump;                          // A random number generator
  inputfile in("twobody.conf");         // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  cell cell(in);                        // We want a spherical, hard cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg(in);                    // Setup pair potential (default values)
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  ioxyz xyz(cell);                      // xyz output for VMD etc.
  distributions dst;                    // Distance dep. averages
  iopqr pqr(cell);                      // PQR output (pos, charge, radius)
  FAUrdf saltrdf(particle::NA,particle::CL, .5, cell.r);
  twostatebinding bind(20.);            // Two state binding model

  vector<macromolecule> g;              // Group for proteins
  dualmove dm(nvt, cell, pot);          //   Class for 1D macromolecular translation
  dm.load( in, g, 40.);                 //   Load proteins and separate them 
  dm.rmax=60;
  salt salt;                            // Group for mobile ions
  salt.add(cell, in);                   //   Add salt particles
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  dm.dp=6;                              // Set displacement parameters
  sm.dp=90;                             // Set displacement parameters
  chargereg tit(nvt,cell,pot,salt,4.0); // Prepare titration.

  ioaam aam(cell);                      // Protein input file format is AAM
  if (aam.load(cell,"confout.aam")) {
    g[0].masscenter(cell);              // Load old config (if present)
    g[1].masscenter(cell);              // ...and recalc mass centers
  }
  
  systemenergy sys(pot.energy(cell.p)); // System energy analysis
  cout << cell.info() << pot.info()     // Print information to screen
       << tit.info();

  #ifdef GROMACS
  ioxtc xtc(cell, cell.r);              // Gromacs xtc output (if installed)
  #endif

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
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
          if (tit.du!=0) {                      // Average charges and dipoles
            dst.add("Q1",dm.r, g[0].charge(cell.p));
            dst.add("Q2",dm.r, g[1].charge(cell.p));
            dst.add("MU1",dm.r,g[0].dipole(cell.p));
            dst.add("MU2",dm.r,g[1].dipole(cell.p));
          }
          break;
      }
      if (slump.random_one()>.5 && macro>1) {
        //saltrdf.update(cell);                   // Analyse salt g(r)
        bind.update(cell, cell.p[g[0].beg], g[1]);
      }

      #ifdef GROMACS
      if (slump.random_one()>.95 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif
    } // End of inner loop

    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector

    dm.gofr.write("rdfprot.dat");               // Write interprotein g(r)
    saltrdf.write("rdfsalt.dat");               // Write salt g(r)
    dst.write("distributions.dat");             // Write other distributions
    xyz.save("coord.xyz", cell.p);              // Write .xyz coordinate file
    aam.save("confout.aam", cell.p);            // Save config. for next run
    pqr.save("confout.pqr", cell.p);            // ...also save a PQR file
    cout << loop.timing(macro);                 // Show progress
  } // End of outer loop

  cout << salt.info(cell)                       // Final information...
       << sm.info() << mr.info() << dm.info()
       << sys.info() << g[0].info() << g[1].info() << cell.info()
       << tit.info() << bind.info( 1/cell.getvolume() );

  #ifdef GROMACS
  xtc.close();
  #endif
}

