/*
 * This program will simulate two protein molecules
 * in a dielectric solvent with explicit mobile ions.
 * The container is SPHERICAL w. hard walls.
 *
 * \author Mikael Lund
 * \date Prague 2007
 * \todo Maybe use a cylindrical cell?
 */

#include "analysis.h"
#include "profile.h"
#include "mcloop.h"
#include "pot_netz.h"
typedef pot_netz T_pairpot;             // Specific pair interaction function
#include "markovmove.h"

using namespace std;

int main() {
  slump slump;                          // A random number generator
  inputfile in("twobody.conf");         // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  cell cell(in);                        // We want a spherical, hard cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg(in);                    // Setup pair potential (default values)
  int_hydrophobic<T_pairpot> pot(cfg);  // Interactions incl. hydrophobic surfaces
  ioxyz xyz(cell);                      // xyz output for VMD etc.
  distributions dst;                    // Distance dep. averages
  iopqr pqr(cell);                      // PQR output (pos, charge, radius)
  rdf saltrdf(particle::NA,particle::I, .5, cell.r);
  cylindric_profile cyl(16,particle::I,-40,40,.5);

  vector<macromolecule> g;              // Group for proteins
  dualmove dm(nvt, cell, pot);          //   Class for 1D macromolecular translation
  dm.setup(in);                         //   Read input params. (optional)
  dm.load( in, g );                     //   Load proteins and separate them 
  salt salt;                            // Group for mobile ions
  salt.add(cell, in);                   //   Add salt particles
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  sm.dp=90;                             // Set displacement parameters
  ioaam aam(cell);                      // Protein input file format is AAM
  if (aam.load(cell,"confout.aam")) {
    g[0].masscenter(cell);              // Load old config (if present)
    g[1].masscenter(cell);              // ...and recalc mass centers
  }
  pot.search(cell.p);                   // Find hydrophobic particles
  pot.end_of_protein_one=g[0].end;      // Hydrophobic interactions w. BOTH proteins
  systemenergy sys(pot.energy(cell.p)); // System energy analysis
  cout << sys.info();

  cout << "# ------ INITIAL INFORMATION ------" << endl
       << cell.info() << pot.info()     // Print information to screen
       << loop.info() << endl
       << "# ------ RUNTIME INFORMATION ------" << endl;

  #ifdef GROMACS
  ioxtc xtc(cell, cell.r);              // Gromacs xtc output (if installed)
  #endif
  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,n;
      switch (rand() % 3) {                     // Pick a random MC move
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
      }
      if (slump.random_one()>.2 && macro>1) {
        saltrdf.update(cell);                   // Analyse salt g(r)
        cyl.update(cell.p);
      }

      #ifdef GROMACS
      if (slump.random_one()>.95 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif
    } // End of inner loop

    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector
    dm.gofr.write("rdfprot.dat");               // Write interprotein g(r)
    //saltrdf.write("rdfsalt.dat");               // Write salt g(r)
    //dst.write("distributions.dat");             // Write other distributions
    //xyz.save("coord.xyz", cell.p);              // Write .xyz coordinate file
    aam.save("confout.aam", cell.p);            // Save config. for next run
    pqr.save("confout.pqr", cell.p);            // ...also save a PQR file
    cout << loop.timing(macro);                 // Show progres;s
  } // End of outer loop

  cout << "# ------ FINAL INFORMATION ------" << endl
       << salt.info(cell)                       // Final information...
       << sm.info() << mr.info() << dm.info()
       << sys.info() << g[0].info() << g[1].info() << cell.info()
       << loop.info();
  cyl.show();

  #ifdef GROMACS
  xtc.close();
  #endif
}

