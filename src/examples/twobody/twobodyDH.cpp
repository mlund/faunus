/*!\page test_twobody Twobody
 * This program will simulate two protein molecules
 * in a dielectric solvent with explicit mobile ions.
 * The container is SPHERICAL w. hard walls.
 *
 * \author Mikael Lund
 * \date Prague 2007
 * \todo Maybe use a cylindrical cell?
 * \include twobody.cpp
 */


#include "faunus/faunus.h"

#ifdef LIBXTC
#define GROMACS
#endif

using namespace Faunus;
using namespace std;

int main() {
  cout << faunus_splash();              // Faunus info
  slump slump;                          // A random number generator
  inputfile in("twobodyDH.conf");       // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  cell cell(in);                        // We want a spherical, hard cell
  canonical nvt;                        // Use the canonical ensemble
  interaction<pot_debyehuckel> pot(in);  // Functions for interactions
  ioxyz xyz;                            // xyz output for VMD etc.
  distributions dst;                    // Distance dep. averages
//  angularcorr  angcorr;                 // Dipole cross-correlation
  iopqr pqr;                            // PQR output (pos, charge, radius)

  vector<macromolecule> g;              // Group for proteins
  dualmove dm(nvt, cell, pot);          //   Class for 1D macromolecular translation
  dm.load( in, g, 50.);                 //   Load proteins and separate them 
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  dm.setup(in);                         // Set displacement parameters
  mr.dp=in.getflt("rot_dp");            //  ---//---
  ioaam aam;                            // Protein input file format is AAM
  if (aam.load(cell,"confout.aam")) {
    g[0].masscenter(cell);              // Load old config (if present)
    g[1].masscenter(cell);              // ...and recalc mass centers
  }

  
  systemenergy sys(pot.energy(cell.p)-pot.internal(cell.p,g[0])-pot.internal(cell.p,g[1])); // System energy analysis
  cout << in.info() << cell.info()
       << pot.info();                   // Print information to screen

  #ifdef GROMACS
  ioxtc xtc(cell, cell.r);              // Gromacs xtc output (if installed)
  #endif

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,n;
      switch (rand() % 2) {                     // Pick a random MC move
        case 0:                                 // Rotate proteins
          for (n=0; n<2; n++) {                 //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            sys+=mr.move(g[i]);                 //   Do the move.
          }
          break;
        case 1:                                 // Translate proteins
          sys+=dm.move(g[0], g[1]);             //   Do the move.
          break;
      }
      if (slump.random_one() < .1 ) {
        double z, l0, l1;
        z=g[0].cm.dist(g[1].cm); 
        l0=g[0].mu.len();
        l1=g[1].mu.len();
        g[0].dipole(cell.p);
        g[1].dipole(cell.p);
        dst.add("Left  protein z-dipcomp", z , g[0].mu.z/l0); 
        dst.add("Right protein z-dipcomp", z , g[1].mu.z/l1); 
      }

      #ifdef GROMACS
      if (slump.random_one()>.95 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif
    } // End of inner loop

    sys.update(pot.energy(cell.p)-pot.internal(cell.p, g[0])-pot.internal(cell.p, g[1]));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector

    dm.gofr.write("rdfprot.dat");               // Write interprotein g(r)
    dst.write("distributions.dat");             // Write other distributions
    xyz.save("coord.xyz", cell.p);              // Write .xyz coordinate file
    aam.save("confout.aam", cell.p);            // Save config. for next run
    //pqr.save("confout.pqr", cell.p);            // ...also save a PQR file
    cout << loop.timing(macro);                 // Show progress
  } // End of outer loop

  cout <<                                       // Final information...
        mr.info() << dm.info()
       << sys.info() << g[0].info() << g[1].info() << cell.info();

  #ifdef GROMACS
  xtc.close();
  #endif
}

