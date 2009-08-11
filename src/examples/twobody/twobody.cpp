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

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  cout << faunus_splash();              // Faunus info
  string config = "twobody.conf";       // Default input (parameter) file
  if (argc==2) config = argv[1];        // ..also try to get it from the command line
  inputfile in(config);                 // Read input file
  checkValue test(in);                  // Test output
  mcloop loop(in);                      // Set Markov chain loop lengths
  cell cell(in);                        // We want a spherical, hard cell
  canonical nvt;                        // Use the canonical ensemble
  interaction<pot_coulomb> pot(in);     // Functions for interactions
  distributions dst;                    // Distance dep. averages
  FAUrdf saltrdf(atom["NA"].id,atom["CL"].id, .5, cell.r);
  twostatebinding bind(20.);            // Two state binding model

  vector<macromolecule> g;              // Group for proteins

  dualmove dm(nvt, cell, pot);          // Class for 1D macromolecular translation
  dm.load( in, g, 80.);                 // Load proteins and separate them 
  dm.rmax=60;
  salt salt;                            // Group for mobile ions
  salt.add(cell, in);                   //   Add salt particles
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  dm.dp=6;                              // Set displacement parameters
  sm.dp=90;                             // Set displacement parameters

  iopqr pqr;                            // PQR output (pos, charge, radius)
  ioxtc xtc(1000.);                     // Gromacs xtc output
  ioaam aam;                            // Protein input file format is AAM
  //if (aam.load(cell,"confout.aam")) {
  aam.load(cell,"confout.aam");
    g[0].masscenter(cell);              // Load old config (if present)
    g[1].masscenter(cell);              // ...and recalc mass centers
  //}

  chargereg tit(nvt,cell,pot,salt,4.0); // Prepare titration.
  
  systemenergy sys(pot.energy(cell.p)); // System energy analysis

  cout << in.info() << cell.info()
       << tit.info() << pot.info();     // Print information to screen

  while ( loop.macroCnt() ) {                   //Markov chain 
    while (loop.microCnt() ) {
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
      if (slp.random_one()>.5) {
        //saltrdf.update(cell);                   // Analyse salt g(r)
        //bind.update(cell, cell.p[g[0].beg], g[1]);
      }

      if (slp.random_one()>.95 && loop.macro>1)
        xtc.save("coord.xtc", cell.p);          // Save trajectory
    } // End of inner loop

    sys.update(pot.energy(cell.p));             // System energy analysis
    cell.check_vector();                        // Check sanity of particle vector

    dm.gofr.write("rdfprot.dat");               // Write interprotein g(r)
    saltrdf.write("rdfsalt.dat");               // Write salt g(r)
    dst.write("distributions.dat");             // Write other distributions
    aam.save("confout.aam", cell.p);            // Save config. for next run
    pqr.save("confout.pqr", cell.p);            // ...also save a PQR file (no salt)
    cout << loop.timing();                      // Show progress

  } // End of outer loop

  xtc.close();                                  // Close xtc file for writing

  cout << salt.info(cell)                       // Final information...
       << sm.info() << mr.info() << dm.info()
       << sys.info() << g[0].info() << g[1].info() << cell.info()
       << tit.info() << bind.info( 1/cell.getvolume() );

  // Output tests
  sys.check(test);
  test.check("Protein1_charge", g[0].Q.avg());
  test.check("Protein2_charge", g[1].Q.avg());

  cout << test.report();

  return test.returnCode();
}

