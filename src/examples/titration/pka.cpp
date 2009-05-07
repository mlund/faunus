/*!\page test_pka Titration
 * Proton titration of a single macromolecule in an explicit salt
 * envirinment. Protons are moved in and out of the protein to
 * sample the average protonation state, pka values and dipole
 * moment as well as charge fluctuations.
 *
 * \author Mikael Lund and Bjorn Persson
 * \include pka.cpp
 */
#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  cout << faunus_splash();             // Faunus spam
  string config = "pka.conf";          // Default input (parameter) file
  if (argc==2) config = argv[1];       // ..also try to get it from the command line
  inputfile in(config);                // Read input file
  mcloop loop(in);                     // Set Markov chain loop lengths
  cell con(in);                        // Use a spherical simulation container
  canonical nvt;                       // Use the canonical ensemble
  interaction<pot_hscoulomb> pot(in);  // Specify pair potential
  macromolecule protein;               // Group for the protein
  ioaam aam(con.atom);                 // Protein input file format is AAM
  iopqr pqr(con.atom);                 // PQR coordinate output
  protein.add( con,
      aam.load(in.getstr("protein"))); // Load protein structure
  protein.move(con, -protein.cm);      // ..translate it to origo (0,0,0)
  protein.accept(con);                 // ..accept translation
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions
  saltmove sm(nvt, con, pot);          // Class for salt movements
  aam.load(con, "confout.aam");        // Load old config (if present)
  widom wid(10);
  wid.add(con.atom("NA"));
  wid.add(con.atom("CL"));
  wid.runfraction=0.05;

#ifdef GCPKA // "Grand Canonical" titration
  HAchargereg tit(nvt,con,pot,salt,in.getflt("pH", 7.),in.getflt("catpot"));
#else        // "Normal" titration
  chargereg tit(nvt,con,pot,salt,in.getflt("pH",7.));
#endif

  systemenergy sys(pot.energy(con.p)); // System energy analysis
  cout << con.info() << tit.info()     // Some information
       << pot.info() << con.atom.info();

  while ( loop.macroCnt() ) {            // Markov chain 
    while ( loop.microCnt() ) {
      switch (rand() % 2) {              // Randomly chose move
        case 0:
          sys+=sm.move(salt);            // Displace salt particles
          break;
        case 1:
          sys+=tit.titrateall();         // Titrate protein sites
          protein.charge(con.p);         // Re-calc. protein charge
          protein.dipole(con.p);         // Re-calc. dipole moment
          break;
      }
      wid.insert(con,pot);
    }                                    // END of micro loop
    sys.update(pot.energy(con.p));       // Update system energy
    aam.save("confout.aam", con.p);      // Save config. to disk
    pqr.save("confout.pqr", con.p, tit); // Save PQR file to disk - cool in VMD!
    cout << loop.timing();               // Show progress
  }                                      // END of macro loop
  cout << sys.info() << sm.info() << wid.info()
       << tit.info() << salt.info(con)
       << protein.info() << loop.info(); // Print final results
}

