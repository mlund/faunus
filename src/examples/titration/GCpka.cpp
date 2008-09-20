/*\page test_GCpka Grand Canonical Titration
 * Program for perform semi-Grand Caninical titration
 * \author Bjorn Persson
 * \include GCpka.cpp
 */
#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  string config = "GCpka.conf";
  if (argc==2) config = argv[1];
  inputfile in(config);                // Read input file
  mcloop loop(in);                     // Set Markov chain loop lengths
  cell con(in);                        // Use a spherical container
  canonical nvt;                       // Use the canonical ensemble
  interaction<pot_coulomb> pot(in);    // Functions for interactions
  macromolecule protein;               // Group for the protein
  ioaam aam(con);                      // Protein input file format is AAM
  iopqr pqr(con);                      // PQR coordinate output
  protein.add( con,
      aam.load( in.getstr("protein")) );
  protein.move(con, -protein.cm);      // ..translate it to origo (0,0,0)
  protein.accept(con);                 // ..accept translation
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions
  saltmove sm(nvt, con, pot);          // Class for salt movements
  aam.load(con, "confout.aam");        // Load old config (if present)
  HAchargereg tit(nvt,con,pot,salt,in.getflt("pH"),in.getflt("catpot")); // Prepare titration. pH 7.6
  systemenergy sys(pot.energy(con.p)); // System energy analysis
  cout << con.info() << tit.info()     // Some information
    << pot.info();

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      switch (rand() % 2) {                 // Randomly chose move
        case 0:
          sys+=sm.move(salt);               // Displace salt particles
          break;
        case 1:
          sys+=tit.titrateall();            // Titrate protein sites
          protein.charge(con.p);            // Re-calc. protein charge
          protein.dipole(con.p);            // Re-calc. dipole moment
          break;
      }
    }                                       // END of micro loop
    sys.update(pot.energy(con.p));          // Update system energy
    aam.save("confout.aam", con.p);         // Save config. to disk
    pqr.save("confout.pqr", con.p, tit);    // Save PQR file to disk - cool in VMD!
    cout << loop.timing(macro);             // Show progress
  }                                         // END of macro loop
  cout << sys.info() << sm.info()           // Print final results
    << tit.info() << salt.info(con) << protein.info();
}

