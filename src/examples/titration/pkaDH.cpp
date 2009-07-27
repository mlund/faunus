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
  interaction<pot_debyehuckel> pot(in);// Specify pair potential
  macromolecule protein;               // Group for the protein
  ioaam aam;                           // Protein input file format is AAM
  iopqr pqr;                           // PQR coordinate output
  protein.add( con,
      aam.load(in.getstr("protein"))); // Load protein structure
  protein.move(con, -protein.cm);      // ..translate it to origo (0,0,0)
  protein.accept(con);                 // ..accept translation
  aam.load(con, "confout.aam");        // Load old config (if present)

#ifdef DHTEIXEIRA
  //Andre: instantiate you class here!
  //ATchargereg tit(...);
  DHchargereg tit(nvt,con,pot,in.getflt("pH", 7.),in.getflt("mu_proton")); // UNCOMMENT THIS!!
#else
  DHchargereg tit(nvt,con,pot,in.getflt("pH", 7.),in.getflt("mu_proton"));
#endif  

  systemenergy sys(pot.energy(con.p)); // System energy analysis
  cout << con.info() << tit.info()     // Some information
       << pot.info() << atom.info();

  while ( loop.macroCnt() ) {            // Markov chain 
    while ( loop.microCnt() ) {
      sys+=tit.titrateall();             // Titrate protein sites
      protein.charge(con.p);             // Re-calc. protein charge
      protein.dipole(con.p);             // Re-calc. dipole moment
    }                                    // END of micro loop
    sys.update(pot.energy(con.p));       // Update system energy
    aam.save("confout.aam", con.p);      // Save config. to disk
    //pqr.save("confout.pqr", con.p, tit); // Save PQR file to disk - cool in VMD!
    cout << loop.timing();               // Show progress
  }                                      // END of macro loop
  cout << sys.info() << tit.info()
       << protein.info() << loop.info(); // Print final results
}

