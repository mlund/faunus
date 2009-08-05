#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  cout << faunus_splash();             // Faunus spam
#ifdef DHTEIXEIRA
  string config = "pkaAT.conf";          // Default input (parameter) file
#else
  string config = "pka.conf";          // Default input (parameter) file
#endif
  if (argc==2) config = argv[1];       // ..also try to get it from the command line
  inputfile in(config);                // Read input file
  mcloop loop(in);                     // Set Markov chain loop lengths
  cell con(in);                        // Use a spherical simulation container
  canonical nvt;                       // Use the canonical ensemble
  interaction<pot_debyehuckel> pot(in);// Specify pair potential
  vector <macromolecule> protein(1);      // Group for the protein
  ioaam aam;                           // Protein input file format is AAM
  iopqr pqr;                           // PQR coordinate output
  protein[0].add( con,
      aam.load(in.getstr("protein"))); // Load protein structure
  protein[0].move(con, -protein[0].cm);      // ..translate it to origo (0,0,0)
  protein[0].accept(con);                 // ..accept translation
  aam.load(con, "confout.aam");        // Load old config (if present)

#ifdef DHTEIXEIRA
  ATchargereg tit(nvt,con,pot,in.getflt("pH", 7.),in,pot.pair);
  protein[0].conc = in.getflt("ProteinConc", 0.0001);
  protein[0].cm.radius = protein[0].vradius(con.p);
#else
  DHchargereg tit(nvt,con,pot,in.getflt("pH", 7.),in.getflt("mu_proton"));
#endif

  systemenergy sys(pot.energy(con.p)); // System energy analysis
  cout << con.info() << tit.info()     // Some information
       << pot.info() << atom.info();

  while ( loop.macroCnt() ) {            // Markov chain
    while ( loop.microCnt() ) {
      #ifdef DHTEIXEIRA
        sys+=tit.titrateall( protein );
      #else
        sys+=tit.titrateall();             // Titrate protein sites
      #endif
      protein[0].charge(con.p);             // Re-calc. protein charge
      protein[0].dipole(con.p);             // Re-calc. dipole moment
    }                                    // END of micro loop
    sys.update(pot.energy(con.p));       // Update system energy
    aam.save("confout.aam", con.p);      // Save config. to disk
    //pqr.save("confout.pqr", con.p, tit); // Save PQR file to disk - cool in VMD!
    cout << loop.timing();               // Show progress
  }                                      // END of macro loop
  cout << sys.info() << tit.info()
       << protein[0].info(con) << loop.info(); // Print final results
}

