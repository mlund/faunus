/*! \page test_rosenbluth Rosenbluth
 *
 * Protein proton titration using grand canonical
 * salt moves.
 *
 * \author Bjorn Persson and Mikael Lund
 * \date Lund, 2010
 */
#include "faunus/faunus.h"
#include "faunus/potentials/pot_hscoulomb.h"
#include "faunus/potentials/pot_hsminimage.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main(int argc, char* argv[]) {
  cout << faunus_splash();              // Show Faunus information
  string config = "gcpka.conf";         // Default input (parameter) file
  if (argc==2) config = argv[1];        // ..also try to get it from the command line
  inputfile in(config);                 // Read input file
  checkValue test(in);                  // Enable unittesting
  mcloop loop(in);                      // Set Markov chain loop lengths
  grandcanonical nmt;                   // Use the canonical ensemble
  canonical nvt;
#ifdef MI
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_hsminimage> pot(in);  // ...and a plain Coulomb/HS potential
#else
  cell cell(in);                        // We want a spherical simulation container
  interaction<pot_hscoulomb> pot(in);   // ...and a Coulomb/HS pot. w. minimum image
#endif
  io io;                                // File i/o
  iopqr pqr;                            // ...pqr files
  ioaam aam;                            // ...aam files
  macromolecule protein;                // Group for the protein
  protein.add( cell,
      aam.load(in.getstr("protein")));  // Load protein structure
  protein.move(cell, -protein.cm);      // ..translate it to origo (0,0,0)
  protein.accept(cell);                 // ..accept translation
  saltmove sm(nmt,cell,pot,in);         // Class for salt movements
  salt salt;                            // Define some groups for mobile ions
  salt.add(cell,in);                    // Insert some ions

  saltbath sb(nmt,cell,pot,in,salt);    // Class for salt movements
  GCchargereg tit(nmt,cell,pot,in);

  widomSW wid2(10);                     // Class for single particle insertion w. charge scaling
  wid2.add( atom("NA") );
  wid2.add( atom("CA") );
  wid2.add( atom("LA") );
  wid2.add( atom("CL") );

  if(nmt.load(cell, "gcgroup.conf")==true)
    aam.load(cell,"confout.aam");        // Read initial config. from disk (if present)

  systemenergy sys(pot.energy(cell.p));  // Track system energy

  cout << cell.info() << atom.info()
       << pot.info() << salt.info(cell)
       << in.info() << tit.info()
       << endl;                         // Print initial information

  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
      wid2.insert(cell,pot);            // - // -
      int m=salt.size();
      for (int i=0; i<salt.size(); i++)
        sys+=sb.move();                 // Grand Canonical salt move
      sys+=tit.titrateall();
      sys.update(pot.energy(cell.p));   // Update system energy
      protein.charge(cell.p);           // Re-calc. protein charge
      protein.dipole(cell.p);           // Re-calc. dipole moment
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p));     // Update system energy
    aam.save("confout.aam",cell.p);     // Save particle configuration to disk
    cout << loop.timing();              // Show progres
  }                                     // END of macro loop and simulation

  io.writefile("gcgroup.conf", nmt.print());
  pqr.save("confout.pqr", cell.p);

  // Unit testing
  sm.check(test);
  sb.check(test);
  sys.check(test);
  test.check("ProteinCharge", protein.Q.avg() );
  test.check("ProteinDipole", protein.dip.avg() );

  cout << cell.info() << sys.info() << sm.info() << sb.info() << tit.info()
    << loop.info() << wid2.info() << protein.info() << tit.info()
    << test.report();

  return test.returnCode();
}

