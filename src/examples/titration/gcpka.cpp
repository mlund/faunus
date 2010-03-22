/*! \page test_rosenbluth Rosenbluth
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method.
 *
 * \author Christophe Labbez and Mikael Lund
 * \date Lund, 2009
 * \include rb.cpp
 */
#include "faunus/faunus.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  cout << faunus_splash();              // Show Faunus information
  inputfile in("gcsalt.conf");          // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  grandcanonical nmt;                   // Use the canonical ensemble
  canonical nvt;
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_hsminimage> pot(in);  // ...and a Coulomb/HS pot. w. minimum image
  io io;
  ioaam aam;                            // File I/O class
  macromolecule protein;               // Group for the protein
  protein.add( cell,
      aam.load(in.getstr("protein"))); // Load protein structure
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
    aam.load(cell,"rb.aam");           // Read initial config. from disk (if present)

  systemenergy sys(pot.energy(cell.p)); // - pot.internal(cell.p, protein)); // Track system energy

  cout << cell.info() << atom.info()
       << pot.info() << salt.info(cell)
       << in.info() << sb.info() <<tit.info();       // Print initial information


  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
      wid2.insert(cell,pot);            // - // -
      int m=salt.size();
      for (int i=0; i<salt.size(); i++)
        sys+=sb.move();                 // Grand Canonical salt move
      sys+=tit.titrateall();
        sys.update(pot.energy(cell.p)); // Update system energy
        protein.charge(cell.p);         // Re-calc. protein charge
        protein.dipole(cell.p);         // Re-calc. dipole moment
          
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p)); //-pot.internal(cell.p, protein));     // Update system energy
    aam.save("rb.aam",cell.p);          // Save particle configuration to disk
    cout << loop.timing()<<sys.info();  // Show progres
  }                                     // END of macro loop and simulation

  io.writefile("gcgroup.conf", nmt.print());

  cout << cell.info() << sys.info() << sm.info() << sb.info()<<tit.info()
       << loop.info()<<wid2.info()<<protein.info()<<tit.info();
  iopqr pqr;
  pqr.save("saltbath.pqr", cell.p);
}

