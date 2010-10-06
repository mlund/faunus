/*! \page test_rosenbluth Rosenbluth
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method.
 *
 * \author Christophe Labbez and Mikael Lund
 * \date Lund, 2009
 */
#include "faunus/faunus.h"
#include "faunus/potentials/pot_hsminimage.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  cout << faunus_splash();              // Show Faunus information
  io io;
  inputfile in("gcsalt.conf");          // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  grandcanonical nmt;                   // Use the canonical ensemble
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_hsminimage> pot(in);  // ...and a Coulomb/HS pot. w. minimum image

  ioaam aam;                            // File I/O class

  salt salt;                            // Define some groups for mobile ions
  salt.add(cell,in);                    // Insert some ions

  saltmove sm(nmt,cell,pot,in);         // Class for salt movements
  saltbath sb(nmt,cell,pot,in,salt);    // Class for salt movements

  if(nmt.load(cell, "gcgroup.conf")==true)
    aam.load(cell,"rb.aam");           // Read initial config. from disk (if present)

  widom wid(10);                        // Class for multiple particle insertion
  wid.add(cell);                        // Detect all species in the cell

  systemenergy sys(pot.energy(cell.p)); // Track system energy

  cout << cell.info() << atom.info()
       << pot.info() << salt.info(cell)
       << in.info() << sb.info();       // Print initial information

  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
      int m=salt.size();
      for (int i=0; i<m; i++){
        sys+=sb.move();                 // Grand Canonical salt move
  /*      if(cell.charge()!=0) {
          std::cout<<cell.charge()<<endl;
          std::cout<<sb.thispair->i<<" "sb.thispiar->j<<endl;
        } */
      }
      wid.insert(cell,pot);             // Widom particle insertion analysis
          
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p));     // Update system energy
    aam.save("rb.aam",cell.p);          // Save particle configuration to disk
    cout << loop.timing()<<sys.info();              // Show progres
  }                                     // END of macro loop and simulation

  io.writefile("gcgroup.conf", nmt.print());

  cout << cell.info() << sys.info() << sm.info() << sb.info()
       << wid.info() << loop.info();
  iopqr pqr;
  pqr.save("saltbath.pqr", cell.p);
}

