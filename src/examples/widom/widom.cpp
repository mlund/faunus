/*! \page test_widom Widom
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method.
 *
 * \author Mikael Lund
 * \date Dejvice, 2007
 * \include widom.cpp
 */
#include "faunus/faunus.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  cout << faunus_splash();              // Show Faunus information
  inputfile in("widom.conf");           // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  canonical nvt;                        // Use the canonical ensemble
#ifdef WIDOM_SPHERE
  cell cell(in);                        // We want a spherical simulation container
  interaction<pot_hscoulomb> pot(in);   // ...and a Coulomb + Hard sphere potential
#else
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_hsminimage> pot(in);  // ...and a Coulomb + Hard sphere potential w. minimum image distances
#endif
  saltmove sm(nvt,cell,pot,in);         // Class for salt movements
  salt salt;                            // Define some groups for mobile ions
  salt.add(cell,in);                    // Insert some ions

  ioaam aam(cell);                      // File I/O class
  aam.load(cell,"widom.aam");           // Read initial config. from disk (if present)

  widom wid1(10);                       // Class for multiple particle insertion
  widomSW wid2(10);                     // Class for single particle insertion w. charge scaling
  wid1.add(cell);                       // Detect all species in the cell
  wid2.add(cell);                       // - // -

  systemenergy sys(pot.energy(cell.p)); // Track system energy

  cout << cell.info() << pot.info()
       << salt.info(cell)
       << in.info();                    // Print initial information

  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
      wid1.insert(cell,pot);            // Widom particle insertion analysis
      wid2.insert(cell,pot);            // - // -
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p));     // Update system energy
    aam.save("widom.aam",cell.p);       // Save particle configuration to disk
    cout << loop.timing();              // Show progres
  }                                     // END of macro loop and simulation

  cout << sys.info() << sm.info()
       << wid1.info() << wid2.info()
       << loop.info();                  // Final information and results!
}

