/*! \page test_widom Widom
 *
 * bla...
 *
 * \author Mikael Lund
 * \date Lund, 2009
 * \include wp.cpp
 */
#include "faunus/faunus.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  cout << faunus_splash();              // Show Faunus information
  inputfile in("wp.conf");              // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  canonical nvt;                        // Use the canonical ensemble
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_hsminimage> pot(in);  // ...and a Coulomb/HS pot. w. minimum image
  saltmove sm(nvt,cell,pot,in);         // Class for salt movements
  salt salt;                            // Define some groups for mobile ions
  salt.add(cell,in);                    // Insert some ions

  polymer pol;
  pol.babeladd( cell, in );
  cout << pol.info() ;
  return 0;
  iopqr pqr(cell.atom);                 // File I/O class
  ioaam aam(cell.atom);                 // File I/O class
  aam.load(cell,"conf.aam");           // Read initial config. from disk (if present)

  systemenergy sys(pot.energy(cell.p)); // Track system energy

  cout << cell.info() << cell.atom.info()
       << pot.info() << salt.info(cell)
       << in.info();                    // Print initial information

  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p));     // Update system energy
    aam.save("conf.aam",cell.p);       // Save particle configuration to disk
    pqr.save("conf.pqr",cell.p);
    cout << loop.timing();              // Show progres
  }                                     // END of macro loop and simulation

  cout << sys.info() << sm.info()
       << loop.info(); // Final information and results!
}

