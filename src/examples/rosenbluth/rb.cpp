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
  inputfile in("rb.conf");              // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  grandcanonical nmt;                   // Use the canonical ensemble
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_hsminimage> pot(in);  // ...and a Coulomb/HS pot. w. minimum image
  saltmove sm(nmt,cell,pot,in);         // Class for salt movements
  salt salt;                            // Define some groups for mobile ions
  salt.add(cell,in);                    // Insert some ions

  nmt.searchsalt(cell,salt);

  vector<rosenbluth> rb;
  int ndx=1;
  while (ndx<6)
    rb.push_back( 
        rosenbluth(nmt,cell,pot,in,ndx++) );

  ioaam aam;                            // File I/O class
  aam.load(cell,"widom.aam");           // Read initial config. from disk (if present)

  widom wid(10);                        // Class for multiple particle insertion
  wid.add(cell);                        // Detect all species in the cell

  systemenergy sys(pot.energy(cell.p)); // Track system energy

  cout << cell.info() << atom.info()
       << pot.info() << salt.info(cell)
       << in.info();                    // Print initial information

  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
      wid.insert(cell,pot);             // Widom particle insertion analysis
      int i=rand() % rb.size();         // Pick a random...
      sys+=rb[i].move();                // ...Rosenbluth pair
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p));     // Update system energy
    aam.save("rb.aam",cell.p);          // Save particle configuration to disk
    cout << loop.timing();              // Show progres
  }                                     // END of macro loop and simulation

  cout << sys.info() << sm.info()
       << wid.info() << loop.info();
  for (int i=0; i<rb.size(); i++)
    cout << rb[i].info();               // Final information and results!
}

