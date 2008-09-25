/*! \page test_widom Widom
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method.
 *
 * \note: No equilibration run is incorporated.
 * \author Mikael Lund
 * \date Dejvice, 2007
 * \include widom.cpp
 */
#include "faunus/faunus.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  slump slump;
  inputfile in("widom.conf");           // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  cell cell(in);                        // We want a spherical cell
  canonical nvt;                        // Use the canonical ensemble
  FAUrdf rdf(particle::NA,particle::CL,.5, 45.);
  interaction<pot_datapmf> pot(in);     // Pair potential
  widom widom(cell, pot,
      particle::NA, particle::CL);      // Class for Widom particle insertions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  salt salt;                            // Define some groups for mobile ions
  salt.add( cell, in );                 // Insert some sodium ions

  pot.pair.loadpmf(cell);               // Load pmf's for the particles in the system
  cout << pot.pair.info(cell);
  systemenergy sys(pot.energy(cell.p)); // Track system energy

  while ( loop.macroCnt() ) {           //Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
      widom.insert(10);                 // Widom particle insertion analysis
      if (slump.random_one()>0.9)
        rdf.update(cell);               // Update g-of-r
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p));     // Update system energy
    cout << loop.timing();              // Show progres;s
  }                                    // END of macro loop and simulation

  rdf.write("rdf.dat");                 // Write g-of-r to disk

  cout << cell.info() << sys.info() << loop.info()
    << sm.info() << widom.info();       // Print information!
}

