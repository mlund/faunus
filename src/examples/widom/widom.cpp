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
#include "faunus/potentials/pot_coulomb.h"
#include "faunus/moves/markovmove.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  inputfile in("widom.conf");           // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  cell cell(in);                        // We want a spherical cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential - default values
  FAUrdf rdf(particle::NA,particle::CL,.5, 45.);
  interaction<T_pairpot> pot(cfg);      // Energy functions
  widom<T_pairpot> widom(cell, pot,
      particle::NA, particle::CL);      // Class for Widom particle insertions
  widom.runfraction=0.5;
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  salt salt;                            // Define some groups for mobile ions
  salt.add( cell, in );                 // Insert some sodium ions

  systemenergy sys(pot.energy(cell.p)); // Track system energy

  while ( loop.macroCnt()==true ) {           //Markov chain 
    while ( loop.microCnt()==true ) {
      sys+=sm.move(salt);               // Displace salt particles
      widom.insert(10);                 // Widom particle insertion analysis
      rdf.update(cell);                 // Update g-of-r
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p));     // Update system energy
    cout << loop.timing();              // Show progres;s
  }                                    // END of macro loop and simulation

  rdf.write("rdf.dat");                 // Write g-of-r to disk

  cout << cell.info() << sys.info() << loop.info()
    << sm.info() << widom.info();       // Print information!
}

