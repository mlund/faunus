/*! \file Example how to simulate a NaCl solution
 *  \example widom-example.C
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method.
 * Note: No equilibration run is incorporated.
 */
#include "faunus/io.h"
#include "faunus/analysis.h"
#include "faunus/container.h"
#include "faunus/potentials.h"
namespace Faunus {
  typedef pot_coulomb T_pairpot;       // Specify pair potential
}
#include "faunus/markovmove.h"
#include "faunus/matubayashi.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  box cell(90.);                        // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential - default values
  interaction<T_pairpot> pot(cfg);      // Energy functions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  widom<T_pairpot> widom(cell, pot,
      particle::NA, particle::CL);      // Class for Widom particle insertions
  group salt;
  salt.add( cell, particle::NA, 60 );   // Insert some sodium ions
  salt.add( cell, particle::CL, 60 );   // Insert some chloride ions
  systemenergy sys(pot.energy(cell.p)); // Track system energy
  for (int macro=0; macro<10; macro++) {        // Markov chain
    for (int micro=0; micro<1e2; micro++) {
      sys+=sm.move(salt);                       // Displace salt particles
      widom.insert(10);                         // Widom particle insertion analysis
    }                                           // END of micro loop
    sys.update(pot.energy(cell.p));             // Update system energy
  }                                             // END of macro loop and simulation
  cout << cell.info() << sys.info()
    << sm.info() << widom.info();               // Print information!
}

