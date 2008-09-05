/*! \file Example how to simulate a NaCl solution
 *  \example widom-example.C
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method.
 * Note: No equilibration run is incorporated.
 */
#include "io.h"
#include "analysis.h"
#include "container.h"
#include "potentials.h"
namespace Faunus {
  typedef pot_coulomb T_pairpot;       // Specify pair potential
}
#include "markovmove.h"
#include "matubayashi.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  box cell(90.);                        // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential - default values
  FAUrdf rdf(particle::NA,particle::CL,.5, 45.);
  interaction<T_pairpot> pot(cfg);      // Energy functions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  widom<T_pairpot> widom(cell, pot,
      particle::NA, particle::CL);      // Class for Widom particle insertions
  matubayashi uhist;                    // Class for energy histogram analysis
  systemenergy sys(pot.energy(cell.p)); // Track system energy

  group salt, anion, cation;            // Define some groups for mobile ions
  cation.add(cell, particle::NA,1);     // Add a single cation
  cell.p[cation.beg].clear();           // ...and move it to origo
  anion.add(cell, particle::CL, 1);     // Add a single anion
  cell.p[anion.beg].clear();            // ...and move it to origo
  cell.p[anion.beg].z=10;               // ...and then to somewhere on the z-axis
  cell.trial=cell.p;                    // Sync trial particle vector
  salt.add( cell, particle::NA, 60 );   // Insert some sodium ions
  salt.add( cell, particle::CL, 60 );   // Insert some chloride ions

  for (int macro=0; macro<10; macro++) {        // Markov chain
    for (int micro=0; micro<1e4; micro++) {
      sys+=sm.move(salt);                       // Displace salt particles
      widom.insert(10);                         // Widom particle insertion analysis
      rdf.update(cell);                         // Update g-of-r
      uhist.add(pot,cell.p,anion);              // Update energy histogram
    }                                           // END of micro loop
    sys.update(pot.energy(cell.p));             // Update system energy
  }                                             // END of macro loop and simulation

  rdf.write("rdf.dat");                         // Write g-of-r to disk
  uhist.hist.write("energyhist.dat");           // Write energy histogram to disk

  cout << cell.info() << sys.info()
    << sm.info() << widom.info();               // Print information!
}

