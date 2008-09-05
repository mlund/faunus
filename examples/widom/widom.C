/*! \file Example how to simulate a NaCl solution
 *  \example widom-example.C
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method.
 * Note: No equilibration run is incorporated.
 */
#include <iostream>
#include "io.h"
#include "analysis.h"
#include "container.h"
#include "potentials.h"
namespace Faunus {
  typedef pot_minimage T_pairpot;
}
#include "markovmove.h"
#include "analysis.h"
#include "histogram.h"

using namespace std;
using namespace Faunus;

int main() {
  box cell(90.);                        // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential - default values
  FAUrdf rdf(particle::NA,particle::CL,.5, 45.);
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  ioxyz xyz(cell);
  widom<T_pairpot> widom(cell, pot,
      particle::NA, particle::CL);      // Class for Widom particle insertion

  group salt;                                   // Group for mobile ions
  salt.add( cell, particle::NA, 60 );           // Insert sodium ions
  salt.add( cell, particle::CL, 60 );           // Insert chloride ions
  systemenergy sys(pot.energy(cell.p));         // Track system energy

  xyz.save("coord.xyz", cell.p);
  
  for (int macro=0; macro<10; macro++) {        // Markov chain
    for (int micro=0; micro<1e3; micro++) {
      sys+=sm.move(salt);                       // Displace salt particles
      widom.insert(10);                         // Widom particle insertion analysis
      rdf.update(cell);
    }
    sys.update(pot.energy(cell.p));             // Update system energy
  }
  rdf.write("rdf.dat");
  cout << cell.info() << sys.info() << sm.info() << widom.info();
}

