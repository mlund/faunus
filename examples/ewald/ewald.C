/*! \file Example how to simulate a NaCl solution
 *  \example widom-example.C
 *  \author Martin Trulsson?
 *  \todo Not yet done -- Martin?
 *
 * This could be changed into a cubic box with
 * long range electrostatic corrections
 *
 */
#include <iostream>
#include "io.h"
#include "analysis.h"
#include "container.h"
#include "potentials.h"
typedef pot_test T_pairpot;
#include "markovmove.C"
#include "analysis.h"
#include "histogram.h"

using namespace std;

int main() {
  box cell(50.);                        // We want a cubic cell - periodic boundaries
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential - default values
  rdf rdf(particle::NA,particle::CL);
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
    for (int micro=0; micro<1e2; micro++) {
      sm.move(salt);                            // Displace salt particles
      sm.adjust_dp(40,50);                      // Stride to 40-50% acceptance
      sys+=sm.du;                               // Sum system energy changes
      widom.insert(10);                         // Widom particle insertion analysis
      rdf.update(cell.p);
    }
    sys.update(pot.energy(cell.p));             // Update system energy
  }
  rdf.write("rdf.dat");
  cout << cell.info() << sys.info() << sm.info() << widom.info();
}

