/*! \file Example how to simulate a NaCl solution
 *  \example markovmove-example.C
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method. An output POVRAY file will
 * be generated. Note: No equilibration run is incorporated.
 */
#include <iostream>
#include "../io.h"
#include "../widom.h"
#include "../container.h"
#include "../potentials.h"
typedef pot_coulomb T_pairpot;
#include "../markovmove.C"

using namespace std;

int main() {
  cell cell(90.);                       // We want a spherical cell
  iopov povray(cell);                   // We want a POVRAY snapshot
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  widom<T_pairpot> widom(cell, pot,
      particle::NA, particle::CL);      // Class for Widom particle insertion

  group salt;                                   // Group for mobile ions
  salt+=cell.insert( particle::NA, 10 );        // Insert sodium ions
  salt+=cell.insert( particle::CL, 10 );        // Insert chloride ions

  for (int macro=0; macro<10; macro++) {        // Markov chain
    for (int micro=0; micro<1e4; micro++) {
      sm.move(salt);                            // Displace salt particles
      widom.insert(10);                         // Widom analysis
    }
  }
  cout << sm.info() << widom.info();
  povray.save("markov-example.pov", cell.p);    // Save POVRAY file
};

