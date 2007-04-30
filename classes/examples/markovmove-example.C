#include <iostream>
#include "../widom.h"
#include "../container.h"
#include "../potentials.h"
typedef pot_coulomb T_pairpot;
#include "../markovmove.C"

using namespace std;

int main() {
  cell cell(50.);                      // We want a spherical cell
  canonical nvt;                        // Canonical ensemble
  pot_setup cfg;                        // Setup pair potential
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  saltmove sm(nvt, cell, pot);          // Class for salt movements
  widom<T_pairpot> widom(cell, pot,
      particle::NA, particle::CL);      // Class for Widom particle insertion

  group salt;                                   // Group for mobile ions
  salt+=cell.insert( particle::NA, 10 );        // Insert sodium ions
  salt+=cell.insert( particle::CL, 10 );        // Insert chloride ions

  // Markov chain
  for (int macro=0; macro<10; macro++) {
    for (int micro=0; micro<1e3; micro++) {
      sm.move(salt, 30.);
      widom.insert(10);
    }
  }
  cout << "NaCl mean activity coefficient = " << widom.gamma() << endl;
};

