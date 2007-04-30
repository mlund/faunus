#include <iostream>
#include "../container.h"
#include "../potentials.h"
typedef pot_coulomb T_pairpot;
#include "../markovmove.C"

using namespace std;

int main() {
  cell cell(100.);                      // We want a spherical cell
  canonical nvt;                        // Canonical ensemble
  pot_setup cfg;                        // Setup pair potential
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  saltmove sm(nvt, cell, pot);          // Class for salt movements

  group salt;                                   // Group for mobile ions
  salt+=cell.insert( particle::NA, 10 );        // Insert sodium ions
  salt+=cell.insert( particle::CL, 10 );        // Insert chloride ions

  cout << salt << endl;
  cout << pot.energy(cell.p) << endl;

  // Markov chain
  for (int macro=0; macro<10; macro++) {
    for (int micro=0; micro<1e3; micro++) {
      sm.move(salt, 30.);
    }
  }
};

