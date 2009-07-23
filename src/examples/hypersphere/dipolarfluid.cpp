#ifndef HYPERSPHERE
#define HYPERSPHERE
#endif

/*
*/
#include "faunus/faunus.h"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

using namespace Faunus;
using namespace std;

int main() {
  cout << faunus_splash();
  slump slump;                            // A random number generator
  physconst phys;
  inputfile in("dipolefluid.conf");       // Read input file
  hypersphere con(in);                    // We want a cubic cell
  canonical nvt;                          // Use the canonical ensemble
  interaction<pot_coulomb> pot(in);       // 
  mcloop loop(in);                        // Keep track of time and MC loop

  hypergroup solvent;

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      switch (rand() % 3) {
        case 0:            
          break;
        case 1:
          break;
        case 2:
          break;
      }
    } // End of inner loop
  }
}

