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
  inputfile in("dipolarfluid.conf");      // Read input file
  hypersphere con(in);                    // We want a hypersphere
  canonical nvt;                          // Use the canonical ensemble
  interaction<pot_hypersphere> pot(in);   // 
  mcloop loop(in);                        // Keep track of time and MC loop

  hypergroup solvent;                     // Group for dipoles

  cout << con.info();

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

