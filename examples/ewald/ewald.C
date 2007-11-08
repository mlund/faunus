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
typedef pot_minimage T_pairpot;         // Coulomb pot. with minimum image convention.
#include "markovmove.C"
#include "analysis.h"
#include "histogram.h"

using namespace std;

int main() {
  box box(50.);                         // A cubic box - periodic boundaries
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential - default values
  cfg.box = box.len;                    // Specific box len. for minimum image
  rdf rdf(particle::NA,particle::CL);   // Prepare Na-Cl radial distribution, g(r)
  interaction<T_pairpot> pot(cfg);      // Particle energy functions
  ioxyz xyz(box);                       // Class for XYZ output
  saltmove sm(nvt, box, pot);           // Class for salt movements
  widom<T_pairpot> widom(box, pot,
      particle::NA, particle::CL);      // Class for Widom particle insertion
  group salt;                                   // Group for mobile ions
  salt.add( box, particle::NA, 80 );            // Insert sodium ions
  salt.add( box, particle::CL, 80 );            // Insert chloride ions

  systemenergy sys(pot.energy(box.p));          // Track system energy
  
  for (int macro=0; macro<10; macro++) {        // Markov chain
    for (int micro=0; micro<1e2; micro++) {
      sm.move(salt);                            // Displace salt particles
      sm.adjust_dp(40,50);                      // Stride to 40-50% acceptance
      sys+=sm.du;                               // Sum system energy changes
      widom.insert(10);                         // Widom particle insertion analysis
      rdf.update(box.p);                        // Analyse Na-Cl distribution
    }
    sys.update(pot.energy(box.p));              // Update system energy
  }
  xyz.save("coord.xyz", box.p);                 // Write XYZ file
  rdf.write("rdf.dat");                         // Write g(r)
  cout << box.info() << sys.info() << sm.info() << widom.info();
}

