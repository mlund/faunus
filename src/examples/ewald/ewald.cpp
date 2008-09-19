/*! \file Example how to simulate a NaCl solution
 *  \example ewald-example.C
 *  \author Martin Trulsson?
 *  \todo Not yet done -- Martin?
 *
 */
#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main() {
  box box(50.);                         // A cubic box - periodic boundaries
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg;                        // Setup pair potential - default values
  cfg.box = box.len;                    // Specific box len. for minimum image
  FAUrdf rdf(particle::NA,particle::CL);// Prepare Na-Cl radial distribution, g(r)
  interaction<pot_minimage> pot(cfg);   // Particle energy functions
  ioxyz xyz(box);                       // Class for XYZ output
  saltmove sm(nvt, box, pot);           // Class for salt movements
  widom widom(box, pot,
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

