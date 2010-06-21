/*! \page test_widom Widom
 *
 * This example program calculates the excess chemical
 * potential of NaCl in an aqueous solution using Widom's
 * particle insertion method.
 *
 * \author Mikael Lund
 * \date Dejvice, 2007
 * \include widom.cpp
 */
#include "faunus/faunus.h"
#include "faunus/potentials/pot_minimage.h"

using namespace std;
using namespace Faunus;                  // Access to Faunus classes

int main(int argc, char* argv[]) {
  cout << faunus_splash();               // Show Faunus information
  string config = "wp_gaussian.conf";    // Default input (parameter) file
  if (argc==2) config = argv[1];         // ..also try to get it from the command line
  inputfile in(config);                  // Read input file
  mcloop loop(in);                       // Set Markov chain loop lengths
  canonical nvt;                         // Use the canonical ensemble
  checkValue test(in);                   // Test output

  slit cell(in);                         // We want a spherical simulation container
  interaction<pot_r12minimageXY> pot(in);// ...and a Coulomb + Hard sphere potential

  saltmove sm(nvt,cell,pot,in);          // Class for salt movements
  salt salt;                             // Define some groups for mobile ions
  salt.add(cell,in);                     // Insert some ions
  
  atom["NA"].variance=0.2;
  atom["NA"].mean=1.0;
  chargeregGaussian tit(nvt,cell,pot);   // Class for charge fluctuations

  ioaam aam;                             // File I/O class
  cell.loadFromDisk("wp_gaussian.dump"); // Read initial config. from disk (if present)
  systemenergy sys(pot.energy(cell.p)+ tit.gaussianEnergyTotal(cell.p));  // Track system energy

  cout << cell.info() << atom.info()
    << pot.info() << salt.info(cell)
    << in.info();                        // Print initial information

  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      switch (rand() % 2) {
        case 0:
          sys+=sm.move(salt);               // Displace salt particles
          break;
        case 1:
          sys+=tit.titrateall();
          break;
      }
    } // END of micro loop
    
    sys.update(pot.energy(cell.p) + tit.gaussianEnergyTotal(cell.p));      // Update system energy
    cell.saveToDisk("wp_gaussian.dump"); // Save particle configuration to disk
    cout << loop.timing();               // Show progres
  }                                      // END of macro loop and simulation

  cout << sys.info() << sm.info()
       << loop.info() << tit.info();     // Final information and results!

  // Now do some testing!
  sm.check(test);
  sys.check(test);
  cout << test.report();
  return test.returnCode();
}

