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

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main(int argc, char* argv[]) {
  cout << faunus_splash();              // Show Faunus information
  string config = "widom.conf";         // Default input (parameter) file
  if (argc==2) config = argv[1];        // ..also try to get it from the command line
  inputfile in(config);                 // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  canonical nvt;                        // Use the canonical ensemble
  checkValue test(in);                  // Test output

#ifdef WIDOM_SPHERE
  cell cell(in);                        // We want a spherical simulation container
  interaction<pot_hscoulomb> pot(in);   // ...and a Coulomb + Hard sphere potential
#else
#ifdef WIDOM_ELHCVDW
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_elhcvdw> pot(in);     // ...and a Coulomb/HS/vdW pot. w. minimum image
  pot.pair.init(atom);
#else
#ifdef WIDOM_ELLJ
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_testminim> pot(in);   // ...and a Coulomb + LJ pot. w. minimum image
  pot.pair.init(atom);
#else
  box cell(in);                         // We want a cubic simulation container
  interaction<pot_hsminimage> pot(in);  // ...and a Coulomb/HS pot. w. minimum image
#endif
#endif
#endif

  saltmove sm(nvt,cell,pot,in);         // Class for salt movements
  salt salt;                            // Define some groups for mobile ions
  salt.add(cell,in);                    // Insert some ions

  ioaam aam;                            // File I/O class
  aam.load(cell,"widom.aam");           // Read initial config. from disk (if present)

  virial virial(cell);                  // Virial analysis
  widom wid1(10);                       // Class for multiple particle insertion
  widomSW wid2(10);                     // Class for single particle insertion w. charge scaling
  wid1.add(cell);                       // Detect all species in the cell
  wid2.add(cell);                       // - // -

  systemenergy sys(pot.energy(cell.p)); // Track system energy

  cout << cell.info() << atom.info()
    << pot.info() << salt.info(cell)
    << in.info();                    // Print initial information

  wid1.runfraction = in.getflt("widom_pair_runfraction",1.0);
  wid2.runfraction = in.getflt("widom_single_runfraction",1.0);
  virial.runfraction = in.getflt("virial_runfraction",1.0);

  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
      wid1.insert(cell,pot);            // Widom particle insertion analysis
      wid2.insert(cell,pot);            // - // -
      virial.sample(cell,pot);          // Virial sampling (NOT really appropriate for this pot.!)
    }                                   // END of micro loop
    sys.update(pot.energy(cell.p));     // Update system energy
    aam.save("widom.aam",cell.p);       // Save particle configuration to disk
    cout << loop.timing();              // Show progres
  }                                     // END of macro loop and simulation

  cout << sys.info() << sm.info()
    << wid1.info() << wid2.info()
    << virial.info() << loop.info(); // Final information and results!

  // Now do some testing!
  sm.check(test);
  sys.check(test);
  wid1.check(test);
  virial.check(test);
  cout << test.report();
  return test.returnCode();
}

