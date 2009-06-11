/*! \page test_wallprotein Wallprotein
 *
 * Simulate a number of flexible polymers in a salt
 * solution.
 *
 * \author Mikael Lund
 * \date Lund, 2009
 * \include wp.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  // General setup
  cout << faunus_splash();
  inputfile in("wp.conf");
  mcloop loop(in);
  canonical nvt;
  box cell(in);
  springinteraction<pot_hsminimage> pot(in);

  // Handle polymers
  polymer pol;
  pol.babeladd( cell, in );
  monomermove mm(nvt,cell,pot,in);
  cout << pol.info();

  // Handle wall particles
  group wall;
  //wall.add(cell, cell.atom["NA"].id, 20);
  saltmove wm(nvt,cell,pot,in);
  wm.dpv.z=0;

  // Handle salt particles
  salt salt(1,1); 
  salt.add(cell,in);           
  saltmove sm(nvt,cell,pot,in);

  // File I/O
  iopqr pqr;
  ioaam aam;
  //aam.load(cell,"conf.aam");

  systemenergy sys(
        pot.internal(cell.p, salt)  +
        pot.energy(cell.p, salt, pol) +
        pot.uself_polymer(cell.p,pol) );

  cout << cell.info() << atom.info()
       << pot.info() << salt.info(cell)
       << in.info();

  while ( loop.macroCnt() ) {           // Markov chain 
    while ( loop.microCnt() ) {
      sys+=sm.move(salt);               // Displace salt particles
      sys+=mm.move(pol);                // Displace monomers
    }                                   // END of micro loop
    sys.update(
        pot.internal(cell.p, salt) +
        pot.energy(cell.p, salt, pol) +
        pot.uself_polymer(cell.p,pol) );     // Update system energy
    aam.save("conf.aam",cell.p);       // Save particle configuration to disk
    pqr.save("conf.pqr",cell.p);
    cout << loop.timing();              // Show progres
  }                                     // END of macro loop and simulation

  cout << sys.info() << sm.info() << mm.info()
       << loop.info(); // Final information and results!
}

