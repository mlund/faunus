/*!\page test_spce SPCE/E
 *
 * \author Mikael Lund and Bjorn Persson
 * \include spce.cpp
 */
#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  slump slump;
  cout << faunus_splash();             // Faunus spam
  inputfile in("spce.conf");           // Read input file
  cell con(in);                        // Use a spherical simulation container
  ioaam aam(con.atom);                 // Protein input file format is AAM
  iopqr pqr(con.atom);                 // PQR coordinate output
  mcloop loop(in);                     // Set Markov chain loop lengths
  canonical nvt;                       // Use the canonical ensemble
  interaction<pot_test> pot(in);       // Specify pair potential
  pot.pair.init(con.atom);
  macrorot mr(nvt, con, pot);          // Class for molecular rotation
  translate mt(nvt, con, pot);         // Class for molecular translation

  // Load central protein
  macromolecule protein;               // Group for the protein
  protein.add( con,
      aam.load(in.getstr("protein"))); // Load protein structure
  protein.move(con, -protein.cm);      // ..translate it to origo (0,0,0)
  protein.accept(con);                 // ..accept translation

  // Add salt
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions
  saltmove sm(nvt, con, pot);          // Class for salt movements

  // Add water
  molecules sol(3);                    // We want a three point water model
  vector<particle>
    water = aam.load("water.aam");     // Load bulk water from disk (typically from MD)
  con.atom.reset_properties(water);    // Set particle parameters according to Faunus

  sol.add(con,water,sol.numatom);      // Inject water into the cell - avoid salt and protein overlap
  water.clear();                       // Free the (large) bulk water reservoir
  macromolecule m;
  mr.dp=0.1;
  mt.dp=0.3;

  int io=con.p[ sol.beg ].id;
  int ih=con.p[ sol.beg+1 ].id;
  cout <<  con.atom[io].sigma << " " << con.atom[io].eps << endl;
  cout <<  con.atom[ih].sigma << " " << con.atom[ih].eps << endl;

  aam.load(con, "confout.aam");        // Load old config (if present)
  systemenergy sys(pot.energy(con.p)); // System energy analysis

  cout << con.info() << con.atom.info()
       << pot.info() << sol.info();

  while ( loop.macroCnt() ) {            // Markov chain 
    while ( loop.microCnt() ) {
      switch (rand() % 2) {              // Randomly chose move
        case 0:
          sys+=sm.move(salt);            // Displace salt particles
          break;
        case 1:
          for (int i=0; i<sol.size()/sol.numatom; i++) {
            m=sol[ sol.random() ];
            if (slump.random_one()<0.5)
              sys+=mr.move(m);
            else
              sys+=mt.move(m);
          }
          break;
      }
    }                                    // END of micro loop
    sys.update(pot.energy(con.p));       // Update system energy
    aam.save("confout.aam", con.p);      // Save config. to disk
    pqr.save("confout.pqr", con.p);
    cout << loop.timing();               // Show progress
  }                                      // END of macro loop
  cout << sys.info() << sm.info()
       << salt.info(con)
       << protein.info() << loop.info()  // Print final results
       << mr.info() << mt.info();
}

