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
  iogro gro(con.atom,in);
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

  // Distribution functions
  FAUrdf saltrdf(con.atom["NA"].id,con.atom["CL"].id,0.2,20.);
  FAUrdf catcat( con.atom["NA"].id,con.atom["NA"].id,0.2,20.);
  FAUrdf spcrdf( con.atom["OW"].id,con.atom["OW"].id,0.2,10.);

  mr.dp=0.6;
  mt.dp=1.0; //1.0
  sm.dp=0.6; //0.6

  aam.load(con, "confout.aam");        // Load old config (if present)
  systemenergy sys(pot.energy(con.p)); // System energy analysis

  cout << in.info()
    << con.info() << con.atom.info()
    << pot.info() << sol.info();

  group head;                            // Group of first three particles [0:2]
  head.set(0,2);                         // (Used to swap w. SPC for parallization reasons)
  macromolecule m;

  while ( loop.macroCnt() ) {            // Markov chain 
    while ( loop.microCnt() ) {
      m=sol[ sol.random() ];
      m.cm=m.cm_trial=con.p[m.beg];
      switch (rand() % 3) {              // Randomly chose move
        case 0:
          sys+=sm.move(salt,1);          // Displace a salt particle
          break;
        case 1:
          if (m.swap(con,head)) {
            sys+=mr.move(m);             // Rotate solvent
            m.swap(con,head);
          }
          break;
        case 2:
          if (m.swap(con,head)) {
            sys+=mt.move(m);             // Translate solvent
            m.swap(con,head);
          }
          break;
      }
      if (slump.random_one()<0.05) {
        saltrdf.update(con, salt);
        catcat.update(con, salt);
      }
    }//end of micro-loop 

    spcrdf.update(con, sol);
    spcrdf.write("rdf-OW-OW.dat");
    catcat.write("rdf-Na-Na.dat");
    saltrdf.write("rdf-Na-Cl.dat");
    sys.update(pot.energy(con.p));       // Update system energy
    gro.save("confout.gro", con.p);      // Save config. to disk
    aam.save("confout.aam", con.p);      // Save config. to disk
    cout << loop.timing();               // Show progress

  }//end of macro-loop
  cout << sys.info()
       << salt.info(con)
       << protein.info() << loop.info()  // Print final results
       << sm.info() << mr.info() << mt.info();
}

