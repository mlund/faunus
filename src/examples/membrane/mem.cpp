/*!\page test_manybody Manybody
 * This program will simulate an arbitrary number of
 * protein molecules in a dielectric solvent with
 * explicit mobile ions.
 *
 * \author Mikael Lund
 * \date Lund 2009
 * \include mem.cpp
 */

#include "faunus/faunus.h"
#include "faunus/potentials/pot_minimage.h"

using namespace std;
using namespace Faunus;

int main() {
  cout << faunus_splash();
  slump slump;                          // A random number generator
  inputfile in("mem.conf");             // Read input file
  mcloop loop(in);                      // Keep track of time and MC loop
  box cell(in);                         // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  interaction<pot_minimage> pot(in);    // Functions for interactions

  vector<macromolecule> g;              // PROTEIN groups
  ioaam aam;                            //   Protein input file format is AAM
  aam.load(cell, in, g);                //   Load and insert proteins
  macrorot mr(nvt, cell, pot);          //   Class for macromolecule rotation
  translate mt(nvt, cell, pot);         //   Class for macromolecular translation

  salt salt;                            // SALT group
  salt.add(cell, in);                   //   Add salt particles
  saltmove sm_bulk(nvt, cell, pot);     //   Class for salt movements

  //salt wallsalt;                        // IONS in surface
  saltmove sm_wall(nvt, cell, pot);
  sm_wall.dpv.z=0;                      //   Move only in xy
  //for (int i=wallsalt.beg; i<=wallsalt.end; i++)
  //  pot.p[i].z=0;


  systemenergy sys(pot.energy(cell.p)); // System energy analysis

  cout << cell.info() << pot.info()     // Print information to screen
    << atom.info();

  while (loop.macroCnt()) {                     //Markov chain 
    while (loop.microCnt()) {
      short i,j,n;
      switch (rand() % 4) {                     // Pick a random MC move
        case 0:                                 // Displace salt
          sys+=sm_bulk.move(salt);              //   Do the move.
          break;
        case 1:
          sys+=sm_wall.move(salt);
          break;
        case 2:                                 // Rotate proteins
          for (n=0; n<g.size(); n++) {          //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            sys+=mr.move(g[i]);                 //   Do the move.
          }
          break;
        case 3:                                 // Translate proteins
          for (n=0; n<g.size(); n++) {          //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            if (i>0)                            //   (freeze 1st molecule)
              sys+=mt.move(g[i]);               //   Do the move.
          }
          break;
      }

    } // End of inner loop

    //cout << loop.timing(macro);                 // Show middle time and ETA
    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector

  } // End of outer loop

  cout << "----------- FINAL INFORMATION -----------" << endl ;
  cout << sys.info() << salt.info(cell)             // Final information...
    << sm_bulk.info() << sm_wall.info() << mr.info() << mt.info();
}

