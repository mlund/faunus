/*
 * \author Mikael Lund
 * \date Canberra 2008
 */

#include "analysis.h"
#include "mcloop.h"
#include "pot_hydrophobic.h"
typedef pot_coulomb T_pairpot;      // Specific pair interaction function
#include "markovmove.h"

using namespace std;

int main() {
  cout << "---------- INITIAL PARAMETERS -----------" << endl;
  slump slump;                          // A random number generator
  inputfile in("binding.conf");         // Read input file
  mcloop loop(in);                      // Keep track of time and MC loop
  cell cell(in);                         // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg(in);                    // Setup pair potential (default values)
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  iogro gro(cell, in);                  // Gromacs file output for VMD etc.
  rdf protrdf(0,0,.5,cell.r);           // Protein and salt radial distributions
  twostatebinding bind(20.);            // Two state binding analysis

  vector<macromolecule> g;              // PROTEIN groups
  ioaam aam(cell);                      //   Protein input file format is AAM
  aam.load(cell, in, g);                //   Load and insert proteins
  g[0].center(cell);                    //   Center first protein (will be frozen)
  macrorot mr(nvt, cell, pot);          //   Class for macromolecule rotation
  translate mt(nvt, cell, pot);         //   Class for macromolecular translation
  salt salt;                            // SALT group
  salt.add(cell, in);                   //   Add salt particles
  saltmove sm(nvt, cell, pot);          //   Class for salt movements
  systemenergy sys(pot.energy(cell.p)); // System energy analysis

  cout << cell.info() << pot.info();    // Print information to screen

  #ifdef GROMACS
  ioxtc xtc(cell, cell.len);            // Gromacs xtc output (if installed)
  #endif

  cout << "---------- RUN-TIME INFORMATION  -----------" << endl;
  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,j,n;
      switch (rand() % 3) {                     // Pick a random MC move
        case 0:                                 // Displace salt
          sys+=sm.move(salt);                   //   Do the move.
          break;
        case 1:                                 // Rotate proteins
          for (n=0; n<g.size(); n++) {          //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            if (i>0)                            //   (freeze 1st molecule)
              sys+=mr.move(g[i]);               //   Do the move.
          }
          break;
        case 2:                                 // Translate proteins
          for (n=0; n<g.size(); n++) {          //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            if (i>0)                            //   (freeze 1st molecule)
              sys+=mt.move(g[i]);               //   Do the move.
            for (j=0; j<g.size(); j++)          //   Analyse g(r)...
              if (j!=i && macro>1)
                protrdf.update(cell,g[i].cm,g[j].cm);
          }
          break;
      }
      bind.update(cell, cell.p[g[0].end], g, 1);

      #ifdef GROMACS
      if (slump.random_one()>.96 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif

    } // End of inner loop

    cout << loop.timing(macro);                 // Show middle time and ETA
    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector
    gro.save("confout.gro", cell.p);            // Write GRO output file

  } // End of outer loop

  cout << "----------- FINAL INFORMATION -----------" << endl;
  cout << sys.info() << salt.info(cell)             // Final information...
       << sm.info() << mr.info() << mt.info() << g[0].info(cell) << g[1].info(cell)
       << bind.info( 1./cell.getvolume() );

  #ifdef GROMACS
  xtc.close();
  #endif
}

