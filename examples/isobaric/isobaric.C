/*
 * This program will simulate an arbitrary number of
 * protein molecules in a dielectric solvent with
 * explicit mobile ions.
 *
 * \author Mikael Lund
 * \date Prague 2007
 */

#include <iostream>
#include "../../classes/io.h"
#include "../../classes/analysis.h"
#include "../../classes/potentials.h"
#include "../../classes/container.h"
#include "../../classes/countdown.h"
#include "../../classes/histogram.h"
#include "../../classes/inputfile.h"
typedef pot_debyehuckelP3 T_pairpot;      // Specific pair interaction function
#include "../../classes/markovmove.h"
#include "../../classes/mcloop.h"
#include "../../classes/physconst.h"


using namespace std;

int main() {
  cout << "---------- INITIAL PARAMETERS -----------" << endl;
  slump slump;                          // A random number generator
  inputfile in("isobaric.conf");        // Read input file
  physconst phys;
  phys.e_r=in.getflt("e_r");
  phys.lB_TO_T(in.getflt("bjerrum"));
  mcloop loop(in);                      // Keep track of time and MC loop
  box cell(in);                         // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  pot_setup cfg(in);                    // Setup pair potential (default values)
  interaction<T_pairpot> pot(cfg);      // Functions for interactions
  iogro gro(cell, in);                  // Gromacs file output for VMD etc.
  rdf protrdf(0,0,.5,cell.len/2.);      // Protein and salt radial distributions
  histogram lendist(1,0,500);          // 5555ram over volume 

  int templendist[1000];        
  for (int p=0;p<1000;p++)
    templendist[p]=0;
  int templencnt=0;

  vector<macromolecule> g;                      // PROTEIN groups
  ioxyz xyz(cell);
  ioaam aam(cell);                      //   Protein input file format is AAM
  aam.load(cell, in, g);                //   Load and insert proteins
  g[0].center(cell);                    //   Center first protein (will be frozen)
  if (aam.load(cell,"confout.aam")) {
    for (int i=0;i<g.size();i++) 
      g[i].masscenter(cell);              // Load old config (if present)
                                          // ...and recalc mass centers
  }

  macrorot mr(nvt, cell, pot);          //   Class for macromolecule rotation
  translate mt(nvt, cell, pot);         //   Class for macromolecular translation
  systemenergy sys(pot.energy(cell.p)); // System energy analysis
  isobaric vol(nvt, cell, pot, in.getflt("pressure"));
//  if (in.getflt("voldp"))
    vol.dp=in.getflt("voldp");

  cout << cell.info();// << pot.info();    // Print information to screen

  #ifdef GROMACS
  ioxtc xtc(cell, cell.len);            // Gromacs xtc output (if installed)
  #endif
  cout << "#    Temperature = "<<phys.T<<" K"<<endl<<endl;
  cout << pot.info();
  cout << "---------- RUN-TIME INFORMATION  -----------" << endl;
  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,j,n;
      switch (int(slump.random_one() * 3 )) {   // Pick a random MC move
        case 0:                                 // Fluctuate the volume
          sys+=vol.move(g);                     //   Do the move.
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
      if (macro==1 && micro<1e3) {
        mt.adjust_dp(30,40);                    // Adjust displacement
        mr.adjust_dp(40,50);                    // parameters. Use ONLY
        vol.adjust_dp(20,30);                    // during equillibration!
      }
      templendist[int(cell.len)]++;
      templencnt++;
      lendist.add(cell.len);
      #ifdef GROMACS
      if (slump.random_one()>.80 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif

    } // End of inner loop

    cout << loop.timing(macro);
    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector
    gro.save("confout.gro", cell);              // Write GRO output file
    protrdf.write("rdfprot.dat");               // Write g(r)'s
    lendist.write("length-distribution.dat");                  // Write volume distribution

  } // End of outer loop

  cout << "----------- FINAL INFORMATION -----------" << endl ;
  cout << sys.info() << vol.info()             // Final information...
       << mr.info() << mt.info();
  cout << "#   Final side length = " <<cell.len<<endl;
  aam.save("confout.aam", cell.p);            // Save config. for next run
  xyz.save("confout.xyz", cell.p);

  ofstream ledist("tempvoldist.dat");
  if (ledist) {
    int s=1000;
    for (int t=0;t<s;t++)
      if (templendist[t]!=0)
        ledist << t <<"  "<< float(templendist[t])/templencnt<< endl;
    ledist.close();
  }

  #ifdef GROMACS
  xtc.close();
  #endif
}

