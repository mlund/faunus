/*
 * This program will simulate an arbitrary number of
 * protein molecules in a dielectric solvent with
 * explicit mobile ions.
 *
 * \author Bjorn Persson
 */

#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main() {
  cout << "---------- INITIAL PARAMETERS -----------" << endl;
  slump slump;                            // A random number generator
  physconst phys;
  inputfile in("isobaric.conf");          // Read input file
  phys.e_r=in.getflt("e_r");
  phys.lB_TO_T(in.getflt("bjerrum"));
  mcloop loop(in);                        // Keep track of time and MC loop
  box cell(in);                           // We want a cubic cell
  canonical nvt;                          // Use the canonical ensemble
  interaction<pot_debyehuckelP3> pot(in); // Functions for interactions
  iogro gro(cell.atom, in);               // Gromacs file output for VMD etc.
  FAUrdf protrdf(0,0,.5,cell.len/2.);     // Protein and salt radial distributions

  vector<macromolecule> g;                // PROTEIN groups
  ioxyz xyz(cell.atom);
  ioaam aam(cell.atom);
  if (in.getflt("lattice")==true)         //   Protein input file format is AAM
    aam.loadlattice(cell, in, g);                //   Load and insert proteins
  else                                    //   Center first protein (will be frozen)
    aam.load(cell, in, g);
  if (aam.load(cell,"confout.aam")) {
    for (int i=0;i<g.size();i++) 
      g[i].masscenter(cell);              // Load old config (if present)
                                          // ...and recalc mass centers
  }

  macrorot mr(nvt, cell, pot);            //   Class for macromolecule rotation
  translate mt(nvt, cell, pot);           //   Class for macromolecular translation
  systemenergy sys(pot.energy(cell.p));   // System energy analysis
  isobaric<pot_debyehuckelP3> vol(
      nvt, cell, pot,
      in.getflt("pressure"),
      in.getflt("penalty"),
      int(in.getflt("max")) );
  histogram lendist(1,0,in.getflt("max"));           //  
  
  vol.dp=in.getflt("voldp");
  mt.dp=in.getflt("mvdp");
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
            sys+=mr.move(g[i]);                 //   Do the move.
          }
          break;
        case 2:                                 // Translate proteins
          for (n=0; n<g.size(); n++) {          //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            sys+=mt.move(g[i]);                 //   Do the move.
            for (j=0; j<g.size(); j++)          //   Analyse g(r)...
              if (j!=i && macro>1)
                protrdf.update(cell,g[i].cm,g[j].cm);
          }
          break;
      }
      if (macro==0 && micro<1e3) {
        mt.adjust_dp(30,40);                    // Adjust displacement
        mr.adjust_dp(40,50);                    // parameters. Use ONLY
        vol.adjust_dp(20,30);                    // during equillibration!
      }
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
    vol.printpenalty("penalty.dat");
    lendist.write("length-distribution.dat");                  // Write volume distribution

  } // End of outer loop

  cout << "----------- FINAL INFORMATION -----------" << endl ;
  cout << sys.info() << vol.info()             // Final information...
       << mr.info() << mt.info();
  cout << "#   Final side length = " <<cell.len<<endl;
  aam.save("confout.aam", cell.p);            // Save config. for next run
  xyz.save("confout.xyz", cell.p);


  #ifdef GROMACS
  xtc.close();
  #endif
}

