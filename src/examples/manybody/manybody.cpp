/*!\page test_manybody Manybody
 * This program will simulate an arbitrary number of
 * protein molecules in a dielectric solvent with
 * explicit mobile ions.
 *
 * \author Mikael Lund
 * \date Prague 2007
 * \include manybody.cpp
 */

#include "faunus/faunus.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/energy/coarsegrain.h"

using namespace std;
using namespace Faunus;

int main() {
  cout << faunus_splash();
  slump slump;                          // A random number generator
  inputfile in("manybody.conf");        // Read input file
  mcloop loop(in);                      // Keep track of time and MC loop
  box cell(in);                         // We want a cubic cell
  grandcanonical nmt;                        // Use the canonical ensemble
#ifndef MONOPOLE
  interaction<pot_minimage> pot(in);    // Functions for interactions
#else
  interaction_monopole<pot_minimage> pot(in); // Far-away monopole approximation?
#endif
  iogro gro(in);                        // Gromacs file output for VMD etc.

  FAUrdf protrdf(0,0,.5,cell.len/2.);   // Protein and salt radial distributions
  FAUrdf saltrdf(atom["NA"].id, atom["SO4"].id, .5, cell.len/2.);

  vector<macromolecule> g;              // PROTEIN groups
  ioaam aam;                            //   Protein input file format is AAM
  iopqr pqr;
  io io;
  aam.load(cell, in, g);                //   Load and insert proteins
  g[0].center(cell);                    //   Center first protein (will be frozen)
  macrorot mr(nmt, cell, pot);          //   Class for macromolecule rotation
  translate mt(nmt, cell, pot, in);     //   Class for macromolecular translation
  salt salt;                            // SALT group
  salt.add(cell, in);                   //   Add salt particles
  saltmove sm(nmt, cell, pot);          //   Class for salt movements


  saltbath sb(nmt,cell,pot,in,salt);    // Class for salt movements
  GCchargereg tit(nmt,cell,pot,in);
  
  if(nmt.load(cell, "gcgroup.conf")==true)
    aam.load(cell,"confout.aam");        // Read initial config. from disk (if present)

  systemenergy sys(pot.energy(cell.p)); // System energy analysis

  cout << cell.info() << pot.info()     // Print information to screen
       << atom.info();
  ioxtc xtc(cell.len);                  // Gromacs xtc output (if installed)

  widomSW wid2(10);                     // Class for single particle insertion w. charge scaling
  wid2.add( atom("NA") );
  wid2.add( atom("CA") );
  wid2.add( atom("LA") );
  wid2.add( atom("CL") );



  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,j,n;
      int ntry;
      switch (rand() % 5) {                     // Pick a random MC move
        case 0:                                 // Displace salt
          sys+=sm.move(salt);                   //   Do the move.
          break;
        case 1:                                 // Rotate proteins
          for (n=0; n<g.size(); n++) {          //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            //if (i>0)                            //   (freeze 1st molecule)
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
        case 3:
          ntry=salt.size();
          for (i=0; i<ntry; i++) {
            sys+=sb.move();                 // Grand Canonical salt move
            if(cell.charge()!=0) {
              cout << "Container is charged by saltbath!"<<endl;
            }
          }
          break;
        case 4:
          sys+=tit.titrateall();
          if(cell.charge()!=0) {
            cout << "Container is charged by gctit!"<<endl;
          }
          break;
      }
      if (macro==1 && micro<1e3) {
        mt.adjust_dp(30,40);                    // Adjust displacement
        mr.adjust_dp(40,50);                    // parameters. Use ONLY
        sm.adjust_dp(20,30);                    // during equillibration!
      }
      if (slump.random_one()>.8 && macro>1)
        saltrdf.update(cell);                   // Update salt g(r)

      if (slump.random_one()>.96 && macro>1)
        wid2.insert(cell,pot);          // sample activity coefficients
  //      xtc.save("ignored-name.xtc", cell.p);   // Save trajectory

    } // End of inner loop

    cout << loop.timing(macro);                 // Show middle time and ETA
    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector
    gro.save("confout.gro", cell.p);            // Write GRO output file
    protrdf.write("rdfprot.dat");               // Write g(r)'s
    saltrdf.write("rdfsalt.dat");               //   -//-
    aam.save("confout.aam", cell.p); 
    pqr.save("confout.pqr", cell.p);

  } // End of outer loop

  cout << "----------- FINAL INFORMATION -----------" << endl ;
  cout << cell.info() << sys.info() << salt.info(cell)             // Final information...
       << sm.info() << mr.info() << mt.info() << sb.info() << tit.info() << wid2.info();

  io.writefile("gcgroup.conf", nmt.print());
  xtc.close();
}

