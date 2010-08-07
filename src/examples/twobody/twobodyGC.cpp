/*!\page test_twobody Twobody
 * This program will simulate two protein molecules
 * in a dielectric solvent with explicit mobile ions.
 * The container is SPHERICAL w. hard walls.
 *
 * \author Mikael Lund
 * \date Prague 2007
 * \todo Maybe use a cylindrical cell?
 * \include twobodyGC.cpp
 */

#include "faunus/faunus.h"
#include "faunus/potentials/pot_coulomb.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  cout << faunus_splash();              // Faunus info

  // Input
  string config = "twobody.conf";       // Default input (parameter) file
  if (argc==2) config = argv[1];        // ..also try to get it from the command line
  inputfile in(config);                 // Read input file

  checkValue test(in);                  // Test output
  mcloop loop(in);                      // Set Markov chain loop lengths
  cell cell(in);                        // We want a spherical, hard cell
  grandcanonical nmt;
  interaction<pot_coulomb> pot(in);     // Functions for interactions
  distributions dst;                    // Distance dep. averages
  atomicRdf gofr_pp(0.5,200);           // Amino acid rdf between the two proteins

  // IO
  io io;
  iopqr pqr;                            // PQR output (pos, charge, radius)
  ioxtc xtc(1000.);                     // Gromacs xtc output
  ioaam aam;                            // Protein input file format is AAM

  float gofr_pp_rf=in.getflt("gofr_pp_rf",0);// How often should gofr_pp be sampled?

  vector<macromolecule> g;              // Group for proteins

  // Protein setup
  macrorot mr(nmt, cell, pot);          // Class for macromolecule rotation
  dualmove dm(nmt, cell, pot);          // Class for 1D macromolecular translation
  dm.setup(in);
  dm.load( in, g, 80.);                 // Load proteins and separate them 
  vector<group> vg;
  vg.push_back(g[0]);
  vg.push_back(g[1]);

  // Salt setup
  salt salt;                            // Group for mobile ions
  salt.add(cell, in);                   //   Add salt particles
  saltmove sm(nmt, cell, pot, in);      // Class for salt movements

  // Titration and GC moves
  saltbath sb(nmt,cell,pot,in,salt);    // Class for salt movements
  GCchargereg tit(nmt,cell,pot,in);   
  tit.runfraction=0.5;

  if(nmt.load(cell, "gcgroup.conf")==true) {                  
    aam.load(cell,"confout.aam");       // Read initial config. from disk (if present)
    g[0].masscenter(cell);              // Load old config (if present)
    g[1].masscenter(cell);              // ...and recalc mass centers
  }

  systemenergy sys(
    pot.energy(cell.p, g[0], g[1]) +
    pot.energy(cell.p, salt) +
    pot.internalElectrostatic(cell.p, g[0]) +
    pot.internalElectrostatic(cell.p, g[1]) +
    pot.internal(cell.p, salt));                // System energy analysis

  cout << "# ---- Initial information ----" << endl
       << in.info() << cell.info()
       << tit.info() << pot.info();             // Print information to screen

  cout << "# ---- Runtime output ----" << endl;
  while ( loop.macroCnt() ) {                   // Markov chain 
    while (loop.microCnt() ) {
      switch (rand() % 4) {                     // Pick a random MC move
        case 0:                                 // Salt space!
          if (slp.random_one()>0.5)
            sys+=sm.move(salt);                 // ...translate
          else
            for (int i=0; i<salt.size(); i++)
              sys+=sb.move();                   // ...or GC move
          break;
        case 1:                                 // Rotate proteins
          for (int n=0; n<2; n++) {             //   Loop over all proteins
            int i = rand() % g.size();          //   and pick at random.
            sys+=mr.move(g.at(i));              //   Do the move.
          }
          break;
        case 2:                                 // Translate proteins
          sys+=dm.move(g[0], g[1]);             //   Do the move.
          break;
        case 3:                                 // Fluctuate charges
          sys+=tit.titrateall();                // Titrate sites on the protein
          dst.add("Q1",dm.r, g[0].charge(cell.p));
          dst.add("Q2",dm.r, g[1].charge(cell.p));
          dst.add("MU1",dm.r,g[0].dipole(cell.p));
          dst.add("MU2",dm.r,g[1].dipole(cell.p));
          break;
      }

      if (slp.random_one()>0.95)
        dst.add("Utot", dm.r, pot.energy(cell.p));

      if (slp.random_one()<gofr_pp_rf)
        gofr_pp.update( cell.p, g[0], g[1] );

      if (slp.random_one()>.98 && loop.macro>1)
        xtc.save("coord.xtc", cell.p, vg);      // Save trajectory

    } // End of inner loop
    
    sys.update(
      pot.energy(cell.p, g[0], g[1]) +
      pot.energy(cell.p, salt) +
      pot.internalElectrostatic(cell.p, g[0]) +
      pot.internalElectrostatic(cell.p, g[1]) +
      pot.internal(cell.p, salt));              // System energy analysis
 

    gofr_pp.write("rdfatomic.dat");             // Write interprotein g(r) - Atomic
    dm.gofr.write("rdfprot.dat");               // Write interprotein g(r) - CM
    dst.write("distributions.dat");             // Write other distributions
    aam.save("confout.aam", cell.p);            // Save config. for next run

    tit.applycharges(cell.trial);               // Set average charges on all titratable sites
    pqr.save("confout.pqr", cell.trial);        // ... save PQR file
    pqr.save("confout_ns.pqr", cell.trial,vg);  // ... save PQR file (proteis only)
    cell.trial=cell.p;                          // ... and restore original charges
    io.writefile("gcgroup.conf", nmt.print());

    cell.check_vector();                        // Check sanity of particle vector
    cout << loop.timing() << std::flush;        // Show progress

  } // End of outer loop

  xtc.close();                                  // Close xtc file for writing

  cout << "# ---- Final information ----" << endl
       << loop.info() << salt.info(cell)
       << sm.info() << mr.info() << dm.info()
       << sys.info() << g[0].info() << g[1].info() << cell.info()
       << tit.info() << sb.info();

  // Output tests
  sys.check(test);
  sb.check(test);
  dm.check(test);
  sm.check(test);
  mr.check(test);
  test.check("Protein1_charge", g[0].Q.avg());
  test.check("Protein2_charge", g[1].Q.avg());
  cout << test.report();

  return test.returnCode();
}

