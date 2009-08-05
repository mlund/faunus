/*
 * This program will simulate two protein molecules
 * in a dielectric solvent with implicit salt
 * using Teixeira's titration scheme.
 * The container is SPHERICAL w. hard walls.
 *
 * \author Andre Teixeira
 * \date Jul 2009
 */
#include "faunus/faunus.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  cout << faunus_splash();              // Faunus info
  slump slump;                          // A random number generator
  string config = "twobodyAT.conf";     // Default input (parameter) file
  if (argc==2) config = argv[1];        // ..also try to get it from the command line
  inputfile in(config);                 // Read input file
  mcloop loop(in);                      // Set Markov chain loop lengths
  cell cell(in);                        // We want a spherical, hard cell
  canonical nvt;                        // Use the canonical ensemble
  interaction<pot_debyehuckel> pot(in);  // Functions for interactions
  ioxyz xyz;                            // xyz output for VMD etc.
  distributions dst;                    // Distance dep. averages
  //  angularcorr  angcorr;                 // Dipole cross-correlation
  iopqr pqr;                            // PQR output (pos, charge, radius)

  vector<macromolecule> g;              // Group for proteins
  dualmove dm(nvt, cell, pot);          //   Class for 1D macromolecular translation
  float init_sep = in.getflt("dm_maxsep", cell.r*1.5)*0.5 + in.getflt("dm_minsep", 0. )*0.5;
  dm.load( in, g, init_sep );           //   Load proteins and separate them
  macrorot mr(nvt, cell, pot);          // Class for macromolecule rotation
  dm.setup(in);                         // Set displacement parameters
  mr.dp=in.getflt("rot_dp");            //  ---//---
  ioaam aam;                            // Protein input file format is AAM
  if (aam.load(cell,"confout.aam")) {
    g[0].masscenter(cell);              // Load old config (if present)
    g[1].masscenter(cell);              // ...and recalc mass centers
  }

  ATchargereg tit(nvt,cell,pot,in.getflt("pH", 7.),in,pot.pair);  // Starting titration
  g[0].conc = in.getflt("ProteinConc", 0.0001) * 0.5;  // Setting concentration of each protein
  g[1].conc = in.getflt("ProteinConc", 0.0001) * 0.5;
  g[0].cm.radius = g[0].vradius(con.p);
  g[1].cm.radius = g[1].vradius(con.p);
  
  systemenergy sys(pot.energy(cell.p)-pot.internal(cell.p,g[0])-pot.internal(cell.p,g[1])); // System energy analysis
  cout << in.info()  << cell.info()
       << tit.info() << pot.info();                   // Print information to screen

  ioxtc xtc(cell.r);                    // Gromacs xtc output (if installed)

  // Calculate atom closests to mass center
  for (int j=0; j<2; j++) {
    double mindist=10000000.;
    int iclosest=-1;
    for (int i=g[j].beg; i<=g[j].end; i++) {
      double dist = cell.p[i].dist( g[j].cm );
      if (dist<mindist) {
        mindist=dist;
        iclosest=i;
      }
    }
    cout << "# Closest particle to CM" << j << " = " << iclosest << ", " << mindist << endl;
  }

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,n;
      switch (rand() % 3) {                     // Pick a random MC move
        case 0:                                 // Rotate proteins
          for (n=0; n<2; n++) {                 //   Loop over all proteins
            i = rand() % g.size();              //   and pick at random.
            sys+=mr.move(g[i]);                 //   Do the move.
          }
          break;
        case 1:                                 // Translate proteins
          sys+=dm.move(g[0], g[1]);             //   Do the move.
          break;
        case 2:                                 // Fluctuate charges
          sys+=tit.titrateall(g);               // Titrate sites on the protein
          if (tit.du!=0) {                      // Average charges and dipoles
            dst.add("Q1",dm.r, g[0].charge(cell.p));
            dst.add("Q2",dm.r, g[1].charge(cell.p));
            dst.add("MU1",dm.r,g[0].dipole(cell.p));
            dst.add("MU2",dm.r,g[1].dipole(cell.p));
          }
          break;
      }
      if (slump.random_one() < .1 ) {
        double z, l0, l1;
        z=g[0].cm.dist(g[1].cm); 
        l0=g[0].mu.len();
        l1=g[1].mu.len();
        g[0].dipole(cell.p);
        g[1].dipole(cell.p);
        dst.add("Left  protein z-dipcomp", z , g[0].mu.z/l0); 
        dst.add("Right protein z-dipcomp", z , g[1].mu.z/l1); 
      }

      if (slump.random_one()>.95 && macro>1 && in.getboo("movie", false)==true)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
    } // End of inner loop

    sys.update(pot.energy(cell.p)-pot.internal(cell.p, g[0])-pot.internal(cell.p, g[1]));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector

    dm.gofr.write("rdfprot.dat");               // Write interprotein g(r)
    dst.write("distributions.dat");             // Write other distributions
    xyz.save("coord.xyz", cell.p);              // Write .xyz coordinate file
    aam.save("confout.aam", cell.p);            // Save config. for next run
    pqr.save("confout.pqr", cell.p);            // ...also save a PQR file
    cout << loop.timing(macro);                 // Show progress
  } // End of outer loop

  cout << mr.info()   << dm.info()
       << tit.info()  << sys.info()
       << g[0].info(cell) << g[1].info(cell)
       << cell.info();

  xtc.close();
}
