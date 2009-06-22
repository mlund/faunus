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


using namespace std;
using namespace Faunus;

int main() {
  cout << faunus_splash();
  slump slump;                          // A random number generator
  inputfile in("manybodyDH.conf");      // Read input file
  mcloop loop(in);                      // Keep track of time and MC loop
  box cell(in);                         // We want a cubic cell
  canonical nvt;                        // Use the canonical ensemble
  interaction<pot_debyehuckelP3trunk> pot(in);    // Functions for interactions

  iogro gro(cell.atom, in);             // Gromacs file output for VMD etc.

  FAUrdf protrdf(0,0,.5,cell.len/2.);   // Protein and salt radial distributions
  FAUrdf protrdf11(0,0,.5,cell.len/2.); // Protein and salt radial distributions
  FAUrdf protrdf22(0,0,.5,cell.len/2.); // Protein and salt radial distributions
  FAUrdf protrdf12(0,0,.5,cell.len/2.); // Protein and salt radial distributions

  vector<macromolecule> g;              // PROTEIN groups
  ioaam aam(cell.atom);                 //   Protein input file format is AAM
  aam.load(cell, in, g);                //   Load and insert proteins
  g[0].center(cell);                    //   Center first protein (will be frozen)

  if (aam.load(cell,"confout.aam")) {
    for (int i=0;i<g.size();i++) 
      g[i].masscenter(cell);              // Load old config (if present)
                                          // ...and recalc mass centers
  }
  // Markove moves and parameters

  macrorot mr(nvt, cell, pot);          //   Class for macromolecule rotation
  translate mt(nvt, cell, pot);         //   Class for macromolecular translation
  clustertrans ct(nvt, cell, pot, g);   //   Class for non-rejective cluster translation
  mr.dp=in.getflt("mr_dp");
  mt.dp=in.getflt("mt_dp");
  ct.dp=in.getflt("ct_dp");



  // Analysis
  double usys;
  usys=0;
  for (int i=0; i<g.size()-1; i++)
    for (int j=i+1; j<g.size(); j++)
      usys+=pot.energy(cell.p, g[i], g[j]);

  systemenergy sys(usys);               // System energy analysis
  aggregation agg(cell, g, 1.0);
  virial vir(cell);
  vir.dr=in.getflt("virdr", 0.01);

  // Switch parameters

  double sum, t, clt, r, randy;
  t=in.getflt("t", 1.), clt=in.getflt("clt", 1.), r=in.getflt("r", 1.);
  sum=t+clt+r;
  t/=sum, clt/=sum, r/=sum;

  cout << cell.info() << pot.info()     // Print information to screen
       << cell.atom.info()
       << in.info(); 
  #ifdef GROMACS
  ioxtc xtc(cell, cell.len);            // Gromacs xtc output (if installed)
  #endif

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,j,n;
      randy=slump.random_one();                 // Pick a random MC move

      if (randy<clt)                            // Cluster translation
        sys+=ct.move(g);                        //   Do the move.

      if (randy>clt && randy<clt+r)             // Rotate proteins
        for (n=0; n<g.size(); n++) {            //   Loop over all proteins
          i = rand() % g.size();                //   and pick at random.
          //if (i>0)                            //   (freeze 1st molecule)
            sys+=mr.move(g[i]);                 //   Do the move.
        }

      if (randy >clt+r) {                       // Translate proteins
        for (n=0; n<g.size(); n++) {            //   Loop over all proteins
          i = rand() % g.size();                //   and pick at random.
         // if (i>0)                            //   (freeze 1st molecule)
          sys+=mt.move(g[i]);                   //   Do the move.
        }
          
      }
      if (slump.random_one()>.7) {
        for (i=0;i<g.size();i++) 
          g[i].masscenter(cell);                // Recalculate mass centers

      
        for (i=0; i<g.size()-1; i++)
          for (j=i+1; j<g.size(); j++) {        //   Analyse g(r)...
            protrdf.update(cell,g[i].cm,g[j].cm);
            if ((i<in.getint("nprot2") && j>=in.getint("nprot2")) || 
                (i>=in.getint("nprot2") && j<in.getint("nprot2")))
              protrdf12.update(cell,g[i].cm,g[j].cm);
            else {
              if(i<in.getint("nprot2"))
                protrdf11.update(cell,g[i].cm,g[j].cm);
              else
                protrdf22.update(cell,g[i].cm,g[j].cm);
            }
         }

        agg.count();
        vir.sample(cell,pot);
      }
      #ifdef GROMACS
      if (slump.random_one()>.96 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif

    } // End of inner loop


    cout << loop.timing(macro);                 // Show middle time and ETA
    usys=0;
    for (int i=0; i<g.size()-1; i++)
      for (int j=i+1; j<g.size(); j++)
        usys+=pot.energy(cell.p, g[i], g[j]);
    sys.update(usys);                           // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector
    gro.save("confout.gro", cell.p);            // Write GRO output file
    protrdf.write("rdfprot.dat");               // Write g(r)'s
    protrdf11.write("rdfprot11.dat");               // Write g(r)'s
    protrdf12.write("rdfprot12.dat");               // Write g(r)'s
    protrdf22.write("rdfprot22.dat");               // Write g(r)'s
    agg.write("aggregates.dat");

    cout <<"#   Energy  (cur,avg)   = "<<usys<<" "<<sys.uavg.avg()<<"(kT)"<<endl
         <<"#   Drift               = "<<abs(sys.sum-usys)<<"      (kT)"<<endl;


  } // End of outer loop

  cout << "----------- FINAL INFORMATION -----------" << endl ;
  cout << sys.info() << agg.info() <<vir.info()      // Final information...
       << ct.info() << mr.info() << mt.info() << loop.info();

  aam.save("confout.aam", cell.p);            // Save config. for next run
//  xyz.save("confout.xyz", cell.p);

  #ifdef GROMACS
  xtc.close();
  #endif
}

