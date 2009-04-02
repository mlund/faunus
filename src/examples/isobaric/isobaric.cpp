/*
 * This program will simulate an arbitrary number of
 * protein molecules in a dielectric solvent with
 * explicit mobile ions.
 *
 * \author Bjorn Persson
 */

#include "faunus/faunus.h"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

using namespace Faunus;
using namespace std;

 
class BadConversion : public std::runtime_error {
public:
  BadConversion(const std::string& s)
    : std::runtime_error(s)
    { }
};
 
inline std::string stringify(double x)
{
  std::ostringstream o;
  if (!(o << x))
    throw BadConversion("stringify(double)");
  return o.str();
}

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
  io io;
  ioxyz xyz(cell.atom);
  ioaam aam(cell.atom);
  if (in.getboo("lattice")==true)         //   Protein input file format is AAM
    aam.loadlattice(cell, in, g);                //   Load and insert proteins
  else                                    //   Center first protein (will be frozen)
    aam.load(cell, in, g);
  if (aam.load(cell,"confout.aam")) {
    for (int i=0;i<g.size();i++) 
      g[i].masscenter(cell);              // Load old config (if present)
                                          // ...and recalc mass centers
  }

  // Markovsteps
  macrorot mr(nvt, cell, pot);            //   Class for macromolecule rotation
  translate mt(nvt, cell, pot);           //   Class for macromolecular translation
  isobaric<pot_debyehuckelP3> vol(
      nvt, cell, pot,
      in.getflt("pressure"),              // Set external pressure (in kT)
      in.getflt("penalty"),               // Set penalty           (in kT)
      in.getflt("scalepen"),              // Scale factor
      int(in.getint("maxlen")),           // Set maximum box length
      int(in.getint("minlen")),           // Set minimum box length
      int(in.getint("binlen")));          // Set bin length for penalty function  

  if ( in.getboo("penalize")==true)       // If penalty.dat is present it will be loaded and used as
   vol.loadpenaltyfunction("penalty.dat");// bias. Else a penalty function will only be constructed
                                          // if the penalty != 0.
  // Markovparameters
  vol.dp=in.getflt("voldp");
  mt.dp =in.getflt("mtdp");
  mr.dp =in.getflt("mrdp");

  // Analysis and energy
  systemenergy sys(pot.energy(cell.p));   // System energy analysis
  histogram lendist(1. ,in.getflt("minlen"), in.getflt("maxlen"));             

  FAUrdf pp(cell.atom["LYS"].id, cell.atom["LYS"].id, 1., 100.); 
  FAUrdf pm(cell.atom["LYS"].id, cell.atom["ASP"].id, 1., 100.); 
  FAUrdf mm(cell.atom["ASP"].id, cell.atom["ASP"].id, 1., 100.); 

  // Switch parameters
  double sum, volr, tr, rr, randy=-1;
  volr = in.getflt("volr"), tr=in.getflt("tr"), rr=in.getflt("rr");
  sum  = volr+tr+rr, volr/=sum, tr/=sum, rr/=sum;
  int switcher=-1;

  // Help variables
  int eprint, cnt=0;
  eprint=int(0.001*in.getflt("microsteps")); //loop.macro);

  #ifdef GROMACS
  ioxtc xtc(cell, cell.len);              // Gromacs xtc output (if installed)
  #endif

  cout << cell.info() << pot.info();      // Print information to screen

  cout <<endl<< "#  Temperature = "<<phys.T<<" K"<<endl<<endl;
  cout << "---------- RUN-TIME INFORMATION  -----------" << endl;

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,j,n;
      cnt++;
      randy=slump.random_one();
      if (randy<volr)                           // Pick a random MC move
        sys+=vol.move(g);                       //   Do the move.
      if (randy<volr+rr && randy >volr)         // Rotate proteins
        for (n=0; n<g.size(); n++) {            //   Loop over all proteins
          i = slump.random_one()*g.size();      //   and pick at random.
          sys+=mr.move(g[i]);                   //   Do the move.
        }
      if (randy<volr+rr+tr && randy>volr+rr)    // Translate proteins
        for (n=0; n<g.size(); n++) {            //   Loop over all proteins
          i = slump.random_one()*g.size();      //   and pick at random.
          sys+=mt.move(g[i]);                   //   Do the move.
          for (j=0; j<g.size(); j++)            //   Analyse g(r)...
            if (j!=i && macro>1)
              protrdf.update(cell,g[i].cm,g[j].cm);
          }
      
      if (macro==-1 && micro<1e3) {
        mt.adjust_dp(30,40);                    // Adjust displacement
        mr.adjust_dp(40,50);                    // parameters. Use ONLY
        vol.adjust_dp(20,30);                   // during equillibration!
      }
      lendist.add(cell.len);
      #ifdef GROMACS
      if (slump.random_one()>.80 && macro>1)
        xtc.save("ignored-name.xtc", cell.p);   // Save trajectory
      #endif
      if (cnt%eprint==0)
        sys.track();
      if(slump.random_one()>.99) {
        pp.update(cell);
        pm.update(cell);
        mm.update(cell);
      }
    } // End of inner loop

    cout << loop.timing(macro);
    cout << "#   Energy (cur, avg, std)    = "<<sys.cur<<", "<<sys.uavg.avg()<<", "<<sys.uavg.stdev()<<endl;
    sys.update(pot.energy(cell.p));             // Update system energy averages
    cell.check_vector();                        // Check sanity of particle vector
    gro.save("confout.gro", cell);              // Write GRO output file
    protrdf.write("rdfprot.dat");               // Write g(r)'s
    if (in.getboo("penalize")==true) {
      vol.printpenalty("oldpenalty.dat");       // Print penalty function (used as bias)
      vol.printupdatedpenalty("penalty.dat");   // Print updated penalty function (not used as bias)
    }
    lendist.write("length-distribution.dat");   // Write volume distribution
    sys.write();
    if (in.getboo("penaltyupdate")==true) {     // Dont use for production
      vol.updatepenalty();
      cout << "# Penalty function updated"<<endl;
    }
    // Update inputfile isobaric.conf and print coordinates
    in.updateval("boxlen", stringify(cell.len));
    io.writefile("isobaric.conf", in.print());
    aam.save("confout.aam", cell.p);

    pp.write("pp-rdf.dat");
    pm.write("pm-rdf.dat");
    mm.write("mm-rdf.dat");

  } // End of outer loop
  
  if (in.getboo("penalize")==false) {
    vol.printupdatedpenalty("penalty.dat");
  }

  cout << "----------- FINAL INFORMATION -----------" << endl ;
  cout << loop.info() << sys.info() << vol.info()             // Final information...
       << mr.info() << mt.info();
  cout <<endl << "#   Final      side length  = " <<cell.len<<endl
       << "#   Ideal     <side length> = " <<pow( double( g.size() )/ in.getflt("pressure" ),1./3.)<<endl
       << "#   Simulated <side length> = " <<vol.len.avg()<<" ("<<vol.len.stdev()<<")"<<endl
       << "#   Ideal     <density>     = " <<in.getflt("pressure")<<endl
       << "#   Simulated <density>     = " <<g.size()*vol.ivol.avg()<<" ("<<g.size()*vol.ivol.avg()<<")"<<endl;
  aam.save("confout.aam", cell.p);            // Save config. for next run
  xyz.save("confout.xyz", cell.p);


  #ifdef GROMACS
  xtc.close();
  #endif
}

