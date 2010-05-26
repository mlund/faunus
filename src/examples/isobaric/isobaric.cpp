/*
 * This program will simulate an arbitrary number of
 * rigid macromolecules in a Debye-Huckel salt solution.
 *
 * \author Bjorn Persson, Anil Kurut, Mikael Lund
 */

#include "faunus/faunus.h"
#include "faunus/energy/coarsegrain.h"
#include "faunus/potentials/pot_debyehuckelP3.h"

using namespace Faunus;
using namespace std;
 
#ifdef MONOPOLE
  typedef interaction_monopole<pot_debyehuckelP3> Tpot;
#elif DIPOLE_CUTOFF
  typedef interaction_dipole<pot_debyehuckelP3> Tpot;
  #define MONOPOLE
#elif defined(FASTDH)
  typedef interaction_vector<pot_debyehuckelP3Fast> Tpot;
#else
  typedef interaction<pot_debyehuckelP3> Tpot;
#endif

string stringify(double x) {
  stringstream o;
  o << x;
  return o.str();
}

int main() {
  cout << faunus_splash();
  cout << "---------- INITIAL PARAMETERS -----------" << endl;
  slump slump;                            // A random number generator
  slump.random_seed(9);
  inputfile in("isobaric.conf");          // Read input file
  checkValue test(in);                    // Class for unit testing
  mcloop loop(in);                        // Keep track of time and MC loop
#ifndef XYPLANE
  box cell(in);                           // We want a cubic cell
#else
  xyplane cell(in);                       // ...or a plane!
#endif
  canonical nvt;                          // Use the canonical ensemble
  Tpot pot(in);                           // Initialize energy functions
  vector<macromolecule> g;                // PROTEIN groups
  io io;
  iogro gro(in);                          // Gromacs file output for VMD etc.
  iopqr pqr;
  ioaam aam;
  if (in.getboo("lattice")==true)         // Protein input file format is AAM
    aam.loadlattice(cell, in, g);         // Load and insert proteins
  else                                    // Center first protein (will be frozen)
    aam.load(cell, in, g);
  if ( cell.loadFromDisk("confout.dump") ) { // Load old config (if present)
    for (int i=0;i<g.size();i++)             // ...and recalc mass centers
      g[i].masscenter(cell);
    cout << "# Config. is read from confout.dump" << endl;
  }  

  // Markovsteps
  macrorot mr(nvt, cell, pot);            //   Class for macromolecule rotation
  translate mt(nvt, cell, pot, in);       //   Class for macromolecular translation
  transrot mtr(nvt, cell, pot);           //   Class for simultaneous macromolecular translation and rotation
  transrot mtrL(nvt, cell, pot);          //   Class for simultaneous macromolecular translation and rotation
  clustertrans ct(nvt, cell, pot, g);     //   Class for non-rejective cluster translation
  clusterrotate cr(nvt, cell, pot);       //   Class for flux-rejective cluster rotation
  cr.sep=in.getflt("cluster_sep",2.0);    //   Set cluster definition (move to constructor)
  ct.skipEnergyUpdate=in.getboo("skipEnergyUpdate",false);

  isobaricPenalty<Tpot> vol(
      nvt, cell, pot,
      in.getflt("pressure"),              // Set external pressure (in kT)
      in.getflt("penalty"),               // Set penalty           (in kT)
      in.getflt("scalepen"),              // Scale factor
      int(in.getint("maxlen")),           // Set maximum box length
      int(in.getint("minlen")),           // Set minimum box length
      int(in.getint("binlen")));          // Set bin length for penalty function  

  if ( in.getboo("penalize")==true)       // If penalty.dat is present it will be loaded and used as
    vol.loadpenaltyfunction("penalty.dat");// bias. Else a penalty function will only be constructed

  // Markov parameters
  bool movie=in.getboo("movie",false);
  double dppercent = in.getflt("dppercent", 0.00080);
  vol.dp=in.getflt("voldp");
  mr.dp =in.getflt("mrdp");
  ct.dp=in.getflt("ctdp");
  cr.dp=in.getflt("crdp");
  mtr.dpr=mtrL.dpr=in.getflt("mrdp");

#ifdef XYPLANE
  mt.dpv.z=0;
#endif

  // Analysis and energy
  double usys=0;
#pragma omp parallel for reduction (+:usys) schedule (dynamic)
  for (int i=0; i<g.size()-1; i++)
    for (int j=i+1; j<g.size(); j++)
      usys+=pot.energy(cell.p, g[i], g[j]);
  systemenergy sys(usys);   // System energy analysis
  
  histogram lendist(in.getflt("binlen",1.) ,in.getflt("minlen"), in.getflt("maxlen"));             
  aggregation agg(cell, g, in.getflt("aggdef",1.5) );
  distributions dist(in.getflt("binlen",1.),in.getflt("minlen"), in.getflt("maxlen"));

  FAUrdf protrdf(0,0,2.0,cell.len/2.);   // Protein and salt radial distributions
  FAUrdf protrdf11(0,0,.5,cell.len/2.); // Protein and salt radial distributions
  FAUrdf protrdf22(0,0,.5,cell.len/2.); // Protein and salt radial distributions
  FAUrdf protrdf12(0,0,.5,cell.len/2.); // Protein and salt radial distributions

  // Switch parameters
  double sum, volr, tr, rr, clt, clr, randy=-1;
  volr = in.getflt("volr"), tr=in.getflt("tr"), rr=in.getflt("rr"), clt=in.getflt("clt"), clr=in.getflt("clr");
  sum  = volr+tr+rr+clt+clr, volr/=sum, tr/=sum, rr/=sum, clt/=sum, clr/=sum;
  int switcher=-1;

  // Help variables
  int eprint, cnt=0;
  eprint=int(0.001*in.getflt("microsteps")); //loop.macro);

  ioxtc xtc(cell.len);                                // Gromacs xtc output
  for (int i=0; i<g.size(); i++)
    xtc.g.push_back( &g[i] );

  cout << slump.info()
       << cell.info() << pot.info() <<in.info() << endl << endl // Print information to screen
       << "---------- RUN-TIME INFORMATION  -----------" << endl;

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,j,n;
      cnt++;
      randy=slump.random_one();
      if (randy<volr) {                         // Pick a random MC move
        sys+=vol.move(g);                       // Volume move
        lendist.add(cell.len);                  // Update L distribution
      }
      // if (randy<volr+rr && randy >volr)         // Rotate proteins
      // for (n=0; n<g.size(); n++) {            //   Loop over all proteins
      // i = slump.random_one()*g.size();      //   and pick at random.
      // sys+=mr.move(g[i]);                   //   Do the move.
      // }
      if (randy<volr+rr+tr && randy>volr)       // Translate and rotate proteins simultaneously
        for (n=0; n<g.size(); n++) {            //   Loop over all proteins
          i = slump.random_one()*g.size();      //   and pick at random.
          if (slump.random_one()>0.9) {
            mtrL.dpt = cell.len / 2;
            sys+=mtrL.move(g, i);
          }
          else {
            // sys+=mt.move(g[i]);               //  Do the move.
            mtr.dpt = dppercent * cell.len * cell.len;
            sys+=mtr.move(g, i);                //  Do the move
          }
        }

      if (randy<clt+volr+rr+tr && randy>volr+rr+tr)// Cluster translation
        sys+=ct.move(g);                          //   Do the move.

      if (randy<clt+volr+rr+tr+clr && randy>volr+rr+tr+clt) {
        sys+=cr.move(g);                          //  Cluster rotation
      }
      // Debug-code, leave this for a short while, spurius results in connection to clusterrotation
/*        for (int i=0; i<g.size()-1; i++)
          for (int j=i+1; j<g.size(); j++)
            usys+=pot.energy(cell.p, g[i], g[j]);
        sys.update(usys);                           // Update system energy averages
        if (abs(sys.sum-sys.cur)>0.0001)
          cr.print();
        for (int i=0; i<cr.cluster.size();i++)
          cout <<" "<<g[cr.cluster[i]].cm;
        cout <<endl;

      }*/
      // End debug
      dist.add("aveenergy", cell.len, sys.sum);

      if (movie==true && slump.random_one()>.995 && macro>1)
        xtc.save("confout.xtc", cell);   // Save trajectory

      if (slump.random_one()>.95)
        sys.track();

      if(slump.random_one()>.95) {
        for (i=0;i<g.size();i++) 
          g[i].masscenter(cell);                // Recalculate mass centers

        for (int i=0; i<g.size()-1; i++) {
          int nprot2=in.getint("nprot2");
          for (int j=i+1; j<g.size(); j++) {        //   Analyse g(r)...
            protrdf.update(cell, g[i].cm, g[j].cm);
            if ((i<nprot2 && j>=nprot2) || (i>=nprot2 && j<nprot2))
              protrdf12.update(cell, g[i].cm, g[j].cm);
            else {
              if (i<nprot2)
                protrdf11.update(cell, g[i].cm, g[j].cm);
              else
                protrdf22.update(cell, g[i].cm, g[j].cm);
            }
          }
        }
        agg.count();
      }
    } // End of inner loop

    usys=0;
#pragma omp parallel for reduction (+:usys) schedule (dynamic)
    for (int i=0; i<g.size()-1; i++)
      for (int j=i+1; j<g.size(); j++)
        usys+=pot.energy(cell.p, g[i], g[j]);
    sys.update(usys);                           // Update system energy averages

    cout << loop.timing(macro);
    cout << "#   Energy (cur, avg, std)    = "
      << sys.cur << ", " << sys.uavg.avg()<<", "<<sys.uavg.stdev()<<endl
      << "#   Drift                     = " << sys.cur-sys.sum<<endl;

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
    //io.writefile("isobaric.conf", in.print());
    aam.save("confout.aam", cell.p);

    protrdf.write("rdfprot.dat");               // Write g(r)'s
    protrdf11.write("rdfprot11.dat");           // Write g(r)'s
    protrdf12.write("rdfprot12.dat");           // Write g(r)'s
    protrdf22.write("rdfprot22.dat");           // Write g(r)'s
    agg.write("aggregates.dat");
    dist.write("aveenergy.dat");

    pqr.save("confout.pqr", cell.p);
    cell.saveToDisk("confout.dump");            // Save container to disk

  } // End of outer loop

  xtc.close();

  if (in.getboo("penalize")==false)
    vol.printupdatedpenalty("penalty.dat");

  cout << "----------- FINAL INFORMATION -----------" << endl
    << loop.info() << sys.info() << agg.info() << vol.info()
    << mtr.info() << mtrL.info() << ct.info() << cr.info() << pot.info() //<< mt.info() << mr.info()
    << endl
    << "#   Final      side length  = " << cell.len << endl
    << "#   Ideal     <side length> = " << pow( double( g.size() ) / in.getflt("pressure" ),1./3.)<<endl
    << "#   Simulated <side length> = " << vol.len.avg() << " (" << vol.len.stdev() << ")" << endl
    << "#   Ideal     <density>     = " << in.getflt("pressure") << endl
    << "#   Simulated <density>     = " << g.size()*vol.ivol.avg() << " (" << g.size()*vol.ivol.avg() << ")" << endl;

  // Unit testing
  ct.check(test);
  mtr.check(test);
  vol.check(test);
  sys.check(test);
  cout << test.report();
  return test.returnCode();
}

