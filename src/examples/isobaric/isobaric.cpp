/*
 * This program will simulate an arbitrary number of
 * rigid macromolecules in a Debye-Huckel salt solution.
 *
 * \author Bjorn Persson, Anil Kurut, Mikael Lund
 */

#include "faunus/faunus.h"
#include "faunus/energy/coarsegrain.h"
#include "faunus/potentials/pot_debyehuckelP3.h"
#include "faunus/potentials/pot_hsdebyehuckelP3.h"
#include "faunus/potentials/pot_hsHamakerDH.h"

using namespace Faunus;
using namespace std;
 
#ifdef MONOPOLE
  typedef interaction_monopole<pot_debyehuckelP3> Tpot;
#elif HSDH
  typedef interaction<pot_hsdebyehuckelP3> Tpot;
#elif HAMAKER
  typedef interaction<pot_hsHamakerDH> Tpot;
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
  slp.random_seed(9);
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
  iopqr pqr;
  ioaam aam;
  xyfile boxlen("boxlen.dat");            // Box length as function of MC step
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
  translate mtL(nvt, cell, pot, in);      //   Class for macromolecular translation
  transrot mtr(nvt, cell, pot);           //   Class for simultaneous macromolecular translation and rotation
  transrot mtrL(nvt, cell, pot);          //   Class for simultaneous macromolecular translation and rotation
  clustertrans ct(nvt, cell, pot, g);     //   Class for non-rejective cluster translation
  clusterrotate cr(nvt, cell, pot);       //   Class for flux-rejective cluster rotation
  cr.sep=in.getflt("cluster_sep",0.05);   //   Set cluster definition (move to constructor)
  ct.skipEnergyUpdate=in.getboo("skipEnergyUpdate",false);

  isobaricPenalty<Tpot> vol(nvt, cell, pot, in);

  if ( vol.penalize==true)                 // If penalty.dat is present it will be loaded and used as
    vol.loadpenaltyfunction("penalty.dat");// bias. Else a penalty function will only be constructed
  
  isobaricPenalty<Tpot> volL(nvt, cell, pot, in); 
  
  // Markov parameters
  bool movie=in.getboo("movie",false);
  bool tr_rot=in.getboo("tr_rot");
  double dppercent = in.getflt("dppercent", 0.00080);
  mr.dp =in.getflt("mrdp",2);
  ct.dp=in.getflt("ctdp");
  cr.dp=in.getflt("crdp");
  mtr.dpr=mtrL.dpr=in.getflt("mrdp",2);

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
  
  histogram lendist(in.getflt("isobar_binlen",0.05) ,in.getflt("isobar_minlen"), in.getflt("isobar_maxlen"));             
  aggregation agg(cell, g, in.getflt("aggdef",1.5) );
  distributions dist(in.getflt("isobar_binlen",0.05),in.getflt("isobar_minlen"), in.getflt("isobar_maxlen"));

  FAUrdf protrdf(0,0,1,cell.len/2.);   // Protein and salt radial distributions
  FAUrdf protrdf11(0,0,1,cell.len/2.); // Protein and salt radial distributions
  FAUrdf protrdf22(0,0,1,cell.len/2.); // Protein and salt radial distributions
  FAUrdf protrdf12(0,0,1,cell.len/2.); // Protein and salt radial distributions

  // Switch parameters
  double sum, volr, tr, rr, clt, clr, randy=-1;
  volr = in.getflt("volr"), tr=in.getflt("tr"), rr=in.getflt("rr"), clt=in.getflt("clt"), clr=in.getflt("clr");
  sum  = volr+tr+rr+clt+clr, volr/=sum, tr/=sum, rr/=sum, clt/=sum, clr/=sum;
  int switcher=-1;

  ioxtc xtc(cell.len);                                // Gromacs xtc output
  for (int i=0; i<g.size(); i++)
    xtc.g.push_back( &g[i] );

  cout << slp.info() << loop.info()
       << cell.info() << pot.info() << in.info() << endl << endl // Print information to screen
       << "---------- RUN-TIME INFORMATION  -----------" << endl;

  int cnt=0;
  while ( loop.macroCnt() ) {                     //Markov chain 
    while ( loop.microCnt() ) {
      cnt++;
      randy=slp.random_one();
      if (randy<volr) {                         // Pick a random MC move
        if (slp.random_one()>0.9) {
          volL.dp=2;                            // set large Volume displacement 
          sys+=volL.move(g);                    // Large Volume move
          lendist.add(cell.len);                // Update L distribution
        }
        else {
          sys+=vol.move(g);                       // Volume move
          lendist.add(cell.len);                  // Update L distribution
        }
        if (slp.random_one()>0.9)
          boxlen.add(loop.count(), cell.len);
      }


      if (tr_rot==true) {                         //   Combined Translate and Rotation
        if (randy<volr+rr+tr && randy>volr)       //   Translate and rotate proteins simultaneously
          for (int n=0; n<g.size(); n++) {        //   Loop over all proteins
            int i = slp.random_one()*g.size();    //   and pick at random.
            if (slp.random_one()>0.9) {           //   10 % is long translation
              mtrL.dpt = cell.len * 2;
              sys+=mtrL.move(g, i);
            }
            else {
              mtr.dpt = dppercent * cell.len * cell.len;
              sys+=mtr.move(g, i);                //   Do the move.
            }
          }
      }
      else {                                      //   Independent translate and rotation
        if (randy<volr+rr && randy >volr) {       //   Rotate proteins
          for (int n=0; n<g.size(); n++) {        //   Loop over all proteins
            int i = slp.random_one()*g.size();    //   and pick at random.
            sys+=mr.move(g[i]);                   //   Do the move.
          }
        }

        if (randy<volr+rr+tr && randy >volr+rr) { //   Translate proteins
          for (int n=0; n<g.size(); n++) {        //   Loop over all proteins
            int i = slp.random_one()*g.size();    //   and pick at random.
            if (slp.random_one()>0.0) {           //   10 % is long translation
              mtL.dp = cell.len / 2;
              sys+=mtL.move(g.at(i));             //   Do the move.
            }
            else {
              mt.dp = dppercent * cell.len * cell.len;
              sys+=mt.move(g[i]);                 //   Do the move.
            }
          }
        }
      }

      if (randy<clt+volr+rr+tr && randy>volr+rr+tr)//  Cluster translation
        sys+=ct.move(g);                           //  Do the move.

      if (randy<clt+volr+rr+tr+clr && randy>volr+rr+tr+clt)
        sys+=cr.move(g);                          //  Cluster rotation

      dist.add("aveenergy", cell.len, sys.sum);

      if (movie==true && slp.random_one()>.995)
        xtc.save("confout.xtc", cell);   // Save trajectory

      if (slp.random_one()>.95)
        sys.track();

      // Analysze g(r)'s
      if(slp.random_one()>0.9) {
        int n=g.size();
        int n1=in.getint("nprot1",0);
        for (int i=0; i<n; i++) 
          g[i].masscenter(cell);
        for (int i=0; i<n-1; i++) {
          for (int j=i+1; j<n; j++) {
            double r=sqrt( pot.pair.sqdist(g[i].cm,g[j].cm) );
            protrdf.add(r);
            if ((i<n1 && j>=n1) || (i>=n1 && j<n1))
              protrdf12.add(r);
            else {
              if (i<n1)
                protrdf11.add(r);
              else
                protrdf22.add(r);
            }
          }
        }
      //agg.count();
      }
    } // End of inner loop

    usys=0;
#pragma omp parallel for reduction (+:usys) schedule (dynamic)
    for (int i=0; i<g.size()-1; i++)
      for (int j=i+1; j<g.size(); j++)
        usys+=pot.energy(cell.p, g[i], g[j]);
    sys.update(usys);                           // Update system energy averages
    sys.write();

    cout << loop.timing()
         << "#   Energy drift (kT) = " << sys.cur-sys.sum << endl;
    
    cell.check_vector();                        // Check sanity of particle vector
    protrdf.write("rdfprot.dat");               // Write g(r)'s
    if (vol.penalize==true) {
      vol.printpenalty("oldpenalty.dat");       // Print penalty function (used as bias)
      vol.printupdatedpenalty("penalty.dat");   // Print updated penalty function (not used as bias)
    }

    lendist.write("length-distribution.dat");   // Write volume distribution
    if (in.getboo("penaltyupdate")==true) {     // Dont use for production
      vol.updatepenalty();
      cout << "# Penalty function updated"<<endl;
    }

    in.updateval("boxlen", stringify(cell.len));
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
  boxlen.close();

  if (vol.penalize==true)
    vol.printupdatedpenalty("penalty.dat");

  cout << "----------- FINAL INFORMATION -----------" << endl
       << loop.info() << sys.info() << agg.info()
       << vol.info() << volL.info() << mr.info()
       << mtr.info() << mtrL.info() << mt.info() << mtL.info()
       << ct.info() << cr.info();

  // Unit testing
  ct.check(test);
  mtr.check(test);
  vol.check(test);
  sys.check(test);
  cout << test.report();
  return test.returnCode();
}

