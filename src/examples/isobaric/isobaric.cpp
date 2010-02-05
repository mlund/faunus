/*
 * This program will simulate an arbitrary number of
 * protein molecules in a dielectric solvent with
 * explicit mobile ions.
 *
 * \author Bjorn Persson and Anil Kurut
 */
#include "faunus/faunus.h"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

using namespace Faunus;
using namespace std;

#ifdef MONOPOLE
  typedef interaction_monopole<pot_debyehuckelP3> Tpot;
#elif defined(FASTDH)
  typedef interaction_vector<pot_debyehuckelP3Fast> Tpot;
#else
  typedef interaction<pot_debyehuckelP3> Tpot;
#endif

int main() {
  cout << faunus_splash();
  cout << "---------- INITIAL PARAMETERS -----------" << endl;
  slump slump;                            // A random number generator
  physconst phys;
  inputfile in("isobaric.conf");          // Read input file
  phys.e_r=in.getflt("e_r");
  phys.lB_TO_T(in.getflt("bjerrum"));
  mcloop loop(in);                        // Keep track of time and MC loop
#ifndef XYPLANE
  box cell(in);                           // We want a cubic cell
#else
  xyplane cell(in);                       // ...or a plane!
#endif
  canonical nvt;                          // Use the canonical ensemble
#ifdef MONOPOLE
  Tpot pot(in,cell);                      // Far-away monopole approximation
#else
  Tpot pot(in);                           // Fast, approximate Debye-Huckel potential
#endif
  vector<macromolecule> g;                // PROTEIN groups
  io io;
  iopqr pqr;
  ioaam aam;
  if (in.getboo("lattice")==true)         //   Protein input file format is AAM
    aam.loadlattice(cell, in, g);                //   Load and insert proteins
  else                                    //   Center first protein (will be frozen)
    aam.load(cell, in, g);
  if (aam.load(cell,"confout.aam")) {
    for (int i=0;i<g.size();i++) 
      g[i].masscenter(cell);              // Load old config (if present)
  }

  // Markovsteps
  macrorot mr(nvt, cell, pot);            //   Class for macromolecule rotation
  translate mt(nvt, cell, pot);           //   Class for macromolecular translation
  clustertrans ct(nvt, cell, pot, g);     //   Class for non-rejective cluster translation
  isobaric< Tpot > vol(
      nvt, cell, pot,
      in.getflt("pressure"),              // Set external pressure (in kT)
      in.getflt("penalty"),               // Set penalty           (in kT)
      in.getflt("scalepen"),              // Scale factor
      int(in.getint("maxlen")),           // Set maximum box length
      int(in.getint("minlen")),           // Set minimum box length
      int(in.getint("binlen")));          // Set bin length for penalty function  

  if ( in.getboo("penalize")==true)         // If penalty.dat is present it will be loaded and used as
    vol.loadpenaltyfunction("penalty.dat"); // bias. Else a penalty function will only be constructed
                                            // if the penalty != 0.
  // Move details
  vol.dp=in.getflt("voldp");
  mt.dp =in.getflt("mtdp");
  mr.dp =in.getflt("mrdp");
  ct.dp=in.getflt("ctdp");
  double cltr,rr,tr,volr,sum;                // If anyone doesnt fancy this way of controlling the
  cltr=in.getflt("cltr",1);                  // Markov-chain, write a new program, dont remove this. /BearniX
  rr  =in.getflt("rr",1);
  tr  =in.getflt("tr",1);
  volr=in.getflt("volr",1);
  sum=1./(cltr+rr+tr+volr);
  cltr*=sum,rr*=sum,tr*=sum,volr*=sum;
#ifdef XYPLANE
  mt.dpv.z=0;
#endif

  // Analysis and energy
  double usys=0;
  for (int i=0; i<g.size()-1; i++)
    for (int j=i+1; j<g.size(); j++)
      usys+=pot.energy(cell.p, g[i], g[j]);

  systemenergy sys(usys);   // System energy analysis
  histogram lendist(in.getflt("binlen",1.) ,in.getflt("minlen"), in.getflt("maxlen"));             
  aggregation agg(cell, g, 1.5);

  FAUrdf protrdf(0,0,.5,cell.len/2.);   // Protein and salt radial distributions
  FAUrdf protrdf11(0,0,.5,cell.len/2.); // Protein and salt radial distributions
  FAUrdf protrdf22(0,0,.5,cell.len/2.); // Protein and salt radial distributions
  FAUrdf protrdf12(0,0,.5,cell.len/2.); // Protein and salt radial distributions

  ioxtc xtc(cell.len);                                // Gromacs xtc output (if installed)
  for (int i=0; i<g.size(); i++)
    xtc.g.push_back(&g[i]);

  cout << cell.info() << pot.info() <<in.info() << endl
       << "#  Temperature = " << phys.T << " K" << endl << endl
       << "#---------- RUN-TIME INFORMATION  -----------" << endl;

  for (int macro=1; macro<=loop.macro; macro++) {//Markov chain 
    for (int micro=1; micro<=loop.micro; micro++) {
      short i,j,n;
      double randy=slump.random_one();
      if(randy<volr)
          sys+=vol.move(g);
      if(volr<randy && randy<volr+rr)
        for (int n=0; n<g.size(); n++) 
          sys+=mr.move(g[slump.random_one()*g.size()]);
      if(volr+rr<randy && randy<volr+rr+tr)
        for (int n=0; n<g.size(); n++) 
          sys+=mt.move(g[slump.random_one()*g.size()]);
      if(volr+rr+tr<randy) {
          sys+=ct.move(g);
          lendist.add(cell.len);
      }
      if (randy>.99 && in.getboo("movie", false)==true)
        xtc.save("coord.xtc", cell);   // Save trajectory
      if (randy>.99) sys.track();
      if (randy>.9) {
        for (i=0;i<g.size();i++) 
          g[i].masscenter(cell);                  // Recalculate mass centers
        for (i=0; i<g.size()-1; i++)
          for (j=i+1; j<g.size(); j++) {          //   Analyse g(r)...
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
      }
    } // End of inner loop

    usys=0;
#pragma omp parallel for reduction (+:usys) schedule (dynamic)
    for (int i=0; i<g.size()-1; i++) {
      for (int j=i+1; j<g.size(); j++)
        usys+=pot.energy(cell.p, g[i], g[j]);
    }
    sys.update(usys);                           // Update system energy averages

    cout << loop.timing(macro)
         << "#   Energy (cur, avg, std)    = "<<sys.cur<<", "<<sys.uavg.avg()<<", "<<sys.uavg.stdev()<<endl
         << "#   Drift                     = "<<sys.cur-sys.sum<<endl;

    cell.check_vector();                        // Check sanity of particle vector
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
    in.updateval("boxlen", cell.len);
    //io.writefile("isobaric.conf", in.print());
    aam.save("confout.aam", cell.p);

    protrdf.write("rdfprot.dat");               // Write g(r)'s
    protrdf11.write("rdfprot11.dat");           // Write g(r)'s
    protrdf12.write("rdfprot12.dat");           // Write g(r)'s
    protrdf22.write("rdfprot22.dat");           // Write g(r)'s
    agg.write("aggregates.dat");

    aam.save("confout.aam", cell.p);            // Save config. for next run
    pqr.save("confout.pqr", cell.p);

  } // End of outer loop

  if (in.getboo("penalize")==false)
    vol.printupdatedpenalty("penalty.dat");

  cout << "----------- FINAL INFORMATION -----------" << endl
       << loop.info() << sys.info() << agg.info() << vol.info()
       << mr.info() << mt.info() << ct.info() << endl
       << "#   Final      side length  = " << cell.len << endl
       << "#   Ideal     <side length> = " << pow( double( g.size() )/ in.getflt("pressure" ),1./3.)<<endl
       << "#   Simulated <side length> = " << vol.len.avg()<<" ("<<vol.len.stdev()<<")"<<endl
       << "#   Ideal     <density>     = " << in.getflt("pressure")<<endl
       << "#   Simulated <density>     = " << g.size()*vol.ivol.avg()<<" ("<<g.size()*vol.ivol.avg()<<")"<<endl;
  xtc.close();
}

