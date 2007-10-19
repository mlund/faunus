/*
 * Monte Carlo simulation program to calculate
 * protonation status of a single protein in a
 * salt solution.
 *
 * (C) M. Lund, 2005.
 *
 * Ref: Biochemistry, 2005(44), 5722-5727.
 *      "On the charge regulation of proteins"
 */

#include <iostream>
#include "point.h"
#include "group.h"
#include "space.h"
#include "slump.h"
#include "reactionfield.h"
#include "peptide.h"
#include "inputfile.h"
#include "montecarlo.h"
#include "vpython.h"
#include "povray.h"
#include "titrate.h"
#include "average.h"
#include "histogram.h"

using namespace std;

int main(int argc, char* argv[] ) {

  /**********************************
    Handle arguments and input file
   **********************************/
  string inputName = "amino.inp";
  if (argc==2)
    inputName = argv[1];
  inputfile par( inputName );
  int randomseed  = par.getInt("randomseed",13);
  int macro       = par.getInt("macro");
  int micro       = par.getInt("micro");
  string jobid    = par.getStr("jobid"); 
  double cell_r   = par.getDbl("cell_r");
  double temp     = par.getDbl("temp"); 
  double dielec   = par.getDbl("dielec",78);
  double dielec_i = par.getDbl("dielec_i",20);
  double r_dielec = par.getDbl("r_dielec",10);
  double pH       = par.getDbl("pH", 7);
  double springk  = par.getDbl("springconst",0.5); 

  string protein1 = par.getStr("protein1");

  int ion_dp      = par.getInt("ion_dp"); //ion displacement
  int nion1       = par.getInt("nion1");  //number of ion 1
  int nion2       = par.getInt("nion2");  //number of ion 2
  int nion3       = par.getInt("nion3");  //number 
  int chion1      = par.getInt("chion1"); //ionic charges
  int chion2      = par.getInt("chion2");
  int chion3      = par.getInt("chion3");

  bool saveaam    = par.getBool("saveAAM", false);
  bool movesites  = par.getBool("movesites", false);

  /**************************************
    Class constructors and declarations
   **************************************/
  enum groupid {PROT=0,SALT}; //convenient
  montecarlo mc( macro, micro );
  physconst phys(temp, dielec);
  geometry geo;
  cell cell(cell_r);
  slump slump;
  space s;
  peptide pep;
  interact_rf pot(phys.beta_ecf, r_dielec, dielec, dielec_i);
  collision col;
  vector<group> g(2);
  average Q, Q2, mu, muex;
  avgmatrix phi(-20,20,-20,20,0.5); //grid of electroc potential

  pot.k   = springk;
  pot.r0  = 0.0;

  //Distribution functions
  vector<histogram::data> h(4);
  histogram hist(0.5, cell_r);
  hist.init(h);
  h[0].name = "Q" ; h[0].type=histogram::AVERAGE;
  h[1].name = "Q2" ; h[1].type=histogram::AVERAGE;
  h[2].name = "ex1" ; h[2].type=histogram::AVERAGE;
  h[3].name = "ex2" ; h[3].type=histogram::AVERAGE;

  /**********************************
    Load proteins, salt and polymers
  ***********************************/
  g[PROT] = s.append( pep.loadAAModel(protein1) );
  g[PROT].name = "protein 1";
  s.addmasscenter(g[PROT]);
  s.adddipole(g[PROT]);
  s.move( g[PROT], -s.center_of_mass(g[PROT]), space::AUTOACCEPT);
  g[PROT].radius=s.radius( g[PROT], s.p[g[PROT].cm] );

  g[SALT].name="mobile ions";
  g[SALT]+=s.insert_salt( nion1, chion1, 2.0, cell);
  g[SALT]+=s.insert_salt( nion2, chion2, 2.0, cell);

  if (s.safetytest_vector()==false) return 1;

  s.load(".coord" + jobid); // load saved coordinates and charges.

  //PREPARE TITRATION
  titrate tit(s.p, g[SALT], pH );
  s.load(".coord" + jobid); //load coords again (to get old charges)
  s.trial=s.p;
  if (s.charge()!=0) {
    cout << "Aborting! System not electroneutral. (Q="<<s.charge()<<")\n";
    return 1;
  };
  
  //INITIAL ENERGY
  double utot = pot.energy(s.p);

  //PRINT INFORMATION
  cout << "\n# --- SYSTEM ----------------------------------------\n"
       << "# Macro loops         = " << mc.macroSteps << endl
       << "# Micro loops         = " << mc.microSteps << endl
       << "# Total nr. of loops  = " << mc.macroSteps*mc.microSteps << endl
       << "# Cell radius         = " << cell_r << endl
       << "# Temperature (K)     = " << phys.T << endl
       << "# Dielectric (out,in) = " << phys.e_r << " " << dielec_i << endl
       << "# Dielec-sphere-radius= " << r_dielec << endl
       << "# Bjerrum length      = " << pot.lB << endl
       << "# Protein 1           = " << protein1 << endl
       << "# Ion 1 (#  chg)      = " << nion1 << " " << chion1 << endl
       << "# Ion 2 (#  chg)      = " << nion2 << " " << chion2 << endl
       << "# Ion 3 (#  chg)      = " << nion3 << " " << chion3 << endl
       << "# System charge       = " << s.charge() << endl
       << "# Initial energy (kT) = " << utot << endl;
  if (movesites==true)
    cout<<"# Spring constant     = " << pot.k << endl;
  tit.info();
  
  cout << "# --- GROUPS ----------------------------------------\n";
  cout << g[PROT]
       << "#   Mass center   = " << s.p[g[PROT].cm] << endl
       << "#   Avg. radius   = " << s.radius(g[PROT], s.p[g[PROT].cm], space::AVERAGE) << endl
       << "#   Interior size = " << s.radius(g[PROT], s.p[g[PROT].cm], space::INTERIOR) << endl
       << g[SALT] << endl;

  /************************
    Monte Carlo Simulation
   ************************/
  for (int macroCnt=1; macroCnt<=mc.macroSteps; macroCnt++) {
    for (int i=1; i<=mc.microSteps; i++) {
      bool rc,coll;
      montecarlo::rejectcause rejectcause;
      double u=0,uold=0,unew=0,du=0;

      //Move salt
      for (int j=0; j<g[SALT].size(); j++) {
	rc=false;
	rejectcause=montecarlo::HC;
	int n=g[SALT].random();  //pick random particle
	s.displace(n, ion_dp);   //displace it...
	if (cell.cellCollision(s.trial[n])==false)
	  if (col.overlap(s.trial, g[PROT], n)==false)
	    if (col.overlap(s.trial, g[SALT], n)==false) {
	      uold = pot.energy_rf(s.p, n); // kT
	      unew = pot.energy_rf(s.trial,n); // kT
	      du = unew - uold;
	      if (pot.metropolis(du)==true)
		rc=true;
	      else
		rejectcause=montecarlo::ENERGY;
	    };
	if (rc==true) {
	  s.p[n]=s.trial[n];
	  mc.accept(montecarlo::ION);
	  utot += du;
	} else {
	  s.trial[n]=s.p[n];
	  mc.reject(montecarlo::ION, rejectcause);
	};
      };

      //Salt (1:1) excess chemical potential
      //as a function of distance to the center
      if (slump.random_one()>0.95) {
	double u_an, u_cat;
        particle anion, cation;
        anion.charge=-1;
        cation.charge=+1;
        anion.radius=cation.radius=2.0;
        for (int k=0; k<100; k++) {
          cell.randomPos(cation);
          cell.randomPos(anion);
          u=u_an=u_cat=1.0e20;
          if (col.overlap(s.p, anion)==false)
            if (col.overlap(s.p, cation)==false) {
              u_cat= pot.energy_rf(s.p, cation);
              u_an = pot.energy_rf(s.p, anion);
	      u=u_cat+u_an;
	    };
	  hist.add( h[2], cation.len(), exp(-u_cat));
	  hist.add( h[3], anion.len(), exp(-u_an));
          muex += exp(-u);
        };
      };

      //Calculate average electric potential in a grid (the ratio, that is).
      if (slump.random_one()>0.8) {
	particle test,origo;
	test.radius=1.0;
	for (double x=phi.xbeg; x<phi.xend; x+=phi.steps) {
	  for (double y=phi.ybeg; y<phi.yend; y+=phi.steps) {
	    test.x = x;
	    test.y = y;
	    if (col.overlap(s.p, g[SALT], test)==false)
	      phi.add(x,y,pot.potential_rf( s.p, test ) );
	  };
	};
      };
      
      //Charge regulation
      if (slump.random_one()>0.7) {
	titrate::action t;
	for (int cnt=0; cnt<tit.sites.size(); cnt++) {
	  t=tit.exchange(s.trial);
	  uold=pot.energy_rf( s.p, t.site ) + pot.energy_rf( s.p, t.proton )
	    - pot.energy_rf(s.p[t.site], s.p[t.proton])*pot.lB*pot.eo;
	  unew=pot.energy_rf( s.trial, t.site ) + pot.energy_rf( s.trial, t.proton )
	    - pot.energy_rf(s.trial[t.site], s.trial[t.proton])*pot.lB*pot.eo;
	  du=unew-uold;
	  rc=pot.metropolis( tit.energy(s.trial, du, t) );
	  if (rc==true) {
	    s.p[t.site]=s.trial[t.site];
	    s.p[t.proton]=s.trial[t.proton];
	    mc.accept(montecarlo::TITRATE);
	    utot += du;
	  } else {
	    tit.exchange(s.trial, t);
	    mc.reject(montecarlo::TITRATE, montecarlo::ENERGY);
	  };
	};
        double q=s.charge(g[PROT]);
        s.calcdipole(g[PROT]);
        Q +=q;
        Q2+=q*q;
        mu+=g[PROT].dip;
        tit.samplesites(s.p);
      };

      //Displace titrateable sites
      if (movesites==true && slump.random_one()>0.7) {
        for (int cnt=0; cnt<tit.sites.size(); cnt++) {
	  rc=false;
	  rejectcause=montecarlo::HC;
          int n=tit.sites[cnt];
          s.displace(n, 2.);
	  if (col.overlap(s.trial, g[SALT], n)==false) {
	    uold = pot.energy_rf(s.p, n) + pot.quadratic(s.p[n], tit.eqpos[cnt]);
	    unew = pot.energy_rf(s.trial,n) + pot.quadratic(s.trial[n], tit.eqpos[cnt]);
	    du = unew - uold;
	    if (pot.metropolis(du)==true)
	      rc=true;
	    else
	      rejectcause=montecarlo::ENERGY;
	  };
	  if (rc==true) {
	    s.p[n]=s.trial[n];
	    mc.accept(montecarlo::SITE);
	    utot += du;
	  } else {
	    s.trial[n]=s.p[n];
	    mc.reject(montecarlo::SITE, rejectcause);
	  };
        }
      }

      //Charge histogram of cell
      if (slump.random_one()<0.8) {
        double q;
        point origo;
        for (double r=0; r<cell_r; r+=0.5) {
          q=s.charge(origo, r);
          hist.add(h[0], r, q);
          hist.add(h[1], r, q*q);
        };
      };
      
    }; //end of micro loop
    
    if (s.safetytest_vector()==false) return 1;

    // System energy check:
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy(s.p);
    cout << "# System energy (kT)       = " << usys << endl
         << "# System energy drift (kT) = " << utot - usys << endl
         << "# Protein charge           = " << s.charge(g[PROT]) << endl
         << "# Total system charge      = " << s.charge() << endl
         << "# Average protein charge   = " << Q << endl
         << "# Protein dipole moment    = " << mu << endl
         << "# Protein capacitance      = " << Q2.avg()-pow(Q.avg(),2) << endl 
         << "# (Alternative Q calc)     = " << tit.sumsites() << endl
         << endl;

    s.save( ".coord" + jobid ); //save coordinates
    phi.save("potgrid" + jobid); //save potential grid

    //Graphical output (Povray)
    povray pov;
    pov.sphere(cell_r);
    pov.zaxis(cell_r);
    pov.add( s.p, g[PROT] );
    pov.add( s.p, g[SALT] );
    pov.save("snapshot.pov" + jobid);

    vrml vrml;
    vrml.add( s.p, g[PROT] );
    vrml.add( s.p, g[SALT] );
    vrml.save("vrml.py");
    
  }; //end of macro step 
    
  /**************************************
    Print results and terminate program
   **************************************/
  hist.show(h,0,cell_r);
  mc.showStatistics();
  tit.showsites(s.p);
  
  cout << "pH=" << pH
       << " <Q>="<< Q.avg()
       << " C=" << Q2.avg()-pow(Q.avg(),2)
       << " <mu>=" << mu.avg() << endl;

  //Smear charges and save AAM file to disk
  tit.applycharges(s.p);
  pep.saveAAModel("saved.aam", s.p, g[PROT]);

  return 0;

};
