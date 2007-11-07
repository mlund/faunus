/*
 * Monte Carlo simulation program to calculate
 * protonation status of a flexible chain in a
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
#include "interact.h"
#include "peptide.h"
#include "inputfile.h"
#include "montecarlo.h"
#include "vpython.h"
#include "povray.h"
#include "titrate.h"
#include "average.h"

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
  double pH       = par.getDbl("pH", 7);
  double springk  = par.getDbl("springconst",0.5);
  double springeq = par.getDbl("springeqdist",4.9);

  string protein1 = par.getStr("protein1");

  int monomer_dp  = par.getInt("monomer_dp");
  int ion_dp      = par.getInt("ion_dp"); //ion displacement
  int nion1       = par.getInt("nion1");  //number of ion 1
  int nion2       = par.getInt("nion2");  //number of ion 2
  int nion3       = par.getInt("nion3");  //number 
  int chion1      = par.getInt("chion1"); //ionic charges
  int chion2      = par.getInt("chion2");
  int chion3      = par.getInt("chion3");

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
  interact pot(phys.beta_ecf);
  collision col;
  vector<group> g(2);
  average Q, Q2, mu;
  average end2end;

  pot.k   = springk;
  pot.r0  = springeq;

  /**********************************
    Load proteins, salt and polymers
  ***********************************/
  g[PROT] = s.append( pep.loadAAModel(protein1, peptide::CHAIN) );
  g[PROT].name = "protein 1";
  g[PROT].chain = true;
  s.accept(g[PROT], space::ALL);
  
  g[SALT].name="mobile ions";
  g[SALT]+=s.insert_salt( nion1, chion1, 2.0, cell);
  g[SALT]+=s.insert_salt( nion2, chion2, 2.0, cell);
  g[SALT]+=s.insert_salt( nion3, chion3, 2.0, cell);

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
       << "# Dielectric Const.   = " << phys.e_r << endl
       << "# Bjerrum length      = " << pot.lB << endl
       << "# Protein 1           = " << protein1 << endl
       << "# Ion 1 (#  chg)      = " << nion1 << " " << chion1 << endl
       << "# Ion 2 (#  chg)      = " << nion2 << " " << chion2 << endl
       << "# Ion 3 (#  chg)      = " << nion3 << " " << chion3 << endl
       << "# System charge       = " << s.charge() << endl
       << "# Initial energy (kT) = " << utot << endl;
  
  tit.info();  
  cout << "# --- GROUPS ----------------------------------------\n";
  cout << g[PROT] << g[SALT] << endl;

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
	      uold = pot.energy(s.p, n);
	      unew = pot.energy(s.trial,n);
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

      //Rattle chains
      for (int i=0; i<g[PROT].size(); i++) {
	rc=false;
	rejectcause=montecarlo::HC;
	int n=g[PROT].random();
	s.displace(n, monomer_dp);
	if (cell.cellCollision(s.trial[n])==false)
	  if (col.overlap(s.trial, n)==false) {
	    uold=pot.chain(s.p, g[PROT], n);
	    unew=pot.chain(s.trial, g[PROT], n);
	    if (s.p[n].charge!=0) {
	      uold+=pot.energy(s.p, n);
	      unew+=pot.energy(s.trial, n);
	    };
	    du = unew-uold;
	    if (pot.metropolis(du))
	      rc=true;
	    else
	      rejectcause=montecarlo::ENERGY;
	  };
	if (rc==true) {
	  s.p[n]=s.trial[n];
	  mc.accept(montecarlo::MONOMER);
	  utot += du;
	} else {
	  s.trial[n]=s.p[n];
	  mc.reject(montecarlo::MONOMER, rejectcause);
	};
	end2end += s.p[g[PROT].beg].dist( s.p[g[PROT].end] );
      };
      
      //Charge regulation
      if (slump.random_one()>0.7) {
	titrate::action t;
	for (int cnt=0; cnt<tit.sites.size(); cnt++) {
	  t=tit.exchange(s.trial);
	  uold=pot.energy( s.p, t.site ) + pot.energy( s.p, t.proton )
	    - pot.energy(s.p[t.site], s.p[t.proton])*pot.lB;
	  unew=pot.energy( s.trial, t.site ) + pot.energy( s.trial, t.proton )
	    - pot.energy(s.trial[t.site], s.trial[t.proton])*pot.lB;
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
        Q +=q;
        Q2+=q*q;
        tit.samplesites(s.p);
      };
    }; //end of micro loop
    
    if (s.safetytest_vector()==false) return 1;

    // System energy check:
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy(s.p);
    cout << "# System energy (kT)       = " << usys << endl
         //<< "# System energy drift (kT) = " << utot - usys << endl
         << "# Protein charge           = " << s.charge(g[PROT]) << endl
         << "# Total system charge      = " << s.charge() << endl
         << "# Average protein charge   = " << Q << endl
         << "# End-2-end distance       = " << end2end.avg() << endl
         << "# Protein capacitance      = " << Q2.avg()-pow(Q.avg(),2) << endl 
         << "# (Alternative Q calc)     = " << tit.sumsites() << endl
         << endl;

    s.save( ".coord" + jobid ); //save coordinates

    //Graphical output (Povray)
    povray pov;
    pov.sphere(cell_r);
    pov.zaxis(cell_r);
    pov.add( s.p, g[PROT] );
    pov.add( s.p, g[SALT] );
    pov.save("snapshot.pov" + jobid);

  }; //end of macro step 
    
  /**************************************
    Print results and terminate program
   **************************************/
  mc.showStatistics();
  tit.showsites(s.p);

  cout << "pH=" << pH
       << " <Q>="<< Q.avg()
       << " C=" << Q2.avg()-pow(Q.avg(),2)
       << " <R>=" << end2end.avg() << endl;
  
  return 0;
};
