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
#include "interact.h"
#include "peptide.h"
#include "inputfile.h"
#include "montecarlo.h"
#include "vpython.h"
#include "povray.h"
#include "titrate.h"
#include "average.h"
#include "histogram.h"
#include "imdwrap.h"
#include "pdbclass.h"

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
  
  string protein1 = par.getStr("protein1");

  double ion_dp   = par.getDbl("ion_dp"); //ion displacement
  int nion1       = par.getInt("nion1");  //number of ion 1
  int nion2       = par.getInt("nion2");  //number of ion 2
  int nion3       = par.getInt("nion3");  //number 
  int chion1      = par.getInt("chion1"); //ionic charges
  int chion2      = par.getInt("chion2");
  int chion3      = par.getInt("chion3");
  int watchsite   = par.getInt("watchsite",-1); //site on which to sample excess chem. pot.

  bool adjust_dp  = par.getBool("adjust_dp", false);
  bool saveaam    = par.getBool("saveAAM", false);
  bool movesites  = par.getBool("movesites", false);
  bool movecharges= par.getBool("movecharges",false);

  /**************************************
    Class constructors and declarations
   **************************************/
  enum groupid {PROT=0,SALT}; //convenient
  montecarlo mc( macro, micro );
  physconst phys(temp, dielec);
  pdbclass pdb;
  geometry geo;
  cell cell(cell_r);
  slump slump;
  space s;
  peptide pep;
  interact pot(phys.beta_ecf);
  collision col;
  vector<group> g(2);
  average Q, Q2, mu, mu2, R, R2;
  average widom_A, widom_Uprot, widom_Ubulk;
  average u_self;

  pot.k   = springk;
  pot.r0  = 0.0;

  //Distribution functions
  vector<histogram::data> h(2);
  histogram hist(0.5, cell_r);
  hist.init(h);
  h[0].name = "Q" ; h[0].type=histogram::AVERAGE;
  h[1].name = "Q2" ; h[1].type=histogram::AVERAGE;

  //Dipole moment distribution
  vector<histogram::data> muhis(2);
  histogram hist2(1, 200);
  hist2.init(muhis);
  muhis[0].name = "mu" ; muhis[0].type=histogram::PROBABILITY;
  muhis[1].name = "Q+30"  ; muhis[1].type=histogram::PROBABILITY;

  /**********************************
    Load proteins, salt and polymers
  ***********************************/
  g[PROT] = s.append( pep.loadAAModel(protein1) );
  g[PROT].name = "protein 1";
  s.addmasscenter(g[PROT]);
  cout << "# CM = " << s.p[g[PROT].cm].x << " "
    << s.p[g[PROT].cm].y << " "
    << s.p[g[PROT].cm].z << endl;
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
       << "# Dielectric Const.   = " << phys.e_r << endl
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
       << g[SALT] << endl;

  //VMD support
  pdb.load_particles(s.p);
  pdb.save("test.pqr");
  
  //imdwrap imd(s.p.size());
  //imd.wait_for_connection();

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
	  if (t.site==watchsite) //follow titration process of a specific site...
	    if (t.action==titrate::PROTONATED)
	      widom_A += exp(-du);
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
        mu2+=pow(g[PROT].dip,2);
        hist2.add( muhis[0], g[PROT].dip );
	hist2.add( muhis[1], q+30.);
        tit.samplesites(s.p);
      };

      //Follow titration process - energy contribution.
      if (watchsite!=-1 && slump.random_one()>0.9) {
	double u=0;
	if (movesites==true) {
	  int i;
	  for (i=0; i<tit.sites.size(); i++) //find particle in titration vector
	    if (tit.sites[i]==watchsite) break;
	  u=pot.quadratic(s.p[watchsite], tit.eqpos[i]); //calc. spring energy
	}
	u += pot.energy(s.p); //total energy
	if (s.p[watchsite].charge>0)
	  widom_Uprot += u; //proton is on protein...
	else
	  widom_Ubulk += u; //no, proton is in bulk.
      };

      //Fluctuate first two charges in macromolecule (dipolar)
      if (movecharges==true && slump.random_one()>0.7) {
	rc=false;
	double z=slump.random_one()*5;
	int i=g[PROT].beg, j=i+1;
	s.trial[i].charge = +z;
	s.trial[j].charge = -z;
	uold = pot.energy(s.p, g[PROT]) + pot.energy(s.p[i], s.p[j] );
	unew = pot.energy(s.trial, g[PROT]) + pot.energy(s.trial[i], s.trial[j]);
	du = unew - uold;
	if (pot.metropolis(du)==true)
	  rc=true;
	else
	  rejectcause=montecarlo::ENERGY;
	if (rc==true) {
	  s.p[i]=s.trial[i];
	  s.p[j]=s.trial[j];
	  mc.accept(montecarlo::SITE);
	  utot += du;
	} else {
	  s.trial[i]=s.p[i];
	  s.trial[j]=s.p[j];
	  mc.reject(montecarlo::SITE, rejectcause);
	};
        s.calcdipole(g[PROT]);
        mu+=g[PROT].dip;
        mu2+=pow(g[PROT].dip,2);
        hist2.add( muhis[0], g[PROT].dip );
      }

      //Displace titrateable sites
      if (movesites==true && slump.random_one()>0.7) {
	double x2;
        for (int cnt=0; cnt<tit.sites.size(); cnt++) {
	  rc=false;
	  rejectcause=montecarlo::HC;
          int n=tit.sites[cnt];
          s.displace(n, 2.);
	  if (col.overlap(s.trial, g[SALT], n)==false)
	    if (col.overlap(s.trial, tit.sites, 2.)==false) {
	      uold = pot.energy(s.p, n) + pot.quadratic(s.p[n], tit.eqpos[cnt]);
	      unew = pot.energy(s.trial,n) + pot.quadratic(s.trial[n], tit.eqpos[cnt]);
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
	  //sample site fluctuations
	  x2=s.p[n].sqdist(tit.eqpos[cnt]);
	  R2+=x2;
	  if (x2!=0)
	    x2=sqrt(x2);
	  R+=x2;
        }
        u_self+=pot.internal(s.p, tit.sites); //sample protein self-energy

        s.calcdipole(g[PROT]);
        mu+=g[PROT].dip;
        mu2+=pow(g[PROT].dip,2);
        hist2.add( muhis[0], g[PROT].dip );
        tit.samplesites(s.p);
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

      //Adjust displacement parameter
      if (adjust_dp==true && slump.random_one()>0.9) {
        mc.adjust_dp(montecarlo::ION, ion_dp);
      };

      //imd.send_particles(s.p);

    }; //end of micro loop
    
    if (s.safetytest_vector()==false) return 1;

    // System energy check:
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy(s.p);
    cout << "# System energy (kT)       = " << usys << endl
         << "# System energy drift (kT) = " << utot - usys << endl
         << "# Ion displacement (A)     = " << ion_dp << endl
         << "# Protein charge           = " << s.charge(g[PROT]) << endl
         << "# Total system charge      = " << s.charge() << endl
         << "# Protein:"                    << endl
         << "#   Charge <Z> <Z2>        = " << Q << " " << Q2.avg() << endl
         << "#   dip.mom <mu> mu2>      = " << mu << " " << mu2.avg() << endl
         << "#   dip.mom fluctuation    = " << mu2.avg() - pow(mu.avg(),2) << endl
         << "#   charge capacitance     = " << Q2.avg()-pow(Q.avg(),2) << endl
         << "#   (Alternative Q calc)   = " << tit.sumsites() << endl
         << "# Internal prot. energy    = " << u_self.avg() << endl;
    if (movesites==true) 
      cout<<"# Site fluct. (R R2 fluc)  = " << R.avg()<<" "<<R2.avg()<<" "<<R2.avg()-pow(R.avg(),2) << endl;
    if (watchsite!=-1) {
      double A,U;
      A = -log( widom_A.avg() );
      U = widom_Uprot.avg() - widom_Ubulk.avg();
      cout<<"# Prot. process (A,U,TS)   = " << A<<" "<<U<<" "<<U-A<< endl;
    };
    cout << endl;

    s.save( ".coord" + jobid ); //save coordinates

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
       << " <mu>=" << mu.avg()
       << " <mu2>-<mu>2=" << mu2.avg() - pow(mu.avg(),2) << endl;

  hist2.show(muhis,0,500);
  
  //Smear charges and save AAM file to disk
  tit.applycharges(s.p);
  pep.saveAAModel("saved.aam", s.p, g[PROT]);
  
  return 0;

};
