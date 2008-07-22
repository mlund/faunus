/* Modified for Grand Canonical simulation
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
#include <vector>
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
#include "widommod2.h"
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
  string output   = par.getStr("output");
  double cell_r   = par.getDbl("cell_r");
  double temp     = par.getDbl("temp"); 
  double dielec   = par.getDbl("dielec",78); 
  double pH       = par.getDbl("pH", 7);
  double springk  = par.getDbl("springconst",0.5);
  double catPot   = par.getDbl("catPot", 0.);                                                
  double SaltPot    = par.getDbl("SaltPot", 0.);
  
  string protein1 = par.getStr("protein1");

  double ion_dp   = par.getDbl("ion_dp"); //ion displacement
  int nion1       = par.getInt("nion1");  //number of ion 1, MUST BE CATION
  int nion2       = par.getInt("nion2");  //number of ion 2
  int nion3       = par.getInt("nion3");  //number 
  int chion1      = par.getInt("chion1"); //ionic charges
  int chion2      = par.getInt("chion2");
  int chion3      = par.getInt("chion3");
  int watchsite   = par.getInt("watchsite",-1); //site on which to sample excess chem. pot.

  bool adjust_dp  = par.getBool("adjust_dp", false);
  bool saveaam    = par.getBool("saveAAM", false);
  bool movesites  = par.getBool("movesites", false);
  bool GCtit      = par.getBool("GCtit", false);
  bool scwidom    = par.getBool("scwidom", false);
  bool GCtitEQ    = par.getBool("GCtitEQ", false);
  bool GC         = par.getBool("GC", false);

  if (GC==true) {
    string CatPot    = "CatPot";
    par.setcatPot(CatPot);
    nion1  = par.getInt("nion1");
    nion2  = par.getInt("nion2");

  };
  SaltPot=catPot;
  catPot=catPot/2.;
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
  average Q, Q2, mu, R, R2, N1, N2;
  average widom_A, widom_Uprot, widom_Ubulk;
  average u_self, muex;
  Widommod scwidomCAT(10, 2.0, chion1, "Cation");
  Widommod scwidomAN(10, 2.0, chion2, "Anion");
  pot.k   = springk;
  pot.r0  = 0.0;
  double volume=pow(cell_r,3)*4*phys.pi/3;
  average checkout=0;
  average checkin=0;
  //Distribution functions
  vector<histogram::data> h(2);
  histogram hist(0.5, cell_r);
  hist.init(h);
  h[0].name = "Q" ; h[0].type=histogram::AVERAGE;
  h[1].name = "Q2" ; h[1].type=histogram::AVERAGE;

  vector<histogram::data> ndist(2);
  histogram nhist(1., 1000);
  nhist.init(ndist);
  ndist[0].name ="Cation" ; ndist[0].type=histogram::PROBABILITY;
  ndist[1].name ="Anion"  ; ndist[1].type=histogram::PROBABILITY;

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
  if (GCtit==true) {
    cout<<"# Initial cation" << endl
        <<"# chemical potential  = " << catPot<< endl;
  };
  tit.info();
  
  cout << "# --- GROUPS ----------------------------------------\n";
  cout << g[PROT]
       << "#   Mass center   = " << s.p[g[PROT].cm] << endl
       << "#   Avg. radius   = " << s.radius(g[PROT], s.p[g[PROT].cm], space::AVERAGE) << endl
       << g[SALT] << endl;

  //VMD support
  pdb.load_particles(s.p);
  pdb.save("test.pqr");
  
//  imdwrap imd(s.p.size());
//  imd.wait_for_connection();

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

      //Salt (1:1) excess chemical potential
      if (slump.random_one()>0.9) {
        double m;
        particle anion, cation;
        anion.charge=-1;
        cation.charge=+1;
        anion.radius=cation.radius=2.0;
        for (int k=0; k<100; k++) {
          cell.randomPos(cation);
          cell.randomPos(anion);
          u=1.0e20;
          if (anion.overlap(cation)==false)
            if (col.overlap(s.p, anion)==false)
              if (col.overlap(s.p, cation)==false)
                u = pot.energy(anion, cation) * pot.lB \
                  + pot.potential(s.p, anion) * anion.charge \
                  + pot.potential(s.p, cation)* cation.charge;
          m=exp(-u);
          muex += m;
        };
      };

      
      //Charge regulation
      if (slump.random_one()>0.7 && GCtitEQ==false) {
	titrate::action t;
	for (int cnt=0; cnt<tit.sites.size(); cnt++) {
	  t=tit.exchange(s.trial);
	  uold=pot.energy( s.p, t.site ) + pot.energy( s.p, t.proton )
	    - pot.energy(s.p[t.site], s.p[t.proton])*pot.lB;
	  unew=pot.energy( s.trial, t.site ) + pot.energy( s.trial, t.proton )
	    - pot.energy(s.trial[t.site], s.trial[t.proton])*pot.lB;
	  du=unew-uold;
          if ( GCtit==true) {//If true, Grand Canonical titration
            rc=pot.metropolis( tit.energy(s.trial, du, catPot, t), tit.idPref(t, volume));
            if (t.action==0)
              checkout+= log(tit.check(t,volume,catPot));
            else
              checkin+= log(tit.check(t,volume,catPot));
          }
          else
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
      }

      if (GC==true ) {//&& slump.random_one()>1.2);
        rejectcause=montecarlo::HC;
        rc=false;
        double density;
        if (slump.random_one()<0.5) {
          vector<particle> Saltpair(2); 
          Saltpair[0].charge=chion1; Saltpair[0].radius=2.0; cell.randomPos(Saltpair[0]);
          Saltpair[1].charge=chion2; Saltpair[1].radius=2.0; cell.randomPos(Saltpair[1]); 
          s.compoChangeIn(Saltpair, g[SALT]);
          if (cell.cellCollision(s.trial[g[SALT].end])==false)
          if (cell.cellCollision(s.trial[g[SALT].end-1])==false)
            if (col.overlap(s.trial, g[PROT], g[SALT].end)==false)
            if (col.overlap(s.trial, g[PROT], g[SALT].end-1)==false)
              if (col.overlap(s.trial, g[SALT], g[SALT].end)==false) 
              if (col.overlap(s.trial, g[SALT], g[SALT].end-1)==false) {
              du = pot.energy(s.trial, g[SALT].end)+ pot.energy(s.trial, g[SALT].end-1)-
                   pot.energy(s.trial[g[SALT].end], s.trial[g[SALT].end-1])*pot.lB;
                density=(s.trialcount(g[SALT], 1))/volume*(s.trialcount(g[SALT], -1))/volume;
                if (pot.metropolis(du, SaltPot, 1, density)==true)
                  rc=true;
                else
                  rejectcause=montecarlo::ENERGY;
              };
          if (rc==true) {
            s.p=s.trial;
            tit.protons.push_back(g[SALT].end-1); //titration
            mc.accept(montecarlo::INFLOW);
            utot += du;
            } else {
            s.trial=s.p;
            g[SALT].end=s.p.size()-1;
            mc.reject(montecarlo::INFLOW, rejectcause);
          };
        }else{
          rejectcause=montecarlo::ENERGY;
          rc=false;
          vector<int> Saltpair(2); //First element cation, second anion
          Saltpair=s.compoChangeOut(g[SALT]);
          int phUpdate = Saltpair[0];
          int numP= tit.protons.size();
          for (int findP=0; findP<numP; findP++) { //Locate the position in protons
            if (tit.protons[findP]==phUpdate)
              phUpdate=findP;
          }
          s.trial.erase(s.trial.begin()+Saltpair[0]); 
          Saltpair=s.compoChangeOut(g[SALT]); //Update after first removal
          s.trial.erase(s.trial.begin()+Saltpair[1]);
          g[SALT].end=s.p.size()-1;
          uold=pot.energy(s.p);  //Fix this!!!
          unew=pot.energy(s.trial);
          du = unew-uold; //-(pot.energy(s.p, Saltpair[0])+pot.energy(s.p, Saltpair[1])
                // -pot.energy(s.p[Saltpair[0]], s.p[Saltpair[1]])*pot.lB );
          density=s.count(g[SALT], 1)/volume*s.count(g[SALT], -1)/volume;
          g[SALT].end=s.trial.size()-1;
          if (pot.metropolis(du, SaltPot, -1, density)==true)
            rc=true;
          if (rc==true) {
            s.p=s.trial;
            tit.update(s.p, g[SALT]);   //titration
            mc.accept(montecarlo::OUTFLOW);
            utot += du;
            } else {
            s.trial=s.p;
            mc.reject(montecarlo::OUTFLOW, rejectcause);
            g[SALT].end=s.p.size()-1;

          };
        };
      };

      N1+=s.count(g[SALT], 1); N2+=s.count(g[SALT], -1);
      //Number density within cell
      if (slump.random_one()<0.8) {
        nhist.add(ndist[0], s.count(g[SALT], 1), 1);
        nhist.add(ndist[1], s.count(g[SALT],-1), 1);
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
     
      if (slump.random_one()>0.8 && scwidom==true){
        scwidomCAT.ghostInsert(s.p, pot, col, cell);
        scwidomAN.ghostInsert(s.p, pot, col, cell);
      }                                    
 
      //Adjust displacement parameter
      if (adjust_dp==true && slump.random_one()>0.9) {
        mc.adjust_dp(montecarlo::ION, ion_dp);
      };

     // imd.send_particles(s.p);

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
         << "# Average protein charge   = " << Q << endl
         << "# Protein dipole moment    = " << mu << endl
         << "# Protein capacitance      = " << Q2.avg()-pow(Q.avg(),2) << endl
         << "# (Alternative Q calc)     = " << tit.sumsites() << endl
         << "# <Number> of cations      = " << N1.avg() <<endl
         << "# <Number> of anions       = " << N2.avg() <<endl
         << "# Excess.chem.pot Salt     = " << -log(muex.avg())<<endl;
       if(scwidom==true) {
         
    cout << "# Excess.chem.pot Cation   = " << scwidomCAT.exessChempot() << endl
         << "# Excess.chem.pot Anion    = " << scwidomAN.exessChempot() << endl;
         };
    cout << "# Internal prot. energy    = " << u_self.avg() << endl;
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
       << " <mu>=" << mu.avg() <<" <N1>="<<N1.avg()
       << " <N2>=" <<N2.avg()<< " checkout="<<checkout.avg()<<" checkin= "<<checkin.avg()<<endl;
  if (scwidom==true) {
    scwidomCAT.getWidomResults();
    scwidomAN.getWidomResults();
  };
  cout << "# Excess.chem.pot Salt     = " << -log(muex.avg())<<endl;
  //Smear charges and save AAM file to disk
  tit.applycharges(s.p);
  pep.saveAAModel("saved.aam", s.p, g[PROT]);
 
  if (scwidom==true ) {
    ofstream f( "CatPot" );
    if (f) {
      f << scwidomCAT.exessChempot() << endl
        << g[SALT].end-g[SALT].beg+1-s.count(g[SALT], -1) <<endl
        << s.count(g[SALT], -1) << endl;
      f.close();
    }
    else { cout << "*** Failed to promt chemical potential **" <<endl;
    };
  };

  nhist.show(ndist, 0, 1000);

  return 0;
};

