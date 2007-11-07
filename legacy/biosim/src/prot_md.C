#include <iostream>
#include "group.h"
#include "point.h"
#include "space.h"
#include "slump.h"
#include "interact.h"
#include "peptide.h"
#include "inputfile.h"
#include "montecarlo.h"
#include "histogram.h"
#include "povray.h"
#include "average.h"
#include "intrinsics.h"

/*
 * Monte Carlo simulation of Protein-Protein Interactions
 * Mikael Lund, 2003-2007.
 */

using namespace std;

int main(int argc, char* argv[] ) { 
  // Handle arguments and input file
  string inputName = "amino.inp";
  if (argc==2)
    inputName = argv[1];
  inputfile par( inputName );
  int randomseed  = par.getInt("randomseed",13);
  int macro       = par.getInt("macro");
  int micro       = par.getInt("micro");
  string jobid    = par.getStr("jobid");            //arbitrary job name
  double cell_r   = par.getDbl("cell_r");
  double maxsep   = par.getDbl("maxsep",cell_r/1.5); //restrict protein separation (max)
  double minsep   = par.getDbl("minsep",0.);        //restrict protein separation (min)
  double temp     = par.getDbl("temp");             //temperature (K)
  double dielec   = par.getDbl("dielec",78); 
  string protein1 = par.getStr("protein1"); //protein 1 coords
  string protein2 = par.getStr("protein2"); //protein 1 coords
  int nion1       = par.getInt("nion1");  //number of ion 1
  int nion2       = par.getInt("nion2");  //number of ion 2
  int nion3       = par.getInt("nion3");  //number 
  int chion1      = par.getInt("chion1"); //ionic charges
  int chion2      = par.getInt("chion2");
  int chion3      = par.getInt("chion3");
  double prot_dp  = par.getDbl("prot_dp");
  double prot_rot = par.getDbl("prot_rot");
  double ion_dp   = par.getDbl("ion_dp"); //ion displacement
  string tion1    = par.getStr("tion1");
  string tion2    = par.getStr("tion2");
  string tion3    = par.getStr("tion3");
  string pmfdir   = par.getStr("pmfdir");
  bool rotation    = par.getBool("rotate", true);
  bool adjust_dp   = par.getBool("adjust_dp", false);
 
  // class constructors and declarations
  enum groupid {P1=0,P2,SALT,PROTEINS,LAST};
  montecarlo mc( macro, micro );
  physconst phys(temp, dielec);
  cell cell(cell_r);
  slump slump;
  space s;
  peptide pep;
  interact pot(phys.beta_ecf);
  collision col;
  vector<group> g(LAST+1);
 
  // prepare distribution functions
  vector<histogram::data> h(1);
  histogram hist( 0.5, 2*cell_r );
  hist.init(h);
  h[0].name = "g(r)"; h[0].type = histogram::PROBABILITY;

  // load proteins, salt
  g[P1] = s.append( pep.loadAAModel(protein1) );
  g[P1].name = "protein 1";
  s.addmasscenter(g[P1]);
  s.adddipole(g[P1]);
  s.move( g[P1], -s.center_of_mass(g[P1]), space::AUTOACCEPT);
  s.zmove( g[P1],  (maxsep/2.-10), space::AUTOACCEPT);
  g[P1].radius=s.radius( g[P1], s.p[g[P1].cm] );
  
  g[P2] = s.append( pep.loadAAModel(protein2) );
  g[P2].name = "protein 2";
  s.addmasscenter( g[P2] );
  s.adddipole(g[P2]);
  s.move( g[P2], -s.center_of_mass(g[P2]), space::AUTOACCEPT);
  s.zmove( g[P2], -(maxsep/2.-10), space::AUTOACCEPT);
  g[P2].radius=s.radius( g[P2], s.p[g[P2].cm] );
  
  g[PROTEINS] = g[P2] + g[P1]; //convenient
  
  g[SALT]  = s.insert_salt( nion3, chion3, 2.0, cell, pep.getId(tion3)); 
  g[SALT] += s.insert_salt( nion1, chion1, 2.0, cell, pep.getId(tion1));
  g[SALT] += s.insert_salt( nion2, chion2, 2.0, cell, pep.getId(tion2));
  g[SALT].name="mobile ions";
 
  pep.loadpmf(pmfdir, tion1);
  pep.loadpmf(pmfdir, tion2);
  pep.loadpmf(pmfdir, tion3);
  pep.showpmf();

  s.load(".coord" + jobid); // load saved coordinates and charges.
  
  // initial energy
  double utot = pep.energy(s.p);
  if (s.safetytest_vector()==false)
    return 1;
  
  // userinfo
  cout << "\n# --- SYSTEM ----------------------------------------\n"
    << "# Steps               = " << mc.macroSteps << " x " << mc.microSteps
    << " = " << mc.macroSteps*mc.microSteps << endl
    << "# Cell radius         = " << cell_r << endl
    << "# Temperature (K)     = " << phys.T << endl
    << "# Jobid               = " << jobid << endl
    << "# Dielectric Const.   = " << phys.e_r << endl
    << "# Bjerrum length      = " << pot.lB << endl
    << "# Protein 1           = " << protein1 << endl
    << "# Protein 2           = " << protein2 << endl
    << "# Max/min separation  = " << minsep << " " << maxsep << endl
    << "# Ion 1 (# z type)    = " << nion1<<" "<< chion1 << " " << tion1 << endl
    << "# Ion 2 (# z type)    = " << nion2<<" "<< chion2 << " " << tion2 << endl
    << "# Ion 3 (# z type)    = " << nion3<<" "<< chion3 << " " << tion3 << endl
    << "# System charge       = " << s.charge() << endl
    << "# Initial energy (kT) = " << utot << endl;
  s.charge(g[P1]);
  s.charge(g[P2]);
  s.charge(g[SALT]);
  cout << "# --- GROUPS ----------------------------------------\n";
  cout << g[P1] 
    << "#   Mass center   = " << s.p[g[P1].cm] << endl << g[P2]
    << "#   Mass center   = " << s.p[g[P2].cm]<<endl<< g[SALT] << g[PROTEINS] << endl;
 
  //
  // Markov chain - main MC loop
  //
  for (unsigned int macroCnt=1; macroCnt<=mc.macroSteps; macroCnt++) {
    for (unsigned int i=1; i<=mc.microSteps; i++) {
      bool rc,coll;
      montecarlo::rejectcause rejectcause;
      double u=0,uold=0,unew=0,du=0;
      
      // move salt
      for (unsigned int j=0; j<g[SALT].size(); j++) {
        rc=false;
        coll=false;
        rejectcause=montecarlo::HC;
        int n=g[SALT].random();  //pick random particle
        s.displace(n, ion_dp);   //displace it...
        if (cell.cellCollision(s.trial[n])==false) {
          #pragma omp parallel
          {
            #pragma omp sections
            { 
              #pragma omp section
              { uold = pep.energy(s.p, n);   }
              #pragma omp section
              { unew = pep.energy(s.trial,n);   }
            }
          }
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
      
      // move proteins
      if (slump.random_one()>0.2) {
        rc=false;
        coll=false;
        bool moved=false;
        rejectcause=montecarlo::HC;
        double oldz,z;
        double dp = prot_dp*slump.random_half();
        oldz = abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
        z    = oldz + 2.*dp;
        if (z <= maxsep && z >= minsep && z<2.*cell_r) {
          s.zmove( g[P1], dp );
          s.zmove( g[P2],-dp );
          moved=true;
        } else coll=true;
        
        if (coll==false) {
          #pragma omp parallel
          {
            #pragma omp sections
            {
              #pragma omp section
              { uold = pep.energy(s.p,g[PROTEINS]) \
                + pep.energy(s.p,g[P1],g[P2]);}
              #pragma omp section
              { unew = pep.energy(s.trial,g[PROTEINS]) \
                + pep.energy(s.trial,g[P1],g[P2]);}
            }
          }
          du=unew-uold;
          if (pot.metropolis(du)==true)
            rc=true;
          else
            rejectcause=montecarlo::ENERGY;
        };
        if (rc==true) {
          s.accept( g[P1] );
          s.accept( g[P2] );
          mc.accept(montecarlo::TRANSLATE);
          utot += du;
        } else {
          mc.reject(montecarlo::TRANSLATE, rejectcause);
          if (moved==true) {
            s.undo( g[P1] );
            s.undo( g[P2] );
          };
        };
        z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
        hist.add( h[0], z); 
      };
      
      // rotate proteins
      for (unsigned int k=P1; k<=P2; k++) {
        double z;
        if (rotation==true && slump.random_one()>0.5) {
          rc=false;
          rejectcause=montecarlo::ENERGY;
          s.rotate( g[k], prot_rot );
          #pragma omp parallel
          {
            #pragma omp sections
            {
              #pragma omp section
              { uold=pep.energy(s.p, g[k]); }
              #pragma omp section
              { unew=pep.energy(s.trial, g[k]); }
            }
          }
          du = unew-uold;
          if (pot.metropolis(du)==true) {
            s.accept(g[k]);
            mc.accept(montecarlo::ROTATE);
            utot += du;
          } else {
            s.undo(g[k]);
            mc.reject(montecarlo::ROTATE, rejectcause);
          };
          z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
          hist.add( h[0], z);
        };
      };
      
      // adjust displacement parameter
      if (adjust_dp==true && slump.random_one()>0.9) {
        mc.adjust_dp(montecarlo::ION, ion_dp);
        mc.adjust_dp(montecarlo::TRANSLATE, prot_dp);
        mc.adjust_dp(montecarlo::ROTATE, prot_rot);
      };    

    }; // end of micro loop
    
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pep.energy(s.p) ;
    cout << "# System energy (kT)    = " << usys << endl
      << "# Energy drift (rel abs)   = " << utot - usys << " " << (utot-usys)/usys << endl
      << "# DP (ion,prot,rot)        = " << ion_dp<<" "<<prot_dp<<" "<<prot_rot<<" "<<endl
      << "# Protein charges          = " << s.charge(g[P1]) << " " << s.charge(g[P2]) << endl
      << "# Total system charge      = " << s.charge() << endl
      << endl;
    
    s.save( ".coord" + jobid); //save coordinates
    
    //Povray output
    povray pov;
    pov.sphere(cell_r);
    pov.zaxis(cell_r);
    pov.add( s.p, g[P1] );
    pov.add( s.p, g[P2] );
    pov.add( s.p, g[SALT] );
    pov.save("snapshot.pov" + jobid);
    cout << endl;
    
  }; //end of macro step 
  
  hist.show(h, minsep, maxsep);
  mc.showStatistics();
  s.save(".coord" + jobid);
  return 0;
};
