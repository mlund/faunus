#include <iostream>
#include "container.h"
#include "potentials.h"
#include "group.h"
#include "point.h"
#include "space.h"
#include "slump.h"
#include "species.h"
#include "inputfile.h"
#include "montecarlo.h"
#include "povray.h"
#include "average.h"
#include "io.h"
#include "histogram.h"
#include "physconst.h"
#include "markovmove.h"

/*
 * Monte Carlo simulation of Protein-Protein Interactions
 * Mikael Lund, 2003-2007.
 */

using namespace std;

int main(int argc, char* argv[] ) { 

  // handle arguments
  string cfgfile="input.cfg";
  if (argc!=2) {
    cout << "Usage: " << argv[0] << " input.cfg\n";
    //return -1;
  } else cfgfile = argv[1];
 
  // class constructors and declarations
  enum groupid {P1=0,P2,SALT,LAST};
  config c(cfgfile);
  montecarlo mc( c.macro, c.micro );
  physconst phys( c.temp, c.dielec );
  cell cell( c.cell_r );
  slump slump;
  io io;
  space s;
  species spc;

  pot_setup potcfg;
  potcfg.lB=7.1;
  potcfg.eps=0.2;
  interaction<T_pairpot> pot(potcfg);

  hardsphere hd;
  vector<group> g(LAST+1);
  average<double> systemenergy;
  rdf rdf_nacl(particle::NA, particle::CL, 0, c.maxsep);

  // prepare distribution functions
  enum distenum {GOFR=0, UTOT, DLAST};
  //vector<histogram [0.5,0,100]> d(2);
  //d[GOFR].set(0.5, "g(r)", true);
  //d[UTOT].set(0.5, "Total energy");

  // load proteins, salt
  //g[P1] = s.append( io.loadaam(spc, c.protein1) );
  g[P1].name = "Protein 1";
  s.mass_center(g[P1]);
  s.move( g[P1], -g[P1].cm, space::AUTOACCEPT);
  s.zmove( g[P1],  (c.maxsep/2.-0.5), space::AUTOACCEPT);
  g[P1].radius=s.radius( g[P1], g[P1].cm );

  //g[P2] = s.append( io.loadaam(spc, c.protein2) );
  g[P2].name = "Protein 2";
  s.mass_center(g[P2]);
  s.move( g[P2], -g[P2].cm, space::AUTOACCEPT);
  s.zmove( g[P2], -(c.maxsep/2.-0.5), space::AUTOACCEPT);
  g[P2].radius=s.radius( g[P2], g[P2].cm );
  
  g[SALT]  = s.insert_salt( c.nion3, spc.get(c.tion3), cell); 
  g[SALT] += s.insert_salt( c.nion1, spc.get(c.tion1), cell);
  g[SALT] += s.insert_salt( c.nion2, spc.get(c.tion2), cell);
  g[SALT].name="Ions";

  #ifdef POT_DATAPMF
  cout << c.pmfdir << endl;
  pot.loadpmf(spc, c.pmfdir, c.tion1);
  pot.loadpmf(spc, c.pmfdir, c.tion2);
  pot.loadpmf(spc, c.pmfdir, c.tion3);
  pot.showpmf(spc);
  #endif

  s.load(".coord" + c.jobid); // load saved coordinates and charges.
  s.mass_center(g[P1]);
  s.mass_center(g[P2]);
  s.mass_center(g[SALT]);
  
  // initial energy
  double utot = pot.energy(s.p);
  if (s.safetytest_vector()==false)
    return 1;
  
  // userinfo
  cout << "\n# --- SYSTEM ----------------------------------------\n"
    << "# Steps               = " << mc.macroSteps << " x " << mc.microSteps
    << " = " << mc.macroSteps*mc.microSteps << endl
    << "# Cell radius         = " << c.cell_r << endl
    << "# Temperature (K)     = " << phys.T << endl
    << "# Jobid               = " << c.jobid << endl
    << "# Dielectric Const.   = " << phys.e_r << endl
    << "# Bjerrum length      = " << pot.pair.f << endl
    << "# Protein 1           = " << c.protein1 << endl
    << "# Protein 2           = " << c.protein2 << endl
    << "# Max/min separation  = " << c.minsep << " " << c.maxsep << endl
    << "# Ion 1 (# type z r)  = " << c.nion1 << " " << c.tion1 << endl
    << "# Ion 2 (# type z r)  = " << c.nion2 << " " << c.tion2 << endl
    << "# Ion 3 (# type z r)  = " << c.nion3 << " " << c.tion3 << endl
    << "# System charge       = " << s.charge() << endl
    << "# Initial energy (kT) = " << utot << endl;
  s.charge(g[P1]);
  s.charge(g[P2]);
  s.charge(g[SALT]);
  cout << "\n# --- GROUPS ----------------------------------------\n";
  cout << g[P1] << g[P2] << g[SALT] << endl;
 
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
        s.displace(n, c.ion_dp);   //displace it...
        if (cell.collision(s.trial[n])==false) {
          #pragma omp parallel
          {
            #pragma omp sections
            { 
              #pragma omp section
              { uold = pot.energy(s.p, n);   }
              #pragma omp section
              { unew = pot.energy(s.trial,n);   }
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
          systemenergy += utot;
        } else {
          s.trial[n]=s.p[n];
          mc.reject(montecarlo::ION, rejectcause);
        };
        //rdf_nacl.update(s.p);
      };
      
      // move proteins
      if (slump.random_one()>0.2) {
        rc=false;
        coll=false;
        bool moved=false;
        rejectcause=montecarlo::HC;
        double oldz,z,dp;
        dp   = c.prot_dp*slump.random_half();
        oldz = abs(g[P1].cm.z-g[P2].cm.z);
        z    = oldz + 2.*dp;
        if (z <= c.maxsep && z >= c.minsep && z<2.*c.cell_r) {
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
              { uold = pot.energy(s.p,g[SALT]) \
                + pot.energy(s.p,g[P1],g[P2]);}
              #pragma omp section
              { unew = pot.energy(s.trial,g[SALT]) \
                + pot.energy(s.trial,g[P1],g[P2]);}
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
          systemenergy+=utot;
        } else {
          mc.reject(montecarlo::TRANSLATE, rejectcause);
          if (moved==true) {
            s.undo( g[P1] );
            s.undo( g[P2] );
          };
        };
        z=abs(g[P1].cm.z-g[P2].cm.z);
        //d[ GOFR ].add(z); 
        //d[ UTOT ].add(z, utot);
      };
      
      // rotate proteins
      for (unsigned int k=P1; k<=P2; k++) {
        double z;
        if (c.rotate==true && slump.random_one()>0.5) {
          rc=false;
          rejectcause=montecarlo::ENERGY;
          s.rotate( g[k], c.prot_rot );
          #pragma omp parallel
          {
            #pragma omp sections
            {
              #pragma omp section
              { uold=pot.energy(s.p, g[k]); }
              #pragma omp section
              { unew=pot.energy(s.trial, g[k]); }
            }
          }
          du = unew-uold;
          if (pot.metropolis(du)==true) {
            s.accept(g[k]);
            s.recalc_dipole(g[k]);
            mc.accept(montecarlo::ROTATE);
            utot += du;
            systemenergy+=utot;
          } else {
            s.undo(g[k]);
            mc.reject(montecarlo::ROTATE, rejectcause);
          };
          z=abs(g[P1].cm.z-g[P2].cm.z);
          //d[ GOFR ].add(z);
          //d[ UTOT ].add(z, utot);
        };
      };
      
      // adjust displacement parameter
      if (c.adjust_dp==true && slump.random_one()>0.9) {
        double min=50, max=60;
        mc.adjust_dp(montecarlo::ION, c.ion_dp, min, max);
        mc.adjust_dp(montecarlo::TRANSLATE, c.prot_dp, min, max);
        mc.adjust_dp(montecarlo::ROTATE, c.prot_rot, min, max);
      };    

    }; // end of micro loop

    if (s.safetytest_vector()==false)
      return 1;
    if (g[P1].cm.z != g[P1].cm_trial.z)
      return 1;

    
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy(s.p) ;
    cout << "# System energy (kT)    = " << usys << endl
      << "# Avg. sys energy (kT)     = " << systemenergy.avg() << endl
      << "# Energy drift (rel abs)   = " << utot - usys << " " << (utot-usys)/usys << endl
      << "# DP (ion,prot,rot)        = " << c.ion_dp<<" "
      << c.prot_dp<<" "<<c.prot_rot<<" "<<endl
      << "# Configs/second           = " << endl//mc.speed() << endl
      << "# Protein charges          = " << s.charge(g[P1]) << " " << s.charge(g[P2]) << endl
      << "#   Mass-centers           = " << g[P1].cm << " " << g[P2].cm << endl
      << "# Total system charge      = " << s.charge() << endl
      << endl;
    s.save( ".coord" + c.jobid); //save coordinates
    
    //Povray output
    povray pov;
    pov.sphere(c.cell_r);
    pov.zaxis(c.cell_r);
    pov.add( s.p, g[P1] );
    pov.add( s.p, g[P2] );
    pov.add( s.p, g[SALT] );
    pov.save("snapshot.pov" + c.jobid);
    cout << endl;
    
  }; //end of macro step 
  
  //cout << d[ GOFR ].show();
  mc.showStatistics();
  s.save(".coord" + c.jobid);
  //rdf_nacl.show();
  return 0;
};
