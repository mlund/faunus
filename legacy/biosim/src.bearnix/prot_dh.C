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
#include "chain.h"
#include "titrate.h"
#include "average.h"
#include "vpython.h"

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
  string jobid    = par.getStr("jobid");            //arbitrary job name
  double cell_r   = par.getDbl("cell_r");
  double maxsep   = par.getDbl("maxsep",cell_r/1.5);//restrict protein separation (max)
  double minsep   = par.getDbl("minsep",0.);        //restrict protein separation (min)
  double temp     = par.getDbl("temp");             //temperature (K)
  double dielec   = par.getDbl("dielec",78); 
  double pH       = par.getDbl("pH", 7);

  double vdw      = par.getDbl("hamaker", 0);
  double debye    = par.getDbl("debye");

  string protein1 = par.getStr("protein1"); //protein 1 coords
  string protein2 = par.getStr("protein2"); //protein 1 coords
  double prot_dp  = par.getDbl("prot_dp");
  double prot_rot = par.getDbl("prot_rot");

  bool minsnapshot = par.getBool("minsnapshot",false);
  bool rotation    = par.getBool("rotate", true);

  /**************************************
    Class constructors and declarations
  **************************************/
  enum groupid {P1=0,P2,PROTEINS}; //convenient
  montecarlo mc( macro, micro );
  physconst phys(temp, dielec);
  geometry geo;
  cell cell(cell_r);
  slump slump;
  space s;
  peptide pep;
  interact pot(phys.beta_ecf);
  collision col;
  vector<group> g(PROTEINS+1);

  pot.kappa = 1./debye;
  pot.vdw = vdw;
  
  /**************************************
    Prepare distribution functions
  **************************************/
  vector<histogram::data> h(13);
  histogram hist( 0.5, 2*cell_r );
  hist.init(h);
  h[0].name = "g(r)";     h[0].type = histogram::PROBABILITY;
  h[1].name = "Qp1(r)";   h[1].type = histogram::AVERAGE;
  h[2].name = "Qp2(r)";   h[2].type = histogram::AVERAGE;
  h[3].name = "Utot(r)";   h[3].type = histogram::AVERAGE;
  h[4].name = "UvdW(r)";  h[4].type = histogram::AVERAGE;
  h[5].name = "Dip_y*y";  h[5].type = histogram::AVERAGE;
  h[6].name = "Dip_z*z";  h[6].type = histogram::AVERAGE;
  h[7].name = "U(ionion)"; h[7].type = histogram::AVERAGE;
  h[8].name = "U(iondip)"; h[8].type = histogram::AVERAGE;
  h[9].name = "U(dipdip)"; h[9].type = histogram::AVERAGE;
  h[10].name= "U(p1p2)"  ; h[10].type= histogram::AVERAGE;
  h[11].name= "Dip_P1"   ; h[11].type= histogram::AVERAGE;
  h[12].name= "Dip_P2"   ; h[12].type= histogram::AVERAGE;

  /**********************************
    Load proteins, salt and polymers
  ***********************************/
  g[P1] = s.append( pep.loadAAModel(protein1) );
  g[P1].name = "protein 1";
  g[P1].vdw = true;
  s.addmasscenter(g[P1]);
  s.adddipole(g[P1]);
  s.move( g[P1], -s.center_of_mass(g[P1]), space::AUTOACCEPT);
  s.zmove( g[P1],  maxsep/2., space::AUTOACCEPT);
  g[P1].radius=s.radius( g[P1], s.p[g[P1].cm] );

  g[P2] = s.append( pep.loadAAModel(protein2) );
  g[P2].name = "protein 2";
  g[P2].vdw = true;
  s.addmasscenter( g[P2] );
  s.adddipole(g[P2]);
  s.move( g[P2], -s.center_of_mass(g[P2]), space::AUTOACCEPT);
  s.zmove( g[P2], -maxsep/2., space::AUTOACCEPT);
  g[P2].radius=s.radius( g[P2], s.p[g[P2].cm] );

  g[PROTEINS] = g[P2] + g[P1]; //convenient
  s.trial = s.p;
  s.load(".coord" + jobid); // load saved coordinates and charges.

  //INITIAL ENERGY
  double utot = pot.energy_dh(s.p, g[P1], g[P2])+pot.energy_vdw(s.p, g[P1], g[P2]); //initial energy
  double umin = 1e6 ;                   //minimum system energy

  if (s.safetytest_vector()==false) return 1;

  //PRINT INFORMATION
  cout << "\n# --- SYSTEM ----------------------------------------\n"
       << "# Macro loops         = " << mc.macroSteps << endl
       << "# Micro loops         = " << mc.microSteps << endl
       << "# Total nr. of loops  = " << mc.macroSteps*mc.microSteps << endl
       << "# Cell radius         = " << cell_r << endl
       << "# Temperature (K)     = " << phys.T << endl
       << "# Dielectric Const.   = " << phys.e_r << endl
       << "# Bjerrum length      = " << pot.lB << endl
       << "# Debye length        = " << 1./pot.kappa << endl
       << "# vdW parameter       = " << pot.vdw << endl
       << "# Protein 1           = " << protein1 << endl
       << "# Protein 2           = " << protein2 << endl
       << "# Protein displ.      = " << prot_dp << endl
       << "# Protein rot.        = " << prot_rot <<endl
       << "# Max protein sep.    = " << maxsep << endl
       << "# Min protein sep.    = " << minsep << endl
       << "# System charge       = " << s.charge() << endl
       << "# Initial energy (kT) = " << utot << endl;

  s.charge(g[P1]);
  s.charge(g[P2]);
  cout << "# --- GROUPS ----------------------------------------\n";
  cout << g[P1] 
       << "#   Mass center   = " << s.p[g[P1].cm] << endl
       << "#   Hydrophobics  = " << s.cntHydrophobic(g[P1]) << endl
       << g[P2]
       << "#   Mass center   = " << s.p[g[P2].cm] << endl
       << "#   Hydrophobics  = " << s.cntHydrophobic(g[P2]) << endl
       << endl;

  /************************
    Monte Carlo Simulation
   ************************/
  for (int macroCnt=1; macroCnt<=mc.macroSteps; macroCnt++) {
    for (int i=1; i<=mc.microSteps; i++) {
      bool rc,coll, moved;
      montecarlo::rejectcause rejectcause;
      double u=0,uold=0,unew=0,du;

      //Move proteins
      rc=false;
      coll=false;
      moved=false;
      rejectcause=montecarlo::HC;
      double newz,z, dp = prot_dp*slump.random_half();
      newz = abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z) + 2*dp;
      if (newz <= maxsep && newz >= minsep) {
        s.zmove( g[P1], dp );
        s.zmove( g[P2],-dp );
	moved=true;
	z=newz;
      } else coll=true;

      if (coll==false)
        if (z < 2*cell_r)
          if (col.overlap(s.trial, g[P1], g[P2])==false) {
            uold = pot.energy_dh(s.p, g[P1],g[P2])+pot.energy_vdw(s.p, g[P1], g[P2]);
            unew = pot.energy_dh(s.trial,g[P1],g[P2])+pot.energy_vdw(s.trial, g[P1], g[P2]);
            du=unew-uold;
            if (pot.metropolis(du)==true)
              rc=true;
            else
              rejectcause=montecarlo::ENERGY;
          };
      if (rc==true) {
	s.accept( g[PROTEINS] );
	mc.accept(montecarlo::TRANSLATE);
	utot += du;
      } else {
	if (moved==true)
	  s.undo( g[PROTEINS] );
	mc.reject(montecarlo::TRANSLATE, rejectcause);
      };
      
      z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
      hist.add( h[0], z);
      
      //Rotate protein 1
      if (rotation==true && slump.random_one()>0.5) {
        rc=false;
        coll=false;
        rejectcause=montecarlo::HC;
        s.rotate( g[P1], prot_rot );
	if (col.overlap(s.trial, g[P1], g[P2])==false) {
	  uold=pot.energy_dh(s.p, g[P1], g[P2])+pot.energy_vdw(s.p, g[P1], g[P2]);
	  unew=pot.energy_dh(s.trial, g[P1], g[P2])+pot.energy_vdw(s.trial, g[P1], g[P2]);
	  du = unew-uold;
	  if (pot.metropolis(du)==true)
	    rc=true;
	  else
	    rejectcause=montecarlo::ENERGY;
	};
        if (rc==true) {
          s.accept(g[P1]);
          mc.accept(montecarlo::ROTATE);
          utot += du;
        } else {
          s.undo(g[P1]);
          mc.reject(montecarlo::ROTATE, rejectcause);
        };
        
        z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
        hist.add( h[0], z);
      };

      //Rotate protein 2
      if (rotation==true && slump.random_one()>0.5) {
        rc=false;
        coll=false;
        rejectcause=montecarlo::HC;
        s.rotate( g[P2], prot_rot );
	if (col.overlap(s.trial, g[P1], g[P2])==false) {
	  uold=pot.energy_dh(s.p, g[P1], g[P2])+pot.energy_vdw(s.p, g[P1], g[P2]);
	  unew=pot.energy_dh(s.trial, g[P1], g[P2])+pot.energy_vdw(s.trial, g[P1], g[P2]);
	  du = unew-uold;
	  if (pot.metropolis(du)==true)
	    rc=true;
	  else
	    rejectcause=montecarlo::ENERGY;
	};
        if (rc==true) {
          s.accept(g[P2]);
	  mc.accept(montecarlo::ROTATE);
          utot += du;
        } else {
          s.undo(g[P2]);
          mc.reject(montecarlo::ROTATE, rejectcause);
        };

        z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
        hist.add( h[0], z);
      };

      if (slump.random_one()>0.9) {
        double q1=s.charge(g[P1]);
        double q2=s.charge(g[P2]);
        point mu1=s.calcdipole(g[P1]);
        point mu2=s.calcdipole(g[P2]);
        hist.add(h[7], z, pot.lB*q1*q2/z);
        hist.add(h[8], z, pot.iondip(mu1,q2,z)+pot.iondip(mu2,q1,z));
        hist.add(h[9], z, pot.dipdip(mu1,mu2,z));
        hist.add( h[11],z,s.p[g[P1].dipv].x);
        hist.add( h[12],z,s.p[g[P2].dipv].x);
        hist.add( h[3], z, utot);
        hist.add( h[4], z, pot.energy_vdw(s.p, g[P1], g[P2]));
        bool savevdw=g[1].vdw; g[1].vdw=false; //turn off vdw interactions
        double up1p2=pot.energy(s.p,g[P1],g[P2]);
        hist.add(h[10],z, up1p2);
        g[1].vdw=savevdw;
      };

      // Adjust displacement parameters
      if (slump.random_one()>0.9) {
	mc.adjust_dp(montecarlo::TRANSLATE, prot_dp);
	mc.adjust_dp(montecarlo::ROTATE, prot_rot);
      };

      //Minimum (electrostatic) energy snapshot
      if (utot<umin && minsnapshot==true) {
        umin=utot;
        povray pov; 
        pov.sphere(cell_r); 
        pov.zaxis(cell_r); 
        pov.add( s.p, g[P1] ); 
        pov.add( s.p, g[P2] ); 
        pov.save( ".min.pov" + jobid );
        s.save( ".min.coord" + jobid );
      };
     
    }; //end of micro loop

    if (s.safetytest_vector()==false)
      return 1;

    // System energy check:
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy_dh(s.p, g[P1], g[P2])+pot.energy_vdw(s.p, g[P1], g[P2]);
    cout << "# System energy (kT)       = " << usys << endl
	 << "# System energy drift (kT) = " << utot - usys << endl
	 << "# DP (prot,rotation)       = " << prot_dp<<" " << prot_rot <<endl
	 << "# Protein charges          = " << s.charge(g[P1]) << " " << s.charge(g[P2]) << endl
	 << "# Total system charge      = " << s.charge() << endl
	 << endl;

    s.save( ".coord" + jobid); //save coordinates

    //Graphical output (Povray and vpython)
    povray pov;
    pov.sphere(cell_r);
    pov.zaxis(cell_r);
    pov.add( s.p, g[P1] );
    pov.add( s.p, g[P2] );
    pov.save("snapshot.pov" + jobid);

  }; //end of macro step 

  hist.show(h, 0, cell_r);
  mc.showStatistics();

  return 0;

};
