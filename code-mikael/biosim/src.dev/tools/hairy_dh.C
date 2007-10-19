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

/*
 * MONTE CARLO ALGORITME:
 * - to vektorer (gml, ny)
 * - forflytning goeres i "ny".
 * - overlap: kopier gml til ny. afslut.
 * - beregn gammel og ny energi
 * - metropolis:
 *   - accepter: kopier aendringer fra ny til gml. afslut.
 *   - forkast: kopier fra gml til ny. afslut
 *
 * 
 * VI GEMMER IKKE COLLISIONS VEKTOREN (sqdist)
 * collision:
 *   collision(group, ny, gml)
 *   collision(particle, ny, gml)
 * energy
 *   energy( salt, Achg)
 *   energy( salt, Bchg)
 *   energy( Achg, Bchg)
 *   energy_vdw( A, B )
 */

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
  double maxsep   = par.getDbl("maxsep",cell_r/1.5); //restrict protein separation (max)
  double minsep   = par.getDbl("minsep",0.);        //restrict protein separation (min)
  double temp     = par.getDbl("temp");             //temperature (K)
  double dielec   = par.getDbl("dielec",78); 
  double debye    = par.getDbl("debye",1000);
  
  double vdw      = par.getDbl("hamaker", 0);
  double springk  = par.getDbl("springconst",0.1);
  double springeq = par.getDbl("springeqdist",10.);
  double u_penalty= par.getDbl("penalty",0.);
  
  string protein1 = par.getStr("protein1"); //protein 1 coords
  string protein2 = par.getStr("protein2"); //protein 1 coords
  double prot_dp  = par.getDbl("prot_dp");
  double prot_rot = par.getDbl("prot_rot");
  double dp_monomer=par.getDbl("monomer_dp",2.); //monomer displacement factor
  
  bool hairy       = par.getBool("hairy", true);
  bool rotation    = par.getBool("rotate", true);
  bool minsnapshot = par.getBool("minsnapshot",false);
  bool adjust_dp   = par.getBool("adjust_dp", false);
  
  /**************************************
    Class constructors and declarations
    **************************************/
  enum groupid {P1=0,P2,SALT,
    CHAIN1,CHAIN2,CHAIN3,CHAIN4,CHAIN5,CHAIN6,
    CHAIN7,CHAIN8,CHAIN9,CHAIN10,CHAIN11,CHAIN12,ENDCHAIN,CHAINS,
    PROTEINS,P1CHAINS,P2CHAINS,LAST};
  montecarlo mc( macro, micro );
  physconst phys(temp, dielec);
  geometry geo;
  cell cell(cell_r);
  slump slump;
  space s;
  peptide pep;
  interact pot(phys.beta_ecf);
  collision col;
  vector<group> g(LAST+1);
  
  pot.vdw = vdw;
  pot.k   = springk;
  pot.r0  = springeq;
  pot.kappa=1./debye;
  
  /**************************************
    Prepare distribution functions
    **************************************/
  vector<histogram::data> h(1);
  histogram hist( 0.5, 2*cell_r );
  hist.init(h);
  h[0].name = "g(r)";     h[0].type = histogram::PROBABILITY;
  
  /**********************************
    Load proteins, salt and polymers
    ***********************************/
  g[P1] = s.append( pep.loadAAModel(protein1) );
  g[P1].name = "protein 1";
  g[P1].vdw = true;
  s.addmasscenter(g[P1]);
  s.adddipole(g[P1]);
  s.move( g[P1], -s.center_of_mass(g[P1]), space::AUTOACCEPT);
  s.zmove( g[P1],  (maxsep/2.-10), space::AUTOACCEPT);
  g[P1].radius=s.radius( g[P1], s.p[g[P1].cm] );
  
  g[P2] = s.append( pep.loadAAModel(protein2) );
  g[P2].name = "protein 2";
  g[P2].vdw = true;
  s.addmasscenter( g[P2] );
  s.adddipole(g[P2]);
  s.move( g[P2], -s.center_of_mass(g[P2]), space::AUTOACCEPT);
  s.zmove( g[P2], -(maxsep/2.-10), space::AUTOACCEPT);
  g[P2].radius=s.radius( g[P2], s.p[g[P2].cm] );
  
  g[PROTEINS] = g[P2] + g[P1]; //convenient
  
  if (hairy==true) {
    int graftpt=0;
    for (int i=CHAIN1; i<ENDCHAIN; i++) {
      graftpt++;
      if (i==CHAIN7)
        graftpt+=3;
      //cout << graftpt << endl;
      //graftpt=-1;
      chain c("6 CATION CATION CATION CATION CATION CATION",graftpt);
      g[i]=s.insert_chain(c);
      g[CHAINS]+=g[i]; //convenient...
    };
    g[P1CHAINS] = g[CHAIN1] + g[CHAIN2] + g[CHAIN3] + g[CHAIN4] + g[CHAIN5] + g[CHAIN6];
    g[P2CHAINS] = g[CHAIN7] + g[CHAIN8] + g[CHAIN9] + g[CHAIN10] + g[CHAIN11] + g[CHAIN12];
  };
  
  s.load(".coord" + jobid); // load saved coordinates and charges.
  
  //INITIAL ENERGY
  double utot = pot.energy_dh(s.p);
  for (int k=CHAIN1; k<ENDCHAIN; k++)
    utot += pot.internal(s.p, g[k]);
  double umin = 1e6 ;                   //minimum system energy
  cout << "Internal = " << utot << endl;
  
  if (s.safetytest_vector()==false) return 1;
  cout << "Pass 3.\n";
  
  //PRINT INFORMATION
  cout << "\n# --- SYSTEM ----------------------------------------\n"
    << "# Steps               = " << mc.macroSteps << " x " << mc.microSteps
    << " = " << mc.macroSteps*mc.microSteps << endl
    << "# Cell radius         = " << cell_r << endl
    << "# Temperature (K)     = " << phys.T << endl
    << "# Jobid               = " << jobid << endl
    << "# Dielectric Const.   = " << phys.e_r << endl
    << "# Bjerrum length      = " << pot.lB << endl
    << "# Debye length        = " << 1./pot.kappa << endl
    << "# vdW parameter       = " << pot.vdw << endl
    << "# Protein 1           = " << protein1 << endl
    << "# Protein 2           = " << protein2 << endl
    << "# Protein displ.      = " << prot_dp << endl
    << "# Protein rot.        = " << prot_rot <<endl
    << "# Max/min separation  = " << minsep << " " << maxsep << endl
    << "# System charge       = " << s.charge() << endl
    << "# Initial energy (kT) = " << utot << endl;
 
  s.charge(g[P1]);
  s.charge(g[P2]);
  s.charge(g[PROTEINS]);
  s.charge(g[CHAINS]);
  
  cout << "# --- GROUPS ----------------------------------------\n";
  cout << g[P1] 
    << "#   Mass center   = " << s.p[g[P1].cm] << endl
    << g[P2]
    << "#   Mass center   = " << s.p[g[P2].cm] << endl
    << g[PROTEINS] << g[CHAINS] << endl;
  for (int k=CHAIN1; k<ENDCHAIN; k++)
    cout << g[k] << endl;
  
  /************************
    Monte Carlo Simulation
    ************************/
  vector<particle> tmp;
  tmp = s.p;
  
  for (int macroCnt=1; macroCnt<=mc.macroSteps; macroCnt++) {
    for (int i=1; i<=mc.microSteps; i++) {
      bool rc,coll;
      montecarlo::rejectcause rejectcause;
      double u=0,uold=0,unew=0,du;
      
      //Rattle chains
      if (hairy==true) {
        for (int j=CHAIN1; j<ENDCHAIN; j++) {
          for (int i=0; i<g[j].size(); i++) {
            rc=false;
            rejectcause=montecarlo::HC;
            
            int n=g[j].random();
            s.displace(n, dp_monomer);
            if (cell.cellCollision(s.trial[n])==false)
              if (col.overlap(s.trial, n)==false) {
                uold=pot.chain(s.p, g[j], n);
                unew=pot.chain(s.trial, g[j], n);
                if (s.p[n].charge!=0) {
                  uold+=pot.energy_dh(s.p, n);
                  unew+=pot.energy_dh(s.trial, n);
                };
                du = unew-uold;
                if (pot.metropolis(du)==true)
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
          };
        };
      };
      
      //Move proteins
      if (slump.random_one()>0.2) {
        rc=false;
        coll=false;
        bool moved=false;
        rejectcause=montecarlo::HC;
        double oldz,z;
        double dp = prot_dp*slump.random_half();
        oldz = abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
        z    = oldz + 2.*dp;
        if (dp!=0 && z <= maxsep && z >= minsep && z<2.*cell_r) {
          s.zmove( g[P1], dp );
          s.zmove( g[P2],-dp );
          moved=true;
        } else coll=true;
        
        if (coll==false)
          if (col.overlap(s.trial, g[PROTEINS])==false)
            if (col.overlap(s.trial, g[P1], g[P2])==false) {
              uold = pot.energy_dh(s.p, g[PROTEINS]) + pot.energy_dh(s.p, g[P1],g[P2]);
              unew = pot.energy_dh(s.trial,g[PROTEINS]) + pot.energy_dh(s.trial,g[P1],g[P2]);
              if (hairy==true) {
                for (int k=CHAIN1; k<ENDCHAIN; k++) {
                  uold += pot.graft(s.p, g[k]);
                  unew += pot.graft(s.trial, g[k]);
                };
              };
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
        hist.add( h[0], z);             //update radial distribution function
      };
      
      //Rotate proteins
      for (int k=P1; k<=P2; k++) {
        double z;
        if (rotation==true && slump.random_one()>0.5) {
          rc=false;
          coll=false;
          rejectcause=montecarlo::HC;
          s.rotate( g[k], prot_rot );
          
          if (hairy==true)
            coll=col.overlap(s.trial, g[k], g[CHAINS]);
          if (coll==false)
            if (col.overlap(s.trial, g[P1], g[P2])==false) {
              uold=pot.energy_dh(s.p, g[k]);
              unew=pot.energy_dh(s.trial, g[k]);
              if (hairy==true)
                for (int l=CHAIN1; l<ENDCHAIN; l++) {
                  uold += pot.graft(s.p, g[l]);
                  unew += pot.graft(s.trial, g[l]);
                };
              du = unew-uold;
              if (pot.metropolis(du)==true)
                rc=true;
              else
                rejectcause=montecarlo::ENERGY;
            };
          if (rc==true) {
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
      
      //Adjust displacement parameter
      if (adjust_dp==true && slump.random_one()>0.9) {
        mc.adjust_dp(montecarlo::TRANSLATE, prot_dp);
        mc.adjust_dp(montecarlo::ROTATE, prot_rot);
        mc.adjust_dp(montecarlo::MONOMER, dp_monomer);
      };    
    }; //end of micro loop
    
    // Check if particle vectors are in sync
    if (s.safetytest_vector()==false)
      return 1; //teminate program
    
    // System energy check:
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy_dh(s.p);
    for (int k=CHAIN1; k<ENDCHAIN; k++)
      usys += pot.internal(s.p, g[k]);
    
    cout << "# System energy (kT)      = " << usys << endl
      << "# System energy drift (kT) = " << utot - usys << endl
      << "# DP (prot,rot,chain)      = " << prot_dp<<" "<<prot_rot<<" "<<dp_monomer<<endl
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
    if (hairy==true)
      for (int i=CHAIN1; i<ENDCHAIN; i++)
        pov.add( s.p, g[i] );
    pov.save("snapshot.pov" + jobid);
    
    vpython vpy;
    vpy.add( s.p, g[P1] );
    vpy.add( s.p, g[P2] );                                           
    if (hairy==true)
      for (int i=CHAIN1; i<ENDCHAIN; i++)
        vpy.add( s.p, g[i] );                                                                  
    vpy.save("vpython" + jobid + ".py");
    
    //test for any static particles...
    for (int i=0; i<tmp.size(); i++)
      if (abs(tmp[i].x-s.p[i].x)<0.01)
        if (abs(tmp[i].y-s.p[i].y)<0.01)
          if (abs(tmp[i].z-s.p[i].z)<0.01)
            cout << "Possible static particle: " << i << endl;
    
    
  }; //end of macro step 
  
  /**************************************
    Print results and terminate program
    **************************************/
  hist.show(h, minsep, maxsep);
  mc.showStatistics();
  
  return 0;
};
