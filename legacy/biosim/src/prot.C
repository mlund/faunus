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
#include "pdbclass.h"
#include "imdwrap.h"
#include "xyz.h"

/*
 * Monte Carlo simulation of Protein-Protein Interactions
 * Mikael Lund, 2003-2007.
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
  double pH       = par.getDbl("pH", 7);
  
  double vdw      = par.getDbl("hamaker", 0);
  double springk  = par.getDbl("springconst",0.1);
  double springeq = par.getDbl("springeqdist",10.);
  double u_penalty= par.getDbl("penalty",0.);
  
  string protein1 = par.getStr("protein1"); //protein 1 coords
  string protein2 = par.getStr("protein2"); //protein 1 coords
  double prot_dp  = par.getDbl("prot_dp");
  double prot_rot = par.getDbl("prot_rot");
  double dp_monomer=par.getDbl("monomer_dp",2.); //monomer displacement factor
  double ion_dp   = par.getDbl("ion_dp"); //ion displacement
  double clust_dp = 10.;
  
  int nion1       = par.getInt("nion1");  //number of ion 1
  int nion2       = par.getInt("nion2");  //number of ion 2
  int nion3       = par.getInt("nion3");  //number 
  int chion1      = par.getInt("chion1"); //ionic charges
  int chion2      = par.getInt("chion2");
  int chion3      = par.getInt("chion3");
  double rion3    = par.getDbl("rion3", 2.);  //radius of ion 3
  
  bool hairy       = par.getBool("hairy", true);
  bool rotation    = par.getBool("rotate", true);
  bool titrateBool = par.getBool("titrate", true);
  bool smear       = par.getBool("smear", false);
  bool minsnapshot = par.getBool("minsnapshot",false);
  bool adjust_dp   = par.getBool("adjust_dp", false);
  bool imdBool     = par.getBool("imdsupport",false);
  
  /**************************************
    Class constructors and declarations
    **************************************/
  enum groupid {P1=0,P2,SALT,
    CHAIN1,CHAIN2,CHAIN3,CHAIN4,CHAIN5,CHAIN6,
    CHAIN7,CHAIN8,CHAIN9,CHAIN10,CHAIN11,CHAIN12,ENDCHAIN,CHAINS,
    PROTEINS,P1CHAINS,P2CHAINS,LAST};
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
  vector<group> g(LAST+1);
  
  average Q1,Q2, dip1,dip2;
  average Q1sq, Q2sq;
  average muex;                 //salt excess chemical potential (1:1)
  average mu1_z, mu2_z;         //dipole z components
  
  pot.vdw = vdw;
  pot.k   = springk;
  pot.r0  = springeq;
  
  /**************************************
    Prepare distribution functions
    **************************************/
  vector<histogram::data> h(23);
  histogram hist( 0.5, 2*cell_r );
  hist.init(h);
  h[0].name = "g(r)";     h[0].type = histogram::PROBABILITY;
  h[1].name = "Qp1(r)";   h[1].type = histogram::AVERAGE;
  h[2].name = "Qp2(r)";   h[2].type = histogram::AVERAGE;
  h[3].name = "Qc1(r)";   h[3].type = histogram::AVERAGE;
  h[4].name = "Dip_x*x";  h[4].type = histogram::AVERAGE;
  h[5].name = "Dip_y*y";  h[5].type = histogram::AVERAGE;
  h[6].name = "Dip_z*z";  h[6].type = histogram::AVERAGE;
  h[7].name = "U(ionion)"; h[7].type = histogram::AVERAGE;
  h[8].name = "U(iondip)"; h[8].type = histogram::AVERAGE;
  h[9].name = "U(dipdip)"; h[9].type = histogram::AVERAGE;
  h[10].name= "U(p1p2)"  ; h[10].type= histogram::AVERAGE;
  h[11].name= "Dip_P1"   ; h[11].type= histogram::AVERAGE;
  h[12].name= "Dip_P2"   ; h[12].type= histogram::AVERAGE;
  h[13].name = "Dip_x1";  h[13].type = histogram::AVERAGE;
  h[14].name = "Dip_y1";  h[14].type = histogram::AVERAGE;
  h[15].name = "Dip_z1";  h[15].type = histogram::AVERAGE;
  h[16].name = "Dip_x2";  h[16].type = histogram::AVERAGE;
  h[17].name = "Dip_y2";  h[17].type = histogram::AVERAGE;
  h[18].name = "Dip_z2";  h[18].type = histogram::AVERAGE;  
  h[19].name = "penalty"; h[19].type = histogram::PENALTY;
  h[20].name = "U(tot)";  h[20].type = histogram::AVERAGE;  
  h[21].name = "muex(salt)";  h[21].type = histogram::AVERAGE;
  h[22].name = "U_el(tot)";   h[22].type = histogram::AVERAGE;
  //h[8].name = "u(p1p2)";  h[8].type = histogram::AVERAGE;
  
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
      chain c("6 CATION CATION CATION CATION CATION CATION",graftpt);
      g[i]=s.insert_chain(c);
      g[CHAINS]+=g[i]; //convenient...
    };
    g[P1CHAINS] = g[CHAIN1] + g[CHAIN2] + g[CHAIN3] + g[CHAIN4] + g[CHAIN5] + g[CHAIN6];
    g[P2CHAINS] = g[CHAIN7] + g[CHAIN8] + g[CHAIN9] + g[CHAIN10] + g[CHAIN11] + g[CHAIN12];
  };
 
  g[SALT]  = s.insert_salt( nion3, chion3, rion3, cell); 
  g[SALT] += s.insert_salt( nion1, chion1, 2.0, cell, peptide::NA);
  g[SALT] += s.insert_salt( nion2, chion2, 2.0, cell, peptide::CL);
  g[SALT].name="mobile ions";
  
  s.load(".coord" + jobid); // load saved coordinates and charges.
  
  //PREPARE TITRATION
  titrate tit( pH );
  if (titrateBool==true) {
    tit.init( s.p, g[SALT] );
    s.trial = s.p;
    s.load(".coord" + jobid)==false; //load coords again (to get old charges)
  };
  
  //INITIAL ENERGY
  double utot = pot.energy(s.p) + pot.energy_vdw(s.p, g[P1], g[P2]);
  if (hairy==true)
	for (int k=CHAIN1; k<ENDCHAIN; k++)
		utot += pot.internal(s.p, g[k]);
  double umin = 1e6 ;                   //minimum system energy
  cout << "Internal = " << utot << endl;
  
  if (s.safetytest_vector()==false) return 1;
  cout << "Pass 3.\n";
 
  //PREPARE TRAJECTORY OUTPUT
  xyz trj("trj20-50");
  
  //PRINT INFORMATION
  cout << "\n# --- SYSTEM ----------------------------------------\n"
    << "# Steps               = " << mc.macroSteps << " x " << mc.microSteps
    << " = " << mc.macroSteps*mc.microSteps << endl
    << "# Cell radius         = " << cell_r << endl
    << "# Temperature (K)     = " << phys.T << endl
    << "# Jobid               = " << jobid << endl
    << "# Dielectric Const.   = " << phys.e_r << endl
    << "# Bjerrum length      = " << pot.lB << endl
    << "# vdW parameter       = " << pot.vdw << endl
    << "# Protein 1           = " << protein1 << endl
    << "# Protein 2           = " << protein2 << endl
    << "# Protein displ.      = " << prot_dp << endl
    << "# Protein rot.        = " << prot_rot <<endl
    << "# Max/min separation  = " << minsep << " " << maxsep << endl
    << "# Ion 1 (#  chg)      = " << nion1 << " " << chion1 << endl
    << "# Ion 2 (#  chg)      = " << nion2 << " " << chion2 << endl
    << "# Ion 3 (#  chg rad)  = " << nion3 << " " << chion3 << " " << rion3 << endl
    << "# System charge       = " << s.charge() << endl
    << "# Initial energy (kT) = " << utot << endl;
  if (titrateBool==true) {
    tit.info();
    cout << endl;
  };
  
  s.charge(g[P1]);
  s.charge(g[P2]);
  s.charge(g[SALT]);
  cout << "# --- GROUPS ----------------------------------------\n";
  cout << g[P1] 
    << "#   Mass center   = " << s.p[g[P1].cm] << endl
    << g[P2]
    << "#   Mass center   = " << s.p[g[P2].cm] << endl
    << g[SALT] << g[PROTEINS] << g[CHAINS] << endl;
  if (hairy==true)
    for (int k=CHAIN1; k<ENDCHAIN; k++)
      cout << g[k] << endl;

  //IMD Support
  imdwrap imd(s.p.size());
  if (imdBool==true) {
    pdb.load_particles(s.p);
    pdb.save("test.pqr");

    imd.wait_for_connection();
  };
  
  /************************
    Monte Carlo Simulation
    ************************/
  vector<particle> tmp;
  tmp = s.p;
  
  for (int macroCnt=1; macroCnt<=mc.macroSteps; macroCnt++) {
    for (int i=1; i<=mc.microSteps; i++) {
      bool rc,coll;
      montecarlo::rejectcause rejectcause;
      double u=0,uold=0,unew=0,du=0;
      
      //Move salt
      for (int j=0; j<g[SALT].size(); j++) {
        rc=false;
        coll=false;
        rejectcause=montecarlo::HC;
        int n=g[SALT].random();  //pick random particle
        s.displace(n, ion_dp);   //displace it...
        if (cell.cellCollision(s.trial[n])==false)
          if (col.overlap(s.trial, n)==false) {
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
        double m,z;
        z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
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
          hist.add(h[21], z, m); 
        };
      };
      
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
                  uold+=pot.energy(s.p, n);
                  unew+=pot.energy(s.trial, n);
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
        if (z <= maxsep && z >= minsep && z<2.*cell_r) {
          s.zmove( g[P1], dp );
          s.zmove( g[P2],-dp );
          moved=true;
        } else coll=true;
        
        if (coll==false)
          if (col.overlap(s.trial, g[PROTEINS])==false)
            if (col.overlap(s.trial, g[P1], g[P2])==false) {
              uold = pot.energy(s.p, g[PROTEINS]) + pot.energy(s.p, g[P1],g[P2]);
              unew = pot.energy(s.trial,g[PROTEINS]) + pot.energy(s.trial,g[P1],g[P2]);
              if (hairy==true) {
                for (int k=CHAIN1; k<ENDCHAIN; k++) {
                  uold += pot.graft(s.p, g[k]);
                  unew += pot.graft(s.trial, g[k]);
                };
              };
              du=unew-uold;
              //if (u_penalty!=0)
              //  du = du + hist.get(h[19],z)-hist.get(h[19],oldz);
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
        hist.add( h[19], z, u_penalty); //update penalty function
        hist.add( h[22],z,-0.37*pot.energy(s.p)+pot.energy_vdw(s.p, g[P1],g[P2]));
      };
      
      //Move proteins AND attached chains
      if (hairy==true && 1>2) {
        rc=false;
        coll=false;
        bool moved=false;
        rejectcause=montecarlo::HC;
        double oldz,z;
        double dp = clust_dp*slump.random_half();
        oldz = abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
        z    = oldz + 2.*dp;
        
        if (z <= maxsep && z >= minsep && z<2.*cell_r) { //check for cell overlap etc.
          s.zmove( g[P1], dp );
          s.zmove( g[P2],-dp );
          s.zmove( g[P1CHAINS], dp );
          s.zmove( g[P2CHAINS],-dp );
          moved=true;
        } else coll=true;
        
        if (coll==false)
          if (col.overlap(s.trial, g[P1CHAINS], g[P2CHAINS])==false)
            if (col.celloverlap(s.trial, g[CHAINS], cell_r)==false)
              if (col.overlap(s.trial, g[SALT])==false)
                if (col.overlap(s.trial, g[P1], g[P2CHAINS])==false)
                  if (col.overlap(s.trial, g[P2], g[P1CHAINS])==false)
                    if (col.overlap(s.trial, g[P1], g[P2])==false) {
                      uold = pot.energy(s.p, g[SALT])
                        + pot.energy(s.p, g[P1],g[P2])
                        + pot.energy(s.p, g[P1CHAINS], g[P2CHAINS])
                        + pot.energy(s.p, g[P1], g[P2CHAINS])
                        + pot.energy(s.p, g[P2], g[P1CHAINS]);
                      unew =  pot.energy(s.trial, g[SALT])
                        + pot.energy(s.trial, g[P1],g[P2])
                        + pot.energy(s.trial, g[P1CHAINS], g[P2CHAINS])
                        + pot.energy(s.trial, g[P1], g[P2CHAINS])
                        + pot.energy(s.trial, g[P2], g[P1CHAINS]);
                      du=unew-uold;
                      if (u_penalty!=0)
                        du = du + hist.get(h[19],z)-hist.get(h[19],oldz);
                      if (pot.metropolis(du)==true)
                        rc=true;
                      else
                        rejectcause=montecarlo::ENERGY;
                    };
        if (rc==true) {
          s.accept( g[PROTEINS] );
          s.accept( g[CHAINS] );
          mc.accept(montecarlo::CLUSTER);
          utot += du;
        } else {
          mc.reject(montecarlo::CLUSTER, rejectcause);
          if (moved==true) {
            s.undo( g[PROTEINS] );
            s.undo( g[CHAINS] );
          };
        };
        z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
        hist.add( h[0], z);             //update radial distribution function
        hist.add( h[19], z, u_penalty); //update penalty function
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
            if (col.overlap(s.trial, g[k], g[SALT])==false)
              if (col.overlap(s.trial, g[P1], g[P2])==false) {
                uold=pot.energy(s.p, g[k]) + pot.energy_vdw(s.p, g[P1], g[P2] );
                unew=pot.energy(s.trial, g[k]) + pot.energy_vdw(s.trial, g[P1], g[P2] );
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
      
      //Charge regulation
      if (titrateBool==true && slump.random_one()>0.7) {
        double z;
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
        tit.samplesites(s.p);
        z=abs(s.trial[g[P1].cm].z-s.trial[g[P2].cm].z);
        if (hairy==true)
          hist.add(h[3], z, s.charge(g[CHAIN1]));
      };
      
      //ANALYSIS
      double z = abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
      hist.add( h[20], z, utot ); //total energy 
 
      //Multipole analysis
          trj.header( g[P1].size()-2+1+g[SALT].size() ); //subtract ghosts, "mu", "cm". Add 1 for P2.cm
          trj.add( s.p, pep, g[P1] );
          trj.add( s.p, pep, g[P2]);
          trj.add( s.p, pep, g[SALT]);
          trj.footer();

      if (slump.random_one()>0.9) {
        double q1=s.charge(g[P1]);
        double q2=s.charge(g[P2]);
        Q1+=q1;
        Q2+=q2;
        Q1sq+=q1*q1;
        Q2sq+=q2*q2;
        point mu1=s.calcdipole(g[P1]);
        point mu2=s.calcdipole(g[P2]);
        dip1+=g[P1].dip;
        dip2+=g[P2].dip;
        hist.add(h[1], z, q1);
        hist.add(h[2], z, q2);
        hist.add(h[4], z, s.p[g[P1].dipv].x * s.p[g[P2].dipv].x);
        hist.add(h[5], z, s.p[g[P1].dipv].y * s.p[g[P2].dipv].y);
        hist.add(h[6], z, (s.p[g[P1].dipv].z-s.p[g[P1].cm].z)
            * (s.p[g[P2].dipv].z-s.p[g[P2].cm].z) );
        hist.add(h[7], z, pot.lB*q1*q2/z);
        hist.add(h[8], z, pot.iondip(mu1,q2,z)+pot.iondip(mu2,q1,z));
        hist.add(h[9], z, pot.dipdip(mu1,mu2,z));

        //single dipole orientation
        hist.add(h[13], z, s.p[g[P1].dipv].x );
        hist.add(h[14], z, s.p[g[P1].dipv].y );
        hist.add(h[15], z, s.p[g[P1].dipv].z-s.p[g[P1].cm].z);
        hist.add(h[16], z, s.p[g[P2].dipv].x );
        hist.add(h[17], z, s.p[g[P2].dipv].y );
        hist.add(h[18], z, s.p[g[P2].dipv].z-s.p[g[P2].cm].z);        
        mu1_z += s.p[g[P1].dipv].z-s.p[g[P1].cm].z;
        mu2_z += s.p[g[P2].dipv].z-s.p[g[P2].cm].z;

        bool savevdw=g[1].vdw; g[1].vdw=false; //turn off vdw interactions
        double up1p2=pot.energy(s.p,g[P1],g[P2]);
        hist.add(h[10],z, up1p2);
        g[1].vdw=savevdw;

        hist.add(h[11],z, g[P1].dip);
        hist.add(h[12],z, g[P2].dip);

        //Minimum (electrostatic) energy snapshot
        if (up1p2<umin && minsnapshot==true) {
          umin=up1p2;
          povray pov; 
          pov.sphere(cell_r); 
          pov.zaxis(cell_r); 
          pov.add( s.p, g[P1] ); 
          pov.add( s.p, g[P2] ); 
          pov.save( ".min.pov" + jobid );
          s.save( ".min.coord" + jobid );
        };
      };

      //Adjust displacement parameter
      //(Detailed balance violation! Use for equilibration, only)
      if (adjust_dp==true && slump.random_one()>0.9) {
        mc.adjust_dp(montecarlo::ION, ion_dp);
        mc.adjust_dp(montecarlo::TRANSLATE, prot_dp);
        mc.adjust_dp(montecarlo::ROTATE, prot_rot);
        mc.adjust_dp(montecarlo::MONOMER, dp_monomer);
        mc.adjust_dp(montecarlo::CLUSTER, clust_dp);
      };    

      if (imdBool==true && slump.random_one()>0.9)
        imd.send_particles(s.p);

    }; //end of micro loop

    // Check if particle vectors are in sync
    if (s.safetytest_vector()==false)
      return 1; //teminate program

    // System energy check:
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy(s.p) + pot.energy_vdw(s.p, g[P1], g[P2]);
    if (hairy==true)
      for (int k=CHAIN1; k<ENDCHAIN; k++)
        usys += pot.internal(s.p, g[k]);

    double mu   = -log(muex.avg());
    cout << "# System energy (kT)       = " << usys << endl
      << "# Energy drift (rel abs)      = " << utot - usys << " " << (utot-usys)/usys << endl
      << "# DP (ion,prot,rot,cha,clu)= " << ion_dp<<" "<<prot_dp<<" "<<prot_rot<<" "<<dp_monomer<<" "<<clust_dp<<endl
      << "# Protein charges          = " << s.charge(g[P1]) << " " << s.charge(g[P2]) << endl
      << "# Total system charge      = " << s.charge() << endl
      << "# Widom analysis (1:1 salt): " << endl
                                            << "#   excess. chem. pot.     = " << mu << endl
                                            << "#   mean activity coeff.   = " << exp(mu) << endl
                                            << "#   debye length           = " << -pot.lB/(2.0*mu) << endl
                                            << "# Protein 1 (Q,C,mu)       = "
                                            << Q1.avg() << " " << Q1sq.avg()-pow(Q1.avg(),2) << " " << dip1.avg() << endl 
                                            << "# Protein 2 (Q,C,mu)       = "
                                            << Q2.avg() << " " << Q2sq.avg()-pow(Q2.avg(),2) << " " << dip2.avg() << endl
                                            << "# Avg. dipole z components = " << mu1_z << " " << mu2_z << endl
                                            << endl;

    s.save( ".coord" + jobid);          //save coordinates
    trj.newfile();                      //dump trj into new file

    //Graphical output (Povray and vpython)
    povray pov;
    pov.sphere(cell_r);
    pov.zaxis(cell_r);
    pov.add( s.p, g[P1] );
    pov.add( s.p, g[P2] );
    pov.add( s.p, g[SALT] );
    if (hairy==true)
      for (int i=CHAIN1; i<ENDCHAIN; i++)
        pov.add( s.p, g[i] );
    pov.save("snapshot.pov" + jobid);

    /*
       vpython vpy;
       vpy.add( s.p, g[SALT]);
       vpy.add( s.p, g[P1] );
       vpy.add( s.p, g[P2] );                                           
       if (hairy==true)
       for (int i=CHAIN1; i<ENDCHAIN; i++)
       vpy.add( s.p, g[i] );                                                                  
       vpy.save("vpython" + jobid + ".py");
       */

    //test for any static particles...
    for (int i=0; i<tmp.size(); i++)
      if (abs(tmp[i].x-s.p[i].x)<0.01)
        if (abs(tmp[i].y-s.p[i].y)<0.01)
          if (abs(tmp[i].z-s.p[i].z)<0.01)
            cout << "Possible static particle: " << i << endl;

    cout << endl;
    //hist.show(h, minsep, maxsep );

  }; //end of macro step 

  /**************************************
    Print results and terminate program
   **************************************/
  trj.close();
  hist.show(h, minsep, maxsep);
  mc.showStatistics();

  if (smear==true)
    tit.applycharges(s.p);    // apply average charges, and smeared out protons
  s.save(".coord" + jobid);

  // Save a PQR file for VMD.
  pqr pqr("out.pqr");
  pqr.add(s.p,pep,g[P1]);
  pqr.add(s.p,pep,g[P2]);
  pqr.add(s.p,pep,g[SALT]);
  pqr.close();

  for (int i=g[CHAIN1].beg; i<g[CHAIN1].end; i++) {
    cout << geo.dist( s.p[i], s.p[i+1]) << endl;
  };

  return 0;
};
