#include <iostream>
#include "group.h"
#include "point.h"
#include "space.h"
#include "slump.h"
#include "potentials.h"
#include "species.h"
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
#include "io.h"

/*
 * Monte Carlo simulation of Protein-Protein Interactions
 * Mikael Lund, 2003-2007.
 */

using namespace std;

int main(int argc, char* argv[] ) {
  
  /**********************************
  Handle arguments and input file
  **********************************/
  string cfgfile = "amino.inp";
  if (argc==2)
    cfgfile = argv[1];
  config c(cfgfile);  
  
  /**************************************
    Class constructors and declarations
    **************************************/
  enum groupid {P1=0,P2,SALT,PROTEINS,LAST};
  montecarlo mc( c.macro, c.micro );
  physconst phys(c.temp, c.dielec);
  pdbclass pdb;
  cell cell(c.cell_r);
  slump slump;
  space s;
  species spc;
  coulombvdw pot;
  hardsphere col;
  io io;
  vector<group> g(LAST+1);
  
  average<double> Q1,Q2, dip1,dip2,Q1sq, Q2sq;
  average<double> muex;                 //salt excess chemical potential (1:1)
  average<double> mu1_z, mu2_z;         //dipole z components
  
  pot.vdw = c.vdw;
  pot.k   = c.springk;
  pot.r0  = c.springeq;
  
  /**************************************
    Prepare distribution functions
    **************************************/
  vector<histogram::data> h(22);
  histogram hist( 0.5, 2*c.cell_r );
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
  //h[8].name = "u(p1p2)";  h[8].type = histogram::AVERAGE;
  
  /**********************************
    Load proteins, salt and polymers
    ***********************************/
  g[P1] = s.append( io.loadaam(spc,c.protein1) );
  g[P1].name = "protein 1";
  g[P1].vdw = true;
  s.addmasscenter(g[P1]);
  s.adddipole(g[P1]);
  s.move( g[P1], -s.center_of_mass(g[P1]), space::AUTOACCEPT);
  s.zmove( g[P1],  (c.maxsep/2.-10), space::AUTOACCEPT);
  g[P1].radius=s.radius( g[P1], s.p[g[P1].cm] );
  
  g[P2] = s.append( io.loadaam(spc,c.protein2) );
  g[P2].name = "protein 2";
  g[P2].vdw = true;
  s.addmasscenter( g[P2] );
  s.adddipole(g[P2]);
  s.move( g[P2], -s.center_of_mass(g[P2]), space::AUTOACCEPT);
  s.zmove( g[P2], -(c.maxsep/2.-10), space::AUTOACCEPT);
  g[P2].radius=s.radius( g[P2], s.p[g[P2].cm] );
  
  g[PROTEINS] = g[P2] + g[P1]; //convenient
  
  g[SALT]  = s.insert_salt( c.nion3, c.zion3, c.rion3, cell); 
  g[SALT] += s.insert_salt( c.nion1, c.zion1, 2.0, cell);
  g[SALT] += s.insert_salt( c.nion2, c.zion2, 2.0, cell);
  g[SALT].name="mobile ions";
  
  s.load(".coord" + c.jobid); // load saved coordinates and charges.
  
  //PREPARE TITRATION
  titrate tit( spc, c.pH );
  if (c.titrateBool==true) {
    tit.init( s.p, g[SALT] );
    s.trial = s.p;
    s.load(".coord" + c.jobid)==false; //load coords again (to get old charges)
  };
  
  //INITIAL ENERGY
  double utot = pot.energy(s.p) + pot.energy_vdw(s.p, g[P1], g[P2]);
  double umin = 1e6 ;                   //minimum system energy
  cout << "Internal = " << utot << endl;
  
  if (s.safetytest_vector()==false) return 1;
  cout << "Pass 3.\n";
  
  //PRINT INFORMATION
  cout << "\n# --- SYSTEM ----------------------------------------\n"
    << "# Steps               = " << mc.macroSteps << " x " << mc.microSteps
    << " = " << mc.macroSteps*mc.microSteps << endl
    << "# Cell radius         = " << c.cell_r << endl
    << "# Temperature (K)     = " << phys.T << endl
    << "# Jobid               = " << c.jobid << endl
    << "# Dielectric Const.   = " << phys.e_r << endl
    << "# Bjerrum length      = " << pot.lB << endl
    << "# vdW parameter       = " << pot.vdw << endl
    << "# Protein 1           = " << c.protein1 << endl
    << "# Protein 2           = " << c.protein2 << endl
    << "# Protein displ.      = " << c.prot_dp << endl
    << "# Protein rot.        = " << c.prot_rot <<endl
    << "# Max/min separation  = " << c.minsep << " " << c.maxsep << endl
    << "# Ion 1 (#  chg)      = " << c.nion1 << " " << c.zion1 << endl
    << "# Ion 2 (#  chg)      = " << c.nion2 << " " << c.zion2 << endl
    << "# Ion 3 (#  chg rad)  = " << c.nion3 << " " << c.zion3 << " " << c.rion3 << endl
    << "# System charge       = " << s.charge() << endl
    << "# Initial energy (kT) = " << utot << endl;
  if (c.titrateBool==true) {
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
    << g[SALT] << g[PROTEINS] << endl;

  //IMD Support
  imdwrap imd(s.p.size());
  if (c.imdBool==true) {
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
        s.displace(n, c.ion_dp);   //displace it...
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
                  + pot.pot(s.p, anion) * anion.charge \
                  + pot.pot(s.p, cation)* cation.charge;
          m=exp(-u);
          muex += m;
          hist.add(h[21], z, m); 
        };
      };
      
      //Move proteins
      if (slump.random_one()>0.2) {
        rc=false;
        coll=false;
        bool moved=false;
        rejectcause=montecarlo::HC;
        double oldz,z;
        double dp = c.prot_dp*slump.random_half();
        oldz = abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
        z    = oldz + 2.*dp;
        if (z < c.maxsep && z > c.minsep && z<2.*c.cell_r) {
          s.zmove( g[P1], dp );
          s.zmove( g[P2],-dp );
          moved=true;
        } else coll=true;
        
        if (coll==false)
          if (col.overlap(s.trial, g[PROTEINS])==false)
            if (col.overlap(s.trial, g[P1], g[P2])==false) {
              uold = pot.energy(s.p, g[PROTEINS]) + pot.energy(s.p, g[P1],g[P2]);
              unew = pot.energy(s.trial,g[PROTEINS]) + pot.energy(s.trial,g[P1],g[P2]);
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
        hist.add( h[19], z, c.u_penalty); //update penalty function
      };
      
      //Rotate proteins
      for (int k=P1; k<=P2; k++) {
        double z;
        if (c.rotate==true && slump.random_one()>0.5) {
          rc=false;
          coll=false;
          rejectcause=montecarlo::HC;
          s.rotate( g[k], c.prot_rot );
          if (coll==false)
            if (col.overlap(s.trial, g[k], g[SALT])==false)
              if (col.overlap(s.trial, g[P1], g[P2])==false) {
                uold=pot.energy(s.p, g[k]) + pot.energy_vdw(s.p, g[P1], g[P2] );
                unew=pot.energy(s.trial, g[k]) + pot.energy_vdw(s.trial, g[P1], g[P2] );
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
      if (c.titrateBool==true && slump.random_one()>0.7) {
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
      };
      
      //ANALYSIS
      double z = abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
      hist.add( h[20], z, utot ); //total energy 
      //Multipole analysis
      if (slump.random_one()>0.5) {
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
        if (up1p2<umin && c.minsnapshot==true) {
          umin=up1p2;
          povray pov; 
          pov.sphere(c.cell_r); 
          pov.zaxis(c.cell_r); 
          pov.add( s.p, g[P1] ); 
          pov.add( s.p, g[P2] ); 
          pov.save( ".min.pov" + c.jobid );
          s.save( ".min.coord" + c.jobid );
        };
      };
      
      //Adjust displacement parameter
      if (c.adjust_dp==true && slump.random_one()>0.9) {
        mc.adjust_dp(montecarlo::ION, c.ion_dp);
        mc.adjust_dp(montecarlo::TRANSLATE, c.prot_dp);
        mc.adjust_dp(montecarlo::ROTATE, c.prot_rot);
      };    

      if (c.imdBool==true && slump.random_one()>0.9)
        imd.send_particles(s.p);

    }; //end of micro loop
    
    // Check if particle vectors are in sync
    if (s.safetytest_vector()==false)
      return 1; //teminate program
    
    // System energy check:
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy(s.p) + pot.energy_vdw(s.p, g[P1], g[P2]);
    
    double mu   = -log(muex.avg());
    cout << "# System energy (kT)       = " << usys << endl
      << "# Energy drift (rel abs)      = " << utot - usys << " " << (utot-usys)/usys << endl
      << "# DP (ion,prot,rot)        = " << c.ion_dp<<" "<<c.prot_dp<<" "<<c.prot_rot<<" "<<endl
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
      << "# Avg. dipole z components = " << mu1_z.avg() << " " << mu2_z.avg() << endl
      << endl;
    
    s.save( ".coord" + c.jobid); //save coordinates
    
    //Graphical output (Povray and vpython)
    povray pov;
    pov.sphere(c.cell_r);
    pov.zaxis(c.cell_r);
    pov.add( s.p, g[P1] );
    pov.add( s.p, g[P2] );
    pov.add( s.p, g[SALT] );
    pov.save("snapshot.pov" + c.jobid);
    
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
  hist.show(h, c.minsep, c.maxsep);
  mc.showStatistics();
  
  if (c.smear==true)
    tit.applycharges(s.p);    // apply average charges, and smeared out protons
  s.save(".coord" + c.jobid);
  
  return 0;
};
