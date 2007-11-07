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
  double maxsep   = par.getDbl("maxsep",cell_r/1.5);//restrict protein separation (max)
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
  
  int nion1       = par.getInt("nion1");  //number of ion 1
  int nion2       = par.getInt("nion2");  //number of ion 2
  int nion3       = par.getInt("nion3");  //number 
  int chion1      = par.getInt("chion1"); //ionic charges
  int chion2      = par.getInt("chion2");
  int chion3      = par.getInt("chion3");

  bool hairy       = par.getBool("hairy", true);
  bool rotation    = par.getBool("rotate", true);
  bool titrateBool = par.getBool("titrate", true);
  bool smear       = par.getBool("smear", false);
  bool minsnapshot = par.getBool("minsnapshot",false);
  bool adjust_dp   = par.getBool("adjust_dp", false);

  /**************************************
    Class constructors and declarations
  **************************************/
  enum groupid {P1=0,P1CHG,P2,P2CHG,SALT,CHAIN1,ENDCHAIN,PROTEINS}; //convenient
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
  vector<histogram::data> h(20);
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

  if (hairy==true) {
    chain ch1("2 UNK UNK",0);
    g[CHAIN1]=s.insert_chain( ch1 );
    g[CHAIN1].name="peptide chain 1";
  };

  g[SALT]  = s.insert_salt( nion1, chion1, 2.0, cell);
  g[SALT] += s.insert_salt( nion2, chion2, 2.0, cell);
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
  double utot = pot.energy(s.p);        //system energy
  double umin = 1e6 ;                   //minimum system energy

  if (s.safetytest_vector()==false) return 1;
  cout << "Pass 3.\n";

  //PRINT INFORMATION
  cout << "\n# --- SYSTEM ----------------------------------------\n"
       << "# Macro loops         = " << mc.macroSteps << endl
       << "# Micro loops         = " << mc.microSteps << endl
       << "# Total nr. of loops  = " << mc.macroSteps*mc.microSteps << endl
       << "# Cell radius         = " << cell_r << endl
       << "# Temperature (K)     = " << phys.T << endl
       << "# Dielectric Const.   = " << phys.e_r << endl
       << "# Bjerrum length      = " << pot.lB << endl
       << "# vdW parameter       = " << pot.vdw << endl
       << "# Protein 1           = " << protein1 << endl
       << "# Protein 2           = " << protein2 << endl
       << "# Protein displ.      = " << prot_dp << endl
       << "# Protein rot.        = " << prot_rot <<endl
       << "# Max protein sep.    = " << maxsep << endl
       << "# Min protein sep.    = " << minsep << endl
       << "# Ion 1 (#  chg)      = " << nion1 << " " << chion1 << endl
       << "# Ion 2 (#  chg)      = " << nion2 << " " << chion2 << endl
       << "# Ion 3 (#  chg)      = " << nion3 << " " << chion3 << endl
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
       << g[CHAIN1] << g[SALT] << g[PROTEINS]
       << endl;

  /************************
    Monte Carlo Simulation
  ************************/
  for (int macroCnt=1; macroCnt<=mc.macroSteps; macroCnt++) {
    for (int i=1; i<=mc.microSteps; i++) {
      bool rc,coll;
      montecarlo::rejectcause rejectcause;
      double u=0,uold=0,unew=0,du;

      //Move salt
      for (int j=0; j<g[SALT].size(); j++) {
	rc=false;
	rejectcause=montecarlo::HC;
	int n=g[SALT].random();  //pick random particle
	s.displace(n, ion_dp);   //displace it...
	if (cell.cellCollision(s.trial[n])==false)
	  if (col.overlap(s.trial, g[P1], n)==false)
	    if (col.overlap(s.trial, g[P2], n)==false)
	      if (col.overlap(s.trial, g[SALT], n)==false) {
		//if (col.overlap(s.trial, g[CHAIN1], n)==false) {
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
      if (slump.random_one()>0.5) {
        particle anion, cation;
        anion.charge=-1;
        cation.charge=+1;
        anion.radius=cation.radius=2.0;
        for (int k=0; k<100; k++) {
          cell.randomPos(cation);
          cell.randomPos(anion);
          u=1.0e20;
          if (col.overlap(s.p, anion)==false)
            if (col.overlap(s.p, cation)==false)
              u=pot.potential(s.p, cation) -
                pot.potential(s.p, anion);
          muex += exp(-u);
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
	    if (col.overlap(s.trial, n)==false) {
	      uold=pot.chain(s.p, g[j], n);
	      unew=pot.chain(s.trial, g[j], n);
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
	  };
	};
      };

      //Move proteins
      rc=false;
      coll=false;
      rejectcause=montecarlo::HC;
      double oldz,z;
      double dp = prot_dp*slump.random_half();
      oldz = abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
      z    = oldz + 2.*dp;
      if (z <= maxsep && z >= minsep) {
        s.zmove( g[P1], dp );
        s.zmove( g[P2],-dp );
        if (z < 2*cell_r)
          if (col.overlap(s.trial, g[PROTEINS])==false)
            if (col.overlap(s.trial, g[P1], g[P2])==false) {
              uold = pot.energy(s.p, g[PROTEINS]) + pot.energy(s.p, g[P1],g[P2]);
              unew = pot.energy(s.trial,g[PROTEINS]) + pot.energy(s.trial,g[P1],g[P2]);
              if (hairy==true) {
                uold += pot.graft(s.p, g[CHAIN1]);
                unew += pot.graft(s.trial, g[CHAIN1]);
              };
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
          mc.accept(montecarlo::TRANSLATE);
          utot += du;
        } else {
          s.undo( g[PROTEINS] );
          mc.reject(montecarlo::TRANSLATE, rejectcause);
        };
      };
      z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
      hist.add( h[0], z);             //update radial distribution function
      hist.add( h[19], z, u_penalty); //update penalty function

      //Rotate protein 1
      if (rotation==true && slump.random_one()>0.5) {
        rc=false;
        coll=false;
        rejectcause=montecarlo::HC;
        s.rotate( g[P1], prot_rot );

        if (hairy==true)
          for (int i=CHAIN1; i<ENDCHAIN; i++) { // chain collision...
            coll=col.overlap(s.trial, g[P1], g[i]);
            if (coll==true) break;
          };
        if (coll==false)
          if (col.overlap(s.trial, g[P1], g[SALT])==false)
            if (col.overlap(s.trial, g[P1], g[P2])==false) {
              uold=pot.energy(s.p, g[P1]);
              unew=pot.energy(s.trial, g[P1]);
              if (hairy==true) {
                uold += pot.graft(s.p, g[CHAIN1]);
                unew += pot.graft(s.trial, g[CHAIN1]);
              };
              du = unew-uold;
              if (pot.metropolis(du))
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
        //hist.add( h[0], z, )
      };

      //Rotate protein 2
      if (rotation==true && slump.random_one()>0.5) {
        rc=false;
        coll=false;
        rejectcause=montecarlo::HC;
        s.rotate( g[P2], prot_rot );
        if (hairy==true)
          for (int i=CHAIN1; i<ENDCHAIN; i++) { // chain collision...
            coll=col.overlap(s.trial, g[P2], g[i]);
            if (coll==true) break;
          };
        if (coll==false)
          if (col.overlap(s.trial, g[P2], g[SALT])==false)
            if (col.overlap(s.trial, g[P1], g[P2])==false) {
              uold=pot.energy(s.p, g[P2]);
              unew=pot.energy(s.trial, g[P2]);
              if (hairy==true) {
                uold += pot.graft(s.p, g[CHAIN1]);
                unew += pot.graft(s.trial, g[CHAIN1]);
              };
              du = unew-uold;
              if (pot.metropolis(du))
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

      //Charge regulation
      if (titrateBool==true && slump.random_one()>0.7) {
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

      //Multipole analysis
      if (slump.random_one()>0.6) {
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
        z=abs(s.p[g[P1].cm].z-s.p[g[P2].cm].z);
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
      if (adjust_dp==true && slump.random_one()>0.9) {
        ion_dp  += mc.adjust_dp(montecarlo::ION);
        prot_dp += mc.adjust_dp(montecarlo::TRANSLATE);
        prot_rot+= mc.adjust_dp(montecarlo::ROTATE);
        if (ion_dp<0) ion_dp=-ion_dp;
        if (prot_dp<0) prot_dp=-prot_dp;
        if (prot_rot<0) prot_rot=-prot_rot;
        if (prot_rot>3.14) prot_rot=3;
      };

      
    }; //end of micro loop

    // Check if particle vectors are in sync
    if (s.safetytest_vector()==false)
      return 1; //teminate program

    // System energy check:
    cout << "# --- MACROSTEP COMPLETED  --------------------------\n";
    mc.showStatus(macroCnt);
    double usys = pot.energy(s.p);
    double mu   = -log(muex.avg());
    cout << "# System energy (kT)       = " << usys << endl
	 << "# System energy drift (kT) = " << utot - usys << endl
         << "# DP (ion,prot,rot)        = " << ion_dp << " " << prot_dp << " "<< prot_rot << endl
	 << "# Protein charges          = " << s.charge(g[P1]) << " " << s.charge(g[P2]) << endl
	 << "# Total system charge      = " << s.charge() << endl
	 << "# Salt (1:1) Excess pot.   = " << mu << endl
	 << "# Salt activity coeff.     = " << muex << endl
	 << "# Debye length             = " << -pot.lB/2.0/mu << endl
	 << "# Protein 1 (Q,C,mu)       = "
	 << Q1.avg() << " " << Q1sq.avg()-pow(Q1.avg(),2) << " " << dip1.avg() << endl 
	 << "# Protein 2 (Q,C,mu)       = "
	 << Q2.avg() << " " << Q2sq.avg()-pow(Q2.avg(),2) << " " << dip2.avg() << endl
         << "# Avg. dipole z components = " << mu1_z << " " << mu2_z << endl
	 << endl;

    s.save( ".coord" + jobid); //save coordinates

    //Graphical output (Povray and vpython)
    povray pov;
    pov.sphere(cell_r);
    pov.zaxis(cell_r);
    pov.add( s.p, g[P1] );
    pov.add( s.p, g[P2] );
    pov.add( s.p, g[SALT] );
    if (hairy==true)
      pov.add( s.p, g[CHAIN1] );
    pov.save("snapshot.pov" + jobid);

    vpython vpy;
    vpy.add( s.p, g[SALT]);
    vpy.add( s.p, g[P1] );
    vpy.add( s.p, g[P2] );

    //vpy.add( s.p, g[P2] );
    //vpy.add( s.p, g[SALT] );                                                                      
    if (hairy==true)
      vpy.add( s.p, g[CHAIN1] );                                                                  
    vpy.save("vpython" + jobid + ".py");

  }; //end of macro step 

  /**************************************
      Print results and terminate program
  **************************************/
  hist.show(h, 0, cell_r);
  mc.showStatistics();

  if (smear==true)
    tit.applycharges(s.p);    // apply average charges, and smeared out protons
  s.save(".coord" + jobid);
    
  for (int i=g[CHAIN1].beg; i<g[CHAIN1].end; i++) {
    cout << geo.dist( s.p[i], s.p[i+1]) << endl;
  };

  return 0;

};
