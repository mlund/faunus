/*! \page test_wallprotein Multiple wallproteins with Debye Huckel interactions
 *
 * Simulate a number of flexible and titratable
 * polymers in a salt solution in the presence
 * of a charged Hydropobic Gouy-Chapman surface 
 * with implicit salt interactions given by Debye Huckel.
 *
 * \author Chris Evers
 * \date Lund, 2011
 * \include wpDH_mult.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/energy/externalpotential.h"
#include "faunus/potentials/pot_r12debyehuckel.h"
#include "faunus/energy/penaltyfunction.h"

using namespace std;
using namespace Faunus;

int main() {
  // General
  cout << faunus_splash();
  inputfile in("wp.conf");
  slump slp;                                // load random number generator
  mcloop loop(in);                          // define simulation loops
  canonical nvt;                            // canonical ensemble

  // Simulation container
#ifdef NOSLIT                               // If wpDH_mult_noslit
  cuboid con(in);
  typedef pot_r12debyehuckel_hydrophobic_tab Tpot;
#elif HYDROPHIC                             // If wpDH_mult_hydrophobic
  cuboidslit con(in);
  typedef pot_r12debyehuckelXY_hydrophobic_tab Tpot;
  typedef expot_hydrophobic Texpot;
#else                                       // If wpDH_mult
  cuboidslit con(in);
  typedef pot_r12debyehuckelXY_hydrophobic_tab Tpot;
  typedef expot_gchydrophobic Texpot;
#endif

  // Polymer
  vector< polymer > pol;
  pol.resize(in.getint("pol_num",1));
  for (int i=0; i<pol.size(); i++) {
    pol[i].babeladd( con, in, i );          //  add from input
    con.trial=con.p;                     //  synchronize particle vector
    pol[i].masscenter(con);              //  update masscenter
    pol[i].move(con, -pol[i].cm + con.randompos() );    // translate to random position
    while ( con.slicecollision(pol[i].cm_trial) )
      pol[i].move(con, -pol[i].cm + con.randompos() );  // translate to random position
    pol[i].accept(con);                  //  accept translation
    if (in.getflt("eqtit_runfrac",0.5)<1e-3) {
      pol[i].loadCharges(in.getstr("pol_charges","q.in"), con.p); // load residue charges from file
      con.trial=con.p;
      cout << "! Loaded polymer charges from disk" << endl;
    }
    cout << pol[i].info();
  }

  // Potentials
#ifdef NOSLIT
  springinteraction<Tpot> pot(in); 
#else
  springinteraction_expot<Tpot, Texpot> pot(in);
  pot.expot.update(con);
#endif

  // Distribution functions
  double* zhalfPtr=&con.len_half.z;         // half length of container in z-direction
  distributions2 dst(0.1,0,(*zhalfPtr)*2.); // distance dependent averages
  distributions pairdst;                    // distance dependent averages
  histogram gofr(0.1,0, (*zhalfPtr)*2. );   // radial distribution function
//   correlation<double> z_cor(200);
  histogram qdist(1,-pol[0].nb.size(), pol[0].nb.size()); // polymer charge distribution function
//   correlation<double> q_cor(100);
  histogram rgdist(.2,0, (*zhalfPtr)*2. );  // radius of gyration distribution function
//   correlation<double> rg_cor(100);
  histogram reedist(.2,0, (*zhalfPtr)*2. ); // end to end distance distribution function
#ifdef NOSLIT
  histogram pairgofr(.2,0, .5*sqrt(3)*(*zhalfPtr)*2. );   // radial distribution function
#else
  histogram pairgofr(.2,0, .5*sqrt(3)*(*zhalfPtr)*2. );   // radial distribution function
//   histogram2D pairgofr2D(.2,0, .5*sqrt(3)*(*zhalfPtr)*2. );   // radial distribution function
#endif

  // Moves
  monomermove mm(nvt,con,pot,in);
  crankShaft cs(nvt,con,pot,in);
  branchRotation br(nvt,con,pot,in);
  macrorot mr(nvt,con,pot,in);
  translate mt(nvt,con,pot,in);
  if ( pol.size() == 1 )
    mt.dpv.x=mt.dpv.y=0;                          // no need to translate in xy direction
  eqtitrate tit(nvt,con,pot,in,"eqtit_");

  // Particle input output
  io io;
  iopqr pqr;
  ioxtc xtc(con.len.z);
  xtc.setbox(con.len.x,con.len.y,con.len.z);
  ioaam aam;
  aam.load(con, "confout.aam");     //   load stored configuration
  for (int i=0; i<pol.size(); i++) 
    pol[i].masscenter(con);              //   update masscenter
  aam.save("confin.aam", con.p);

  // Initial system energy
  double utot=tit.intrinsicenergy(con.p);
  for (int i=0; i<pol.size(); i++)
  utot += pot.energy( con.p, pol[i]) + 
          pot.uself_polymer(con.p, pol[i]);
  if ( pol.size()>1 )
    for (int i=0; i<pol.size()-1; i++)
      for (int j=i+1; j<pol.size(); j++)
        utot -= pot.energy( con.p, pol[i], pol[j] );  // double counting
  systemenergy sys(utot);

  cout << con.info() << atom.info()
    << pot.info() << tit.info() 
    << in.info();

  // Simulation loop
  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      int i (rand() % pol.size());
      switch (rand() % 6) {     // random number between 0 and 5
        case 0:
          for (int j=0; j<pol[i].size(); j++)
            sys+=mm.move(pol[i]);   // move monomers 
          pol[i].masscenter(con);
          gofr.add(*zhalfPtr-pol[i].cm.z);
          break;
        case 1:
          sys+=mt.move(pol[i]);     // translate polymer
          gofr.add(*zhalfPtr-pol[i].cm.z);
          break;
        case 2:
          pol[i].masscenter(con);
          sys+=mr.move(pol[i]);     // rotate polymer
          break;
        case 3:
          sys+=cs.move(pol[i]);     // crankshaft
          gofr.add(*zhalfPtr-pol[i].cm.z);
          break;
        case 4:
          sys+=br.penaltymove(pol[i]);     // branch rotation
          gofr.add(*zhalfPtr-pol[i].cm.z);
          break;
        case 5:
          sys+=tit.move();       // titrate titratable sites
          pol[i].charge(con.p);
          break;
      }

      if (slp.random_one()>0.3) {
        int i (rand() % pol.size());
        pol[i].masscenter(con);
        double z=*zhalfPtr-pol[i].cm.z;  // distance between masscenter and wall
//         z_cor += z;

        double charge=pol[i].charge(con.p);
        dst.add("qtot", z, charge);
        qdist.add(charge);
//         q_cor += charge;

        double rg2=pol[i].sqmassgradius(con);
        double rg=sqrt(rg2);
        dst.add("Rg2", z, rg2);
        dst.add("Rg", z, rg);
        rgdist.add(rg);
//         rg_cor += sqrt(rg2);

        double ree2=pol[i].sqend2enddistance(con);
        double ree=sqrt(ree2);
        dst.add("Ree2", z, ree2); 
        dst.add("Ree", z, ree); 
        reedist.add(ree);

        if ( pol.size() > 1 ) {
          int j (rand() % pol.size());
          while ( i==j )
            j = (rand() % pol.size());
          pol[j].masscenter(con);
          double r12=con.dist(pol[i].cm,pol[j].cm); // distance between masscenters
          pairgofr.add(r12);
          pairdst.add("qtot", r12, charge);
          pairdst.add("Rg2", r12, rg2);
          pairdst.add("Rg", r12, rg);
          pairdst.add("Ree2", r12, ree2);
          pairdst.add("Ree", r12, ree);
        }

        for (int j=pol[i].beg; j<=pol[i].end; j++) {
          std::ostringstream s; s << "q" << j-pol[i].beg;
          dst.add( s.str(), *zhalfPtr-con.p[i].z, con.p[i].charge );    
        }
      }

      if (slp.random_one()>0.95) {
        xtc.save( "traj.xtc", con.p );
      }
    } // END of micro loop

    utot=tit.intrinsicenergy(con.p);
    for (int i=0; i<pol.size(); i++)
      utot += pot.energy( con.p, pol[i]) + 
              pot.uself_polymer(con.p, pol[i]);
    if ( pol.size()>1 )
      for (int i=0; i<pol.size()-1; i++)
        for (int j=i+1; j<pol.size(); j++)
          utot -= pot.energy( con.p, pol[i], pol[j] );
    sys.update(utot);

    cout << loop.timing() << "#   Energy drift = " << sys.cur-sys.sum << " kT. "
      << "System charge = " << con.charge() << ". " << endl;

    // Write files to disk
    aam.save("confout.aam", con.p);
    if ( in.getboo("write_files",true) ) {
      io.writefile("vmdbonds.tcl", pol[0].getVMDBondScript());
      pqr.save("confout.pqr",con.p);
      dst.write("dist.out");
      pairdst.write("pairdist.out");
      dst.cntwrite("cntdist.out");
      gofr.write("gofr.out");
      gofr.dump("gofr.xy");
      pairgofr.dump("pairgofr.xy");
//       z_cor.dump("z_cor.xy");

      qdist.dump("q.xy");
//       q_cor.dump("q_cor.xy");
      rgdist.dump("rg.xy");
//       rg_cor.dump("rg_cor.xy");
      reedist.dump("ree.xy");

//       pot.pair.dump("pairpot.xy");
#ifndef NOSLIT
      pot.expot.dump("expot.xy");
#endif

      tit.applycharges(con.trial);                      // Set average charges on all titratable sites
      pol[0].saveCharges("q.out", con.trial);              // Save average charges to disk
      con.trial=con.p;                                  // Restore original charges
    }

  } // END of macro loop and simulation

cout << sys.info() << loop.info() << mm.info()
  << cs.info() << br.info() << mr.info() 
  << mt.info() << tit.info(con.p);
for (int i=0; i<pol.size(); i++)
  cout << pol[i].info();
cout << pot.info();

}
