/*! \page test_wallprotein WallproteinDH
 *
 * Simulate a number of flexible and titratable
 * polymers in a salt solution in the presence
 * of a charged Gouy-Chapman surface with implicit
 * salt interactions given by Debye Huckel.
 *
 * \author Chris Evers
 * \date Lund, 2011
 * \include wpDH.cpp
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
#ifdef NOSLIT
  cuboid con(in);
  #ifdef NOTAB                              // wpDH_noslit_notab
    typedef pot_r12debyehuckel Tpot;
  #else                                     // wpDH_noslit
    typedef pot_r12debyehuckel_hydrophobic_tab Tpot; 
  #endif
#else
  cuboidslit con(in);
  #ifdef NOTAB
    typedef pot_r12debyehuckelXY Tpot;      // wpDH_notab
  #else
    typedef pot_r12debyehuckelXY_hydrophobic_tab Tpot;
  #endif
  #ifdef HYDROPHIC                          // wpDH_hydrophobic
    typedef expot_hydrophobic Texpot;
  #else                                     // wpDH
//     typedef expot_gouychapman Texpot;
    typedef expot_gchydrophobic Texpot;
  #endif
#endif

  // Polymer
  polymer pol;
  pol.babeladd( con, in );          //  add from input
  con.trial=con.p;                  //  synchronize particle vector
  pol.masscenter(con);              //  update masscenter
  pol.move(con, -pol.cm+(con.slice_min+con.slice_max)*0.5);  // translate to the middle of the slice or the origo (0,0,0) if no slice defined
  pol.accept(con);                  //  accept translation
  cout << pol.info();

  // Potentials
#ifdef NOSLIT
  springinteraction<Tpot> pot(in); 
#else
  #ifdef PENALTY
    springinteraction_expot_penalty<Tpot, Texpot> pot(in, con, pol); 
    pot.pen.load("penalty.xy");  // Load penaltyfunction from disk 
    pot.pen.gofrload("gofr.xy"); // Load penaltyfunction from disk 
    pot.expot.update(con);
  #else
    springinteraction_expot<Tpot, Texpot> pot(in); 
    pot.expot.update(con);
  #endif
#endif

  // Distribution functions
  double* zhalfPtr=&con.len_half.z;         // half length of container in z-direction
  distributions dst;                        // distance dependent averages
  histogram gofr(0.1,0, (*zhalfPtr)*2. );   // radial distribution function
  histogram q(1,-pol.nb.size(), pol.nb.size()); // polymer charge distribution function
  histogram rg(.2,0, *zhalfPtr );      // radius of gyration distribution function
  histogram ree(.2,0, *zhalfPtr );     // end to end distance distribution function

  // Moves
  monomermove mm(nvt,con,pot,in);
  crankShaft cs(nvt,con,pot,in);
  branchRotation br(nvt,con,pot,in);
  macrorot mr(nvt,con,pot,in);
  translate mt(nvt,con,pot,in);
  mt.dpv.x=mt.dpv.y=0;                          // no need to translate in xy direction
  if (in.getflt("eqtit_runfrac",0.5)<1e-3) {
    pol.loadCharges(in.getstr("pol_charges","q.in"), con.p); // load residue charges from file
    con.trial=con.p;
    cout << "! Loaded polymer charges from disk" << endl;
  }
  eqtitrate tit(nvt, con, pot, in, "eqtit_");

  // Particle input output
  io io;
  iopqr pqr;
  ioxtc xtc(con.len.z);
  xtc.setbox(con.len.x,con.len.y,con.len.z);
  ioaam aam;
  aam.load(con, "confout.aam");     //   load stored configuration
  pol.masscenter(con);              //   update masscenter
  aam.save("confin.aam", con.p);

  // Initial system energy
  systemenergy sys(
    pot.energy( con.p, pol) +
    tit.intrinsicenergy(con.p) +
    pot.uself_polymer(con.p, pol) );

  cout << con.info() << atom.info()
    << pot.info() << tit.info() 
    << in.info();

  // Simulation loop
  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      switch (rand() % 6) {     // random number between 0 and 5
        case 0:
          for (int i=0; i<pol.size(); i++)
            sys+=mm.move(pol);   // move monomers 
          pol.masscenter(con);
          gofr.add(*zhalfPtr-pol.cm.z);
#ifdef PENALTY
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#endif
          break;
        case 1:
          sys+=mt.move(pol);     // translate polymer
          gofr.add(*zhalfPtr-pol.cm.z);
#ifdef PENALTY
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#endif
          break;
        case 2:
          pol.masscenter(con);
          sys+=mr.move(pol);     // rotate polymer
          break;
        case 3:
#ifdef PENALTY
          sys+=cs.penaltymove(pol);     // crankshaft
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#else
          sys+=cs.move(pol);     // crankshaft
#endif
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 4:
#ifdef PENALTY
          sys+=br.penaltymove(pol);     // branch rotation
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#else
          sys+=br.move(pol);     // branch rotation
#endif
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 5:
          sys+=tit.move();       // titrate titratable sites
          pol.charge(con.p);
          break;
      }

      if (slp.random_one()>0.3) {
        pol.masscenter(con);
        double z=*zhalfPtr-pol.cm.z;  // distance between masscenter and wall
        double charge=pol.charge(con.p);
        dst.add("qtot", z, charge);
        q.add(charge);

        double rg2=pol.sqmassgradius(con);
        dst.add("Rg2", z, rg2);
        dst.add("Rg", z, sqrt(rg2));
        rg.add(sqrt(rg2));

        double ree2=pol.sqend2enddistance(con);
        dst.add("Ree2", z, ree2); 
        dst.add("Ree", z, sqrt(ree2)); 
        ree.add(sqrt(ree2));

        for (int i=pol.beg; i<=pol.end; i++) {
          std::ostringstream s; s << "q" << i;
          dst.add( s.str(), *zhalfPtr-con.p[i].z, con.p[i].charge );	
        }
      }

      if (slp.random_one()>0.95) {
        xtc.save( "traj.xtc", con.p );
      }
    } // END of micro loop

    sys.update(
      pot.energy( con.p, pol) +
      tit.intrinsicenergy(con.p) +
      pot.uself_polymer(con.p, pol) );

    cout << loop.timing() << "#   Energy drift = " << sys.cur-sys.sum << " kT. "
      << "System charge = " << con.charge() << ". " << endl;

    aam.save("confout.aam", con.p);
    // Write files to disk
    if ( in.getboo("write_files",true) ) {
      io.writefile("vmdbonds.tcl", pol.getVMDBondScript());
      pqr.save("confout.pqr",con.p);
      dst.write("dist.out");
      dst.cntwrite("cntdist.out");
      gofr.write("gofr.out");
      gofr.dump("gofr.xy");

      q.dump("q.xy");
      rg.dump("rg.xy");
      ree.dump("ree.xy");

// #ifndef NOTAB
//       pot.pair.dump("pairpot.xy");
// #endif
#ifndef NOSLIT
      pot.expot.dump("expot.xy");
#endif
#ifdef PENALTY
      pot.pen.dump("penalty", loop.cnt_macro, "xy");
      pot.pen.gofrdump("gofr", loop.cnt_macro, "xy");
#endif

      tit.applycharges(con.trial);                      // Set average charges on all titratable sites
      pol.saveCharges("q.out", con.trial);              // Save average charges to disk
      con.trial=con.p;                                  // Restore original charges
    }

  } // END of macro loop and simulation

cout << sys.info() << loop.info() << mm.info()
  << cs.info() << br.info() << mr.info() 
  << mt.info() << tit.info(con.p) << pol.info() 
  << pot.info();

}
