/*! \page test_wallprotein WallproteinDH
 *
 * Simulate a number of flexible and titratable
 * polymers in a salt solution in the presence
 * of a charged Gouy-Chapman surface with implicit
 * salt interactions given by Debye Huckel.
 *
 * \author Chris Evers
 * \date Lund, 2011
 * \include wp.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/energy/externalpotential.h"
#include "faunus/potentials/pot_r12debyehuckel.h"
//#include "faunus/energy/penaltyfunction.h"

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
#ifdef NOSLIT                               // If wpDH_noslit
  cuboid con(in);
  springinteraction<pot_r12debyehuckel> pot(in);
#else                                       // If wpDH
  cuboidslit con(in);
  springinteraction_expot<pot_r12debyehuckelXY, expot_gouychapman> pot(in);
  pot.expot.update(con);                    // initialize Gouy-Chapman potential
#endif
  double* zhalfPtr=&con.len_half.z;         // half length of container in z-direction
  distributions dst;                        // distance dependent averages
  histogram gofr(0.1,0, (*zhalfPtr)*2. );   // radial distribution function

  // Polymer
  polymer pol;                      // Polymer
  pol.babeladd( con, in );          //  add from input
  con.trial=con.p;                  //  synchronize particle vector
  pol.masscenter(con);              //  update masscenter
  pol.move(con, -pol.cm+(con.slice_min+con.slice_max)*0.5);  // translate to the middle of the slice or the origo (0,0,0) if no slice defined
  pol.accept(con);                  //  accept translation
  cout << pol.info();

  // Moves
  monomermove mm(nvt,con,pot,in);
  crankShaft cs(nvt,con,pot,in);
  macrorot mr(nvt, con, pot,in);
  translate mt(nvt, con, pot, in);
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
  ioaam aam;
  aam.load(con, "confout.aam");     //   load stored configuration
  pol.masscenter(con);              //   update masscenter
  aam.save("confout_init.aam", con.p);

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
      switch (rand() % 5) {     // random number between 0 and 4
        case 0:
          for (int i=0; i<pol.size(); i++)
            sys+=mm.move(pol);   // move monomers 
          pol.masscenter(con);
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 1:
          sys+=mt.move(pol);     // translate polymers
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 2:
          pol.masscenter(con);
          sys+=mr.move(pol);     // rotate polymers
          break;
        case 3:
          sys+=cs.move(pol);     // crankshaft
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 4:
          sys+=tit.move();       // titrate titratable sites
          pol.charge(con.p);
          break;
      }

      if (slp.random_one()>0.3) {
        pol.masscenter(con);
        double z=*zhalfPtr-pol.cm.z;  // distance between masscenter and wall
        dst.add("qtot", z, pol.charge(con.p));
        dst.add("Rg", z, sqrt(pol.sqmassgradius(con)));
        dst.add("Rg2", z, pol.sqmassgradius(con));
        dst.add("Ree", z, sqrt(pol.sqend2enddistance(con))); 
        dst.add("Ree2", z, pol.sqend2enddistance(con)); 

        for (int i=pol.beg; i<=pol.end; i++) {
          std::ostringstream s; s << "q" << i;
          dst.add( s.str(), *zhalfPtr-con.p[i].z, con.p[i].charge );	
        }
      }

      if (slp.random_one()>0.95) {
        xtc.setbox(con.len.x,con.len.y,con.len.z);
        xtc.save( "traj.xtc", con.p );
      }
    } // END of micro loop

    sys.update(
      pot.energy( con.p, pol) +
      tit.intrinsicenergy(con.p) +
      pot.uself_polymer(con.p, pol) );

    cout << loop.timing() << "#   Energy drift = " << sys.cur-sys.sum << " kT. "
      << "System charge = " << con.charge() << ". "
      << endl;

    // Write files to disk
    io.writefile("vmdbonds.tcl", pol.getVMDBondScript());
    pqr.save("confout.pqr",con.p);
    aam.save("confout.aam", con.p);
    dst.write("dist.out");
    dst.cntwrite("cntdist.out");
    gofr.write("gofr.out");
#ifndef NOSLIT
    pot.expot.dump("expot.xy");
#endif

    tit.applycharges(con.trial);                      // Set average charges on all titratable sites
    pol.saveCharges("q.out", con.trial);              // Save average charges to disk
    con.trial=con.p;                                  // Restore original charges

  } // END of macro loop and simulation

cout << sys.info() << loop.info() << mm.info()
  << cs.info() << mr.info() << mt.info()
  << tit.info(con.p) << pol.info() << pot.info();

}
