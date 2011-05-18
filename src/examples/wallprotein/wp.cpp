/*! \page test_wallprotein Wallprotein
 *
 * Simulate a number of flexible and titratable
 * polymers in a salt solution in the presense
 * of a charged surface.
 *
 * \author Mikael Lund / Chris Evers
 * \date Lund, 2009-2011
 * \include wp.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/energy/externalpotential.h"
#include "faunus/energy/penaltyfunction.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/potentials/pot_coulomb.h"

using namespace std;
using namespace Faunus;

int main() {
  // General
  cout << faunus_splash();
  inputfile in("wp.conf");
  slump slp;                                // load random number generator
  mcloop loop(in);                          // define simulation loops
  grandcanonical nmt;                       // grand canonical ensemble

  // Simulation container
#ifdef NOSLIT
  cuboid con(in);
  typedef pot_r12minimage Tpot;
#else
  cuboidslit con(in);
  typedef pot_r12minimageXY Tpot;
#endif

  // Polymers
  polymer pol;
//   ioaam aam;
//   if (in.getstr("polymer")=="aam")          // Protein input file format is AAM
//     if (in.getboo("lattice")==true)         // Protein input file format is AAM
//       aam.loadlattice(con, in, pol);         // Load and insert proteins
//     else                                    // Center first protein (will be frozen)
//       aam.load(con, in, pol);
//   else
  pol.babeladd( con, in );                  //  add from input
  con.trial=con.p;                          //  synchronize particle vector
  pol.masscenter(con);                      //  update masscenter
  pol.move(con, -pol.cm+(con.slice_min+con.slice_max)*0.5);  // translate to the middle of the slice or the origo (0,0,0) if no slice defined
  pol.accept(con);                          //  accept translation
  cout << pol.info();

  // Distribution functions
  double* zhalfPtr=&con.len_half.z;         // half length of container in z-direction
  distributions2 dst(0.1,0,(*zhalfPtr)*2.); // distance dependent averages
  histogram gofr(0.1,0, (*zhalfPtr)*2. );   // radial distribution function
  histogram q(1,-pol.nb.size(), pol.nb.size()); // polymer charge distribution function
  histogram rg(.2,0, (*zhalfPtr)*2. );      // radius of gyration distribution function
  histogram ree(.2,0, (*zhalfPtr)*2. );     // end to end distance distribution function

  // Potentials
#ifdef NOSLIT
  springinteraction<Tpot> pot(in); 
#else
  #ifdef PENALTY
    springinteraction_expot_penalty<Tpot, expot_akesson> pot(in, con, pol);
    pot.pen.load("penalty.xy");               // Load penaltyfunction from disk 
    pot.pen.gofrload("gofr.xy");              // Load penaltyfunction from disk 
    pot.expot.update(con);
  #else
    springinteraction_expot<Tpot, expot_akesson> pot(in); 
    pot.expot.update(con);
  #endif
#endif

  // Handle wall particles
  group wall;
  wall.add(con, atom[in.getstr("wall_tion1")].id, in.getint("wall_nion1") );
  for (int i=wall.beg; i<=wall.end; i++)
    con.p[i].z=*zhalfPtr;                   // move all particles to edge
  con.trial=con.p;

  // Handle salt particles
  salt salt( atom["NA"].id, atom["CL"].id ); 
  salt.add(con,in);

  // Moves
  monomermove mm(nmt,con,pot,in);           // Rattle MC move
  crankShaft cs(nmt,con,pot,in);            // ...crankshaft
  branchRotation br(nmt,con,pot,in);        // ...branch rotation
  macrorot mr(nmt, con, pot,in);            // ...rotate
  translate mt(nmt, con, pot, in);          // ...translate
  mt.dpv.x=mt.dpv.y=0;                      // ...no need to translate in xy direction

  saltmove wm(nmt,con,pot,in,"wall_");      // Wall particle moves
  wm.name="WALL PARTICLE DISPLACEMENTS";
  wm.dpv.z=0;                               // ...constrain displacements to the xy-plane

  saltmove sm(nmt,con,pot,in);              // Salt moves
  saltbath sb(nmt,con,pot,in,salt);         // Grand canonical salt bath

  if (in.getflt("tit_runfrac",0.5)<1e-3)    // No titration?
    for (int i=0; i<atom.list.size(); i++) 
      atom[i].pka = 0.;                     // .. set pKa values to zero
  GCchargereg tit(nmt,con,pot,in);          // Grand canonical proton titration
  if (in.getflt("tit_runfrac",0.5)<1e-3)    // No titration?
    if (pol.loadCharges(in.getstr("pol_charges","q.in"), con.p)) // .. load polymer charges from file
      con.trial=con.p;
    else {
      cerr << "!! charge file not loaded !!" << endl;
      return 0;
    }

  // File I/O
  io io;
  iopqr pqr;
  ioxtc xtc(con.len.z);
  ioqtraj qtraj;
  xtc.setbox(con.len.x,con.len.y,con.len.z);
  vector<group> vg;                         // Particles in xtc-trajectories
  vg.push_back(pol);
  vg.push_back(wall);
  ioaam aam;
  if ( nmt.load(con, "gcgroup.conf")==true ) {
    aam.load(con, "confout.aam");           // Load stored configuration
    pol.masscenter(con);                    // Update masscenter
    aam.save("confout_init.aam", con.p);
  }

  // Neutralize system charge
  if (con.p[salt.beg].id==atom["ghost"].id)
    con.p[salt.beg].charge=0;
  double qtot=con.charge();                         // Total system charge
  int qint=int(floor(qtot));
  if (qint < 0) {                                   // Negative charge
    salt.group::add(con, atom["NA"].id, -qint );    // ..add cations
    cout << "# Added " << -qint << " Na-";
  } else  {                                         // Positive charge
    salt.group::add(con, atom["CL"].id, qint );     // ..add anions
    cout << "# Added " << qint << " Cl-";
  }
  if (con.p[salt.beg].id==atom["ghost"].id) {
    double qghost=qtot-qint;                         // Noninteger charge (always positive)
    con.p[salt.beg].charge=-qghost;                  //   put on ghost particle
    con.trial[salt.beg].charge=con.p[salt.beg].charge;
    cout << " and set the ghost particle charge to " << -qghost << " to obtain electroneutrality" << endl;
  }

  // Initial system energy
  double utot=pot.energy( con.p, wall, salt) +
      pot.energy( con.p, wall, pol) +
      pot.energy( con.p, salt, pol) +
      pot.uself_polymer(con.p, pol) + 
      pot.internal( con.p, wall ) +
      pot.internal( con.p, salt );
#ifndef NOSLIT
  utot+=pot.expot.energy_group(con.p, pol);
#endif
#ifdef PENALTY
  utot+=pot.pen.energy(*zhalfPtr-pol.cm.z);
#endif
  systemenergy sys(utot);

  cout << con.info() << atom.info()
      << pot.info() << salt.info(con)
      << in.info() << tit.info();

  // Simulation loop
  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      switch (rand() % 9) {
        case 0:
          sys+=sm.move(salt);           // salt moves
          break;
        case 1:
          sys+=wm.move(wall);           // move surface charges
          break;
        case 2:
          for (int i=0; i<pol.size(); i++)
            sys+=mm.move(pol);          // move monomers
          pol.masscenter(con);
          gofr.add(*zhalfPtr-pol.cm.z);
#ifdef PENALTY
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#endif
          break;
        case 3:
          sys+=mt.move(pol);            // translate polymers
          gofr.add(*zhalfPtr-pol.cm.z);
#ifdef PENALTY
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#endif
          break;
        case 4:
          pol.masscenter(con);
          sys+=mr.move(pol);            // rotate polymers
          break;
        case 5:
#ifdef PENALTY
          sys+=cs.penaltymove(pol);     // crankshaft
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#else
          sys+=cs.move(pol);
#endif
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 6:
#ifdef PENALTY
          sys+=br.penaltymove(pol);     // branchrot
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
#else
          sys+=br.move(pol);
#endif
          gofr.add(*zhalfPtr-pol.cm.z);
          break;
        case 7:
          sys+=sb.move();               // grand canonical salt move
          break;
        case 8:
          sys+=tit.titrateall();        // titrate titratable sites
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

        for (int i=salt.beg; i<=salt.end; i++) 
          if (con.p[i].id==atom["NA"].id) 
            dst.add("Na", *zhalfPtr-con.p[i].z, *zhalfPtr-con.p[i].z);
          else if (con.p[i].id==atom["CL"].id) 
            dst.add("Cl", *zhalfPtr-con.p[i].z, *zhalfPtr-con.p[i].z);
      }

#ifndef NOSLIT
      if (slp.random_one()>0.8)
        sys += pot.expot.update(con);
#endif

      if (slp.random_one()<in.getflt("traj_runfrac",0.05) ) {
        xtc.save( "traj.xtc", con.p, vg );
        qtraj.save( "q.traj", con.p, vg );
      }

    } // END of micro loop

    utot=pot.energy( con.p, wall, salt) +
        pot.energy( con.p, wall, pol) +
        pot.energy( con.p, salt, pol) +
        pot.uself_polymer(con.p, pol) + 
        pot.internal( con.p, wall ) +
        pot.internal( con.p, salt );
#ifndef NOSLIT
    utot+=pot.expot.energy_group(con.p, pol);
#endif
#ifdef PENALTY
    utot+=pot.pen.energy(*zhalfPtr-pol.cm.z);
#endif
    sys.update(utot);

    cout << loop.timing() << "#   Energy drift = " << sys.cur-sys.sum << " kT. "
        << "System charge = " << con.charge() << ". " << endl;

    // Write files to disk
    if ( in.getboo("write_files",true) ) {
      io.writefile("vmdbonds.tcl", pol.getVMDBondScript());
      pqr.save("confout.pqr",con.p);
      pqr.save("confout_polwall.pqr",con.p,vg);
      aam.save("confout.aam", con.p);
      io.writefile("gcgroup.conf", nmt.print()); 
      dst.write("dist.out");
      dst.cntwrite("cntdist.out");
      gofr.write("gofr.out");
      gofr.dump("gofr.xy");

      q.dump("fluctQ.xy");
      rg.dump("fluctRg.xy");
      ree.dump("fluctRee.xy");

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

  cout << sys.info() << sm.info() << wm.info() 
      << "#   Surf area / charge (AA^2) = " << (con.len.x*con.len.y)/wall.size() << endl 
      << mm.info() << cs.info() << br.info() << mr.info() 
      << mt.info() << tit.info() << sb.info() << pol.info()
      << pot.info();
}
