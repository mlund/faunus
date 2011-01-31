/*! \page test_wallprotein Wallprotein
 *
 * Simulate a number of flexible and titratable
 * polymers in a salt solution in the presense
 * of a charged surface.
 *
 * \author Mikael Lund
 * \date Lund, 2009-2010
 * \include wp.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/energy/externalpotential.h"
#include "faunus/energy/penaltyfunction.h"
#include "faunus/potentials/pot_hsminimage.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/potentials/pot_coulomb.h"

using namespace std;
using namespace Faunus;

#ifdef NOSLIT
typedef pot_r12minimage Tpot;
typedef cuboid Tcon;
#else
typedef pot_r12minimageXY Tpot;
typedef cuboidslit Tcon;
#endif

// Redefine energy function
template<class Tpairpot, class Texpot>
class springinteraction_expot_penalty : public springinteraction_expot<Tpot,Texpot> {
  private:
    container* conPtr;
    group* groupPtr;
    double zlen_half;
  public:
    using springinteraction_expot<Tpairpot,Texpot>::energy;
    penaltyfunction pen;

    springinteraction_expot_penalty(inputfile &in, container &c, polymer &p)
      : springinteraction_expot<Tpot,Texpot>(in),
      pen( 0, in.getflt("cuboid_zlen",100), 0.25, in.getflt("penalty_energy",.2), in.getflt("penalty_scalingfactor",1)) {
          zlen_half=.5*in.getflt("cuboid_zlen",100);
          conPtr=&c;
          groupPtr=&p;
          interaction<Tpairpot>::name+=" + penalty function";
    }

    string info() {
      std::ostringstream o;
      o << springinteraction_expot<Tpot,Texpot>::info() << pen.info();
      return o.str();
    }

    // Energy in momomer and crankshaft functions
    double u_monomer(vector<particle> &p, const polymer &g, unsigned int i ) {
      double u=0;
      if (&g==groupPtr) {
        point cm = g.masscenter(*conPtr, p);
        u=pen.energy(zlen_half-cm.z);
      }
      return springinteraction_expot<Tpairpot,Texpot>::u_monomer(p,g,i) + u;
    }

    // Energy in translation function
    double energy(vector<particle> &p, const group &g) {  
      double u=0;
      if (&g==groupPtr) {
        point cm = g.masscenter(*conPtr, p);
        u=pen.energy(zlen_half-cm.z);
      }
      return springinteraction_expot<Tpairpot,Texpot>::energy(p,g) + u;
    }
};

int main() {

  // General setup
  cout << faunus_splash();
  slump slp;
  inputfile in("wp.conf");
  mcloop loop(in);
  grandcanonical nmt;
  Tcon con(in);
  double* zhalfPtr;
  zhalfPtr=&con.len_half.z;
  distributions dst;    // Distance dep. averages
  histogram gofr(0.1,0, (*zhalfPtr)*2. );

  // Handle polymers
  polymer pol;
#ifdef BABEL
  pol.babeladd( con, in );
#endif
  //  atom.reset_properties(con.p);  //rudimentary with new babeladd
  con.trial=con.p;
  pol.masscenter(con);
  springinteraction_expot_penalty<Tpot, expot_akesson> pot(in, con, pol);
  pot.pen.load("penalty.xy");  // Load penaltyfunction from disk 
  pot.pen.gofrload("gofr.xy"); // Load penaltyfunction from disk 
  //   springinteraction_expot<Tpot, expot_akesson> pot(in);

  pol.move(con, -pol.cm+(con.slice_min+con.slice_max)*0.5);         // Translate polymer to the middle of the slice or the origo (0,0,0) if no slice defined
  pol.accept(con);                // .. accept translation
  monomermove mm(nmt,con,pot,in); // Rattle MC move
  crankShaft cs(nmt,con,pot,in);  // ...crankshaft
  macrorot mr(nmt, con, pot,in);  // ...rotate
  translate mt(nmt, con, pot, in);// ...translate
  mt.dpv.x=mt.dpv.y=0;            // ...no need to translate in xy direction
  cout << pol.info();

  // Handle wall particles
  group wall;
  wall.add(con, atom[in.getstr("wall_tion1")].id, in.getint("wall_nion1") );
  for (int i=wall.beg; i<=wall.end; i++)
    con.p[i].z=*zhalfPtr; // move all particles to edge
  con.trial=con.p;
  saltmove wm(nmt,con,pot,in,"wall_");
  wm.name="WALL PARTICLE DISPLACEMENTS";
  wm.dpv.z=0; // constrain displacements to the xy-plane

  vector<group> vg;
  vg.push_back(pol);
  vg.push_back(wall);

  // Handle salt particles
  salt salt( atom["NA"].id, atom["CL"].id ); 
  salt.add(con,in);
  saltmove sm(nmt,con,pot,in);

  // File I/O
  io io;
  iopqr pqr;
  ioaam aam;
  ioxtc xtc(con.len.z);

  // No titration?
  if (in.getflt("tit_runfrac",0.5)<1e-3) {
    for (int i=0; i<atom.list.size(); i++) 
      atom[i].pka = 0.;
  }

  // Grand canonical stuff
  saltbath sb(nmt,con,pot,in,salt);
  GCchargereg tit(nmt,con,pot,in);

  // Load stored configuration?
  if ( nmt.load(con, "gcgroup.conf")==true ) {
    aam.load(con, "confout.aam");
    pol.masscenter(con);
  }

  // No titration?
  if (in.getflt("tit_runfrac",0.5)<1e-3) {
    if (pol.loadCharges("q.in", con.p))                   // Load polymer charges from file
      con.trial=con.p;
    else {
      cerr << "!! charge file not loaded !!" << endl;
      return 0;
    }
  }

  // Neutralize residual charge with counterions and smearing
  if (con.p[salt.beg].id==atom["ghost"].id)          // If ghost particle
    con.p[salt.beg].charge=0;                       //   set its charge to zero
  double q;
  int qint;
  q = con.charge();                                 // Total system charge
  qint = int(floor(q));                                  // Integer system charge
  if (qint < 0) {                                   // Negative charge
    salt.group::add(con, atom["NA"].id, -qint );    //   add cations
    cout << "# Added " << -qint << " Na-";
  } else  {                                          // Positive charge
    salt.group::add(con, atom["CL"].id, qint );     //   add anions
    cout << "# Added " << qint << " Cl-";
  }
  if (con.p[salt.beg].id==atom["ghost"].id) {        // If ghost particle
    q = q-qint;                                     // Noninteger charge (always positive)
    con.p[salt.beg].charge=-q;                      //   put on ghost particle
    con.trial[salt.beg].charge=con.p[salt.beg].charge;
    cout << " and set the ghost particle charge to " << -q << " to obtain electroneutrality" << endl;
  }

  // Calculate initial system energy
  systemenergy sys(
      pot.energy( con.p, wall, salt) +
      pot.energy( con.p, wall, pol) +
      pot.energy( con.p, salt, pol) +
      pot.uself_polymer(con.p, pol) + pot.expot.energy(con.p, pol) +
      pot.internal( con.p, wall ) +
      pot.internal( con.p, salt ) +
      pot.pen.energy(*zhalfPtr-pol.cm.z) );
  
  cout << con.info() << atom.info()
    << pot.info() << salt.info(con)
    << in.info() << tit.info();

  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      switch (rand() % 8) {
        case 0:
          sys+=sm.move(salt);    // salt moves
          break;
        case 1:
          sys+=wm.move(wall);    // move surface charges
          break;
        case 2:
          for (int i=0; i<pol.size(); i++)
            sys+=mm.move(pol);   // move monomers
          pol.masscenter(con);
          gofr.add(*zhalfPtr-pol.cm.z);
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
          break;
        case 3:
          sys+=mt.move(pol);     // translate polymers
          gofr.add(*zhalfPtr-pol.cm.z);
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
          break;
        case 4:
          pol.masscenter(con);
          sys+=mr.move(pol);     // rotate polymers
          break;
        case 5:
          sys+=cs.move(pol);     // crankshaft
          gofr.add(*zhalfPtr-pol.cm.z);
          sys+=pot.pen.update(*zhalfPtr-pol.cm.z);
          break;
        case 6:
          sys+=sb.move();        // grand canonical salt move
          break;
        case 7:
          sys+=tit.titrateall(); // titrate titratable sites
          pol.charge(con.p);
          break;
      }

      if (slp.random_one()>0.3) {
        pol.masscenter(con);
        double z=*zhalfPtr-pol.cm.z;
        dst.add("qtot", z, pol.charge(con.p));
        dst.add("Rg2", z, pol.sqmassgradius(con));
        dst.add("Ree2", z, pol.sqend2enddistance(con)); 

        for (int i=salt.beg; i<=salt.end; i++) {
          if (con.p[i].id==atom["NA"].id) {
            dst.add("Na", *zhalfPtr-con.p[i].z, *zhalfPtr-con.p[i].z);
          }
          if (con.p[i].id==atom["CL"].id) {
            dst.add("Cl", *zhalfPtr-con.p[i].z, *zhalfPtr-con.p[i].z);
          }
        }
        for (int i=pol.beg; i<=pol.end; i++) {
          std::ostringstream s; s << "q" << i;
          dst.add( s.str(), *zhalfPtr-con.p[i].z, con.p[i].charge );	
        }
      }

      if (slp.random_one()>0.8)
        sys += pot.expot.update(con);

      if (slp.random_one()>0.95) {
        xtc.setbox(con.len.x,con.len.y,con.len.z);
        xtc.save( "traj.xtc", con.p, vg );
      }
    } // END of micro loop

    pot.expot.save(con);

    sys.update(
        pot.energy( con.p, wall, salt) +
        pot.energy( con.p, wall, pol) +
        pot.energy( con.p, salt, pol) +
        pot.uself_polymer(con.p, pol) + pot.expot.energy(con.p, pol) +
        pot.internal( con.p, wall ) +
        pot.internal( con.p, salt ) +
        pot.pen.energy(*zhalfPtr-pol.cm.z) );

    cout << loop.timing() << "#   Energy drift = " << sys.cur-sys.sum << " kT. "
      << "System charge = " << con.charge() << ". "
      << "New penalty energy = " << pot.pen.scaledu() << endl;

    // Write files to disk
    io.writefile("vmdbonds.tcl", pol.getVMDBondScript());
    pqr.save("confout.pqr",con.p);
    pqr.save("confout-polwalll.pqr",con.p,vg);
    aam.save("confout.aam", con.p);
    dst.write("dist.dat");
    dst.cntwrite("cntdist.dat");
    gofr.write("gofr.dat");
    io.writefile("gcgroup.conf", nmt.print());
    pot.pen.write("penalty.dat");
    pot.pen.dump("penalty", loop.cnt_macro, "xy");
    pot.pen.gofrdump("gofr", loop.cnt_macro, "xy");    

    tit.applycharges(con.trial);                      // Set average charges on all titratable sites
    pol.saveCharges("q.out", con.trial);              // Save average charges to disk
    con.trial=con.p;                                  // Restore original charges

  } // END of macro loop and simulation

  cout << sys.info() << sm.info() << wm.info() << loop.info()
    << mm.info()  << cs.info() << mr.info() << mt.info()
    << tit.info() << sb.info() << pol.info()<< pot.info();
}
