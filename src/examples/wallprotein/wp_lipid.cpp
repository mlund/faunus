/*! 
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
#include "faunus/moves/membrane.h"
#include "faunus/potentials/pot_hsminimage.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/potentials/pot_coulomb.h"

using namespace std;
using namespace Faunus;

#ifdef NOSLIT
	typedef pot_r12minimage Tpot;
	typedef box Tcon;
#else
	typedef pot_r12minimageXY Tpot;
	typedef slit Tcon;
#endif

int main() {
  // General setup
  cout << faunus_splash();
  inputfile in("wp.conf");
  mcloop loop(in);
  grandcanonical nmt;
  Tcon con(in);
  springinteraction<Tpot> pot(in);
  distributions dst;    // Distance dep. averages
  histogram gofr(0.1,0,con.len);

  // Handle polymers
  polymer pol;
#ifdef BABEL
  pol.babeladd( con, in );
#endif
  atom.reset_properties(con.p);
  con.trial=con.p;
  pol.move(con, -pol.cm);         // Translate polymer to origo (0,0,0)
  pol.accept(con);                // .. accept translation
  monomermove mm(nmt,con,pot,in); // Rattle MC move
  crankShaft cs(nmt,con,pot,in);  // ...crankshaft
  macrorot mr(nmt, con, pot);     // ...rotate
  translate mt(nmt, con, pot, in);// ...translate
  mt.dpv.x=mt.dpv.y=0;            // ...no need to translate in xy direction

  // Handle lipids
  popscmembrane mem;
  mem.load(in,con);
  membranemove memm(nmt, con, pot, mem);
  memm.dpmm=in.getflt("monomer_dp",3);
  memm.dpgp=in.getflt("graftpoint_dp",2);
  memm.dplatt=in.getflt("latticetranslation_dp",4);
  
  // Handle salt particles
  salt salt( atom["NA"].id, atom["CL"].id ); 
  salt.add(con,in);
  saltmove sm(nmt,con,pot,in);
 
  // File I/O
  io io;
  iopqr pqr;
  ioaam aam;
  ioxtc xtc(con.len);

  // Grand canonical stuff
  saltbath sb(nmt,con,pot,in,salt);
  GCchargereg tit(nmt,con,pot,in);

  // Load stored configuration?
  if ( nmt.load(con, "gcgroup.conf")==true )                                       
    aam.load(con, "confout.aam");

  // Calculate initial system energy
  systemenergy sys(
                   pot.energy( con.p, mem, salt) +
                   pot.energy( con.p, mem, pol) +
                   pot.energy( con.p, salt, pol) +
                   pot.uself_polymer(con.p, pol) +
                   pot.uself_popscmem(con.p, mem) +
                   pot.internal( con.p, salt ));
 
  cout << con.info() << atom.info() << in.info()
       << pot.info() << salt.info(con)
       << tit.info() << mem.info(con) << pol.info()
       << errlog.info() << endl;

  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      switch (rand() % 8) {
        case 0:
          sys+=sm.move(salt);    // salt moves
          break;
        case 1:
          sys+=memm.move(mem);   // move lipids
          break;
        case 2:
          for (int i=0; i<pol.size(); i++)
            sys+=mm.move(pol);   // move monomers
          pol.masscenter(con);
          gofr.add(pol.cm.z+con.len_half);
          break;
        case 3:
          sys+=mt.move(pol);     // translate polymer
          gofr.add(pol.cm.z+con.len_half);
          break;
        case 4:
          sys+=mr.move(pol);     // rotate polymer
          break;
        case 5:
          sys+=cs.move(pol);     // crankshaft
          gofr.add(pol.cm.z+con.len_half);
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
        dst.add("qtot", pol.cm.z+con.len_half, pol.charge(con.p));
        for (int i=pol.beg; i<=pol.end; i++) {
          std::ostringstream s; s << "q" << i;
          dst.add( s.str(), pol.cm.z+con.len_half, con.p[i].charge );
        }
      }

      if (slp.random_one()>0.95)
        xtc.save( "traj.xtc", con.p );
    } // END of micro loop

    sys.update(
        pot.energy( con.p, mem, salt) +
        pot.energy( con.p, mem, pol) +
        pot.energy( con.p, salt, pol) +
        pot.uself_polymer(con.p, pol) +
        pot.uself_popscmem(con.p, mem) +
        pot.internal( con.p, salt ));

    cout << loop.timing() << "# Energy drift = " << sys.cur-sys.sum << " kT" << endl;     

    // Write files to disk
    io.writefile("vmdbonds.tcl", pol.getVMDBondScript() + mem.getVMDBondScript());
    aam.save("confout.aam",con.p);
    pqr.save("confout.pqr",con.p);
    dst.write("dist.dat");
    gofr.write("gofr.dat");
    io.writefile("gcgroup.conf", nmt.print());

  } // END of macro loop and simulation

  cout << sys.info() << loop.info() << sm.info()
       << mm.info()  << cs.info() << mr.info() << mt.info()
       << tit.info() << sb.info() << memm.info() << pol.info();
}
