/*! \page test_wallprotein Wallprotein
 *
 * Simulate a number of flexible and titratable
 * polymers in a salt solution in the presense
 * of a charged surface.
 *
 * \author Mikael Lund
 * \date Lund, 2009
 * \include wp.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/potentials/pot_hsminimage.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  // General setup
  slump slump;
  cout << faunus_splash();
  inputfile in("wp.conf");
  mcloop loop(in);
  canonical nvt;
  slit con(in);
  springinteraction<pot_hsminimageXY> pot(in);
  distributions dst;    // Distance dep. averages

  // Handle polymers
  polymer pol;
#ifdef BABEL
  pol.babeladd( con, in );
#endif
  pol.move(con, -pol.cm);         // Translate polymer to origo (0,0,0)
  pol.accept(con);                // .. accept translation
  monomermove mm(nvt,con,pot,in); // Rattle MC move
  crankShaft cs(nvt,con,pot,in);  // ...crankshaft
  macrorot mr(nvt, con, pot);     // ...rotate
  translate mt(nvt, con, pot, in);// ...translate
  cout << pol.info();

  // Handle wall particles
  group wall;
  wall.add(con, atom[in.getstr("wall_tion1")].id, in.getint("wall_nion1") );
  for (int i=wall.beg; i<=wall.end; i++)
    con.p[i].z=con.len_half; // move all particles to edge
  con.trial=con.p;
  saltmove wm(nvt,con,pot,in);
  wm.name="WALL PARTICLE DISPLACEMENTS:";
  wm.dp=in.getflt("wall_dp",10);
  wm.dpv.z=0; // constrain displacements to the xy-plane

  // Handle salt particles
  salt salt( atom["NA"].id, atom["CL"].id ); 
  salt.add(con,in);
  saltmove sm(nvt,con,pot,in);
  group mobile = wall + salt;

  // File I/O
  iopqr pqr;
  ioaam aam;
  ioxtc xtc(con.len);
  aam.load(con,"conf.aam");

  chargereg tit(nvt,con,pot,salt,in.getflt("pH",7.));

  systemenergy sys(
        pot.internal(con.p, mobile) +
        pot.energy(con.p, mobile, pol) +
        pot.uself_polymer(con.p,pol) );

  cout << con.info() << atom.info()
       << pot.info() << salt.info(con)
       << in.info() << tit.info();

  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      switch (rand() % 7) {
        case 0:
          sys+=sm.move(salt);    // salt moves
          break;
        case 1:
          sys+=wm.move(wall);    // move surface charges
          break;
        case 2:
          sys+=mm.move(pol);     // move monomers
          break;
        case 3:
          sys+=mt.move(pol);     // translate polymers
          break;
        case 4:
          sys+=mr.move(pol);     // rotate polymers
          break;
        case 5:
          sys+=cs.move(pol);     // crankshaft
          break;
        case 6:
          sys+=tit.titrateall(); // titrate titratable sites
          dst.add("poltot", abs(pol.cm.z), pol.cm.charge);
          for (int i=pol.beg; i<=pol.end; i++) {
            std::ostringstream s; s << "pol" << i;
            dst.add( s.str(), abs(con.p[i].z), con.p[i].charge );
          }
      }

      if (slump.random_one()>0.9)
        xtc.save( "traj.xtc", con.p );
    }                                   // END of micro loop

    sys.update(
        pot.internal(con.p, mobile) +
        pot.energy(con.p, mobile, pol) +
        pot.uself_polymer(con.p,pol) );
    aam.save("conf.aam",con.p);
    pqr.save("conf.pqr",con.p);
    dst.write("dist.dat");
    cout << loop.timing();     
  }                                     // END of macro loop and simulation

  cout << sys.info() << sm.info() << wm.info() << loop.info()
       << mm.info() << cs.info() << mr.info() << mt.info() << tit.info();
}
