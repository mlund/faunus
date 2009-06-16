/*! \page test_wallprotein Wallprotein
 *
 * Simulate a number of flexible polymers in a salt
 * solution.
 *
 * \author Mikael Lund
 * \date Lund, 2009
 * \include wp.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  // General setup
  slump slump;
  cout << faunus_splash();
  inputfile in("wp.conf");
  //atom.load(in);
  mcloop loop(in);
  canonical nvt;
  slit con(in);
  springinteraction<pot_hsminimageXY> pot(in);

  // Handle polymers
  polymer pol;
  pol.babeladd( con, in );
  monomermove mm(nvt,con,pot,in); // ...rattle
  macrorot mr(nvt, con, pot);     // ...rotate
  translate mt(nvt, con, pot);    // ...translate
  cout << pol.info();

  // Handle wall particles
  group wall;
  wall.add(con, atom[in.getstr("wall_tion1")].id, in.getint("wall_nion1") );
  for (int i=wall.beg; i<=wall.end; i++)
    con.p[i].z=con.len_half;
  con.trial=con.p;
  saltmove wm(nvt,con,pot,in);
  wm.name="WALL PARTICLE DISPLACEMENTS:";
  wm.dp=in.getflt("wall_dp",10);
  wm.dpv.z=0;

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
      switch (rand() % 5) {
        case 0:
          sys+=sm.move(salt);
          break;
        case 1:
          sys+=wm.move(wall);
          break;
        case 2:
          sys+=mm.move(pol); 
          break;
        case 3:
          sys+=mt.move(pol);
          break;
        case 4:
          sys+=mr.move(pol);
          break;
      }
      if (slump.random_one()>0.5)
        xtc.save( "tis", con.p );

    }                                   // END of micro loop

    sys.update(
        pot.internal(con.p, mobile) +
        pot.energy(con.p, mobile, pol) +
        pot.uself_polymer(con.p,pol) );
    aam.save("conf.aam",con.p);
    pqr.save("conf.pqr",con.p);
    cout << loop.timing();     
  }                                     // END of macro loop and simulation

  cout << sys.info() << sm.info() << wm.info() << loop.info()
       << mm.info() << mr.info() << mt.info();
}
