/*! \page Phospholipid membrane and protein
 *
 * Simulate a number of flexible lipid headgroups in a salt
 * solution and a external protein.
 *
 * \author Bjoern Persson
 * \date Lund, 2009
 * \include memp.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/moves/membrane.h"

using namespace std;
using namespace Faunus;                 // Access to Faunus classes

int main() {
  // General setup
  slump slump;
  cout << faunus_splash();
  inputfile in("memp.conf");
  //atom.load(in);
  mcloop loop(in);
  canonical nvt;
  slit con(in);
  springinteraction<pot_debyehuckelXY> pot(in);

  // File I/O
  iopqr pqr;
  ioaam aam;
  ioxtc xtc(con.len);

  // Phospholipid membrane
  popscmembrane mem;
  mem.load( in , con);
  membranemove mm(nvt, con, pot, mem);
  mm.dpmm=in.getflt("monomer_dp",3);
  mm.dpgp=in.getflt("graftpoint_dp",2);
  mm.dplatt=in.getflt("latticetranslation_dp",4);

  // Protein
  macromolecule prot;
  prot.add(con,
          aam.load(in.getstr("protein")));
  prot.move(con, -prot.cm);
  prot.accept(con);
  if (in.getboo("place")==true) {
    point move;
    move.x=0, move.y=0, move.z=in.getflt("place_dp", 0);
    prot.move(con, -move);
    prot.accept(con);
  }
  zmove z(nvt, con, pot);
  z.dp=in.getflt("prot_zdp",10);
  macrorot r(nvt, con, pot);
  r.dp=in.getflt("prot_rotdp", 1);

  // Analysis
  point slask;
  radial_profile popsend(slask, 0, con.len/2);
  histogram dfpw(1, 0, con.len);

  aam.load(con,"conf.aam");
  prot.masscenter(con);
  for (int i=0; i<mem.pops.size(); i++)
    mem.pops[i].masscenter(con);
  for (int i=0; i<mem.popc.size(); i++)
    mem.popc[i].masscenter(con);
  cout << prot.info();
  cout << mem.info();
  systemenergy sys(pot.uself_popscmem(con.p, mem)
                  +pot.energy(con.p, mem, prot));

  // Input and system info

  cout << con.info() << atom.info()
       << pot.info() 
       << in.info();

  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      switch (int(slump.random_one()*3)) {
        case 0:
          sys+=mm.move(mem);
          break;
        case 1:
          sys+=z.move(prot);
         // slask=prot.cm;
          dfpw.add(prot.cm.z+con.len/2);
          break;
        case 2:
          sys+=r.move(prot);
          break;
      }
      if (slump.random_one()>0.5)
        xtc.save( "tis", con.p );

      }                                   // END of micro loop
      sys.update(pot.uself_popscmem(con.p, mem) 
                +pot.energy(con.p, mem, prot));
      dfpw.write("rdfw.dat");
      cout <<loop.timing(); 
      cout <<"#   Energy    (c, a, d)    = "<<sys.cur<<", "<<sys.uavg.avg()<<", "<<sys.cur-sys.sum<<endl;
      for (int i=0; i<mem.pops.size(); i++){
        slask=con.p[mem.pops[i].end];
        popsend.add(slask);
      }
    }                                     // END of macro loop

    // Object 

    cout <<z.info()<<r.info()<<mm.info()<<sys.info();

    aam.save("conf.aam",con.p);
    pqr.save("conf.pqr",con.p);
    popsend.write("end-dist.dat");
                                       // END of simulation

}
