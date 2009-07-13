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
#include "faunus/energy/springinteractionexfield.h"
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
  springinteractionfield<pot_debyehuckelXYcc> pot(in);

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
  prot.masscenter(con);
  if (in.getboo("place")==true) {
    point move;
    move.x=0, move.y=0, move.z=in.getflt("place_dp", 0);
    prot.move(con, -move);
    prot.accept(con);
  }
  zmove z(nvt, con, pot);
  z.zmax=100;
  z.dp=in.getflt("prot_zdp",10);
  macrorot r(nvt, con, pot);
  r.dp=in.getflt("prot_rotdp", 1);

  // Analysis
  radial_profile popsend( 0, con.zlen/2, 1);
  histogram dfpw(1, 0, con.zlen);
  histogram memz(0.1, 0, con.zlen);
  distributions dist(1.0, 0, con.zlen);

  aam.load(con,"conf.aam");
  prot.masscenter(con);
  for (int i=0; i<mem.pops.size(); i++)
    mem.pops[i].masscenter(con);
  for (int i=0; i<mem.popc.size(); i++)
    mem.popc[i].masscenter(con);
  cout << prot.info();
  cout << mem.info();
  
  pot.pair.setchargedens(mem.charge(con.p)/con.xyarea);
  systemenergy sys(pot.uself_popscmem(con.p, mem)
                  +pot.energy(con.p, mem, prot)
                  +pot.extpot(con.p, prot)
                  +pot.hydrophobic(con.p, prot));

  // Input and system info

  cout << con.info() << atom.info()
       << pot.info() 
       << in.info();

  // Switch parameters
  double sum, tr, rr, mr, randy;
  tr=in.getflt("tr", 1), rr=in.getflt("rr", 1), mr=in.getflt("mr",1);
  sum=tr+rr+mr, tr/=sum, rr/=sum, mr/=sum;

  // Help variables
  point slask;

  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      randy=slump.random_one();
      if(randy<mr){
        sys+=mm.move(mem);
        if (prot.cm.z+con.len*0.5<50.)
          for (int i=0; i<mem.pops.size(); i++){
            slask=con.p[mem.pops[i].end];
            popsend.add(prot.cm, slask);
          }
      }
      if(randy<tr+mr && randy > mr) {
        sys+=z.move(prot);
        if(prot.cm.z<con.zlen-130)  //avoid sampling the hard wall -> +z
          dfpw.add(prot.cm.z+con.zlen/2);
      }
      if(randy<tr+mr+rr && randy>tr+mr) 
          sys+=r.move(prot);
      if (slump.random_one()>0.9 && in.getboo("movie",false)==true)
        xtc.save( "tis", con.p );

        dist.add("Prot-Membrane energy", prot.cm.z+con.zlen/2., pot.energy(con.p, mem, prot)+pot.extpot(con.p, prot)+pot.hydrophobic(con.p, prot)); 
        dist.add("Prot-External Correction", prot.cm.z+con.zlen*0.5, pot.extpot(con.p, prot));
        dist.add("Prot-Membrane LJ-energy", prot.cm.z+con.zlen*0.5, pot.ljenergy(con.p, prot));
        dist.add("Prot-Membrane HP-energy", prot.cm.z+con.zlen*0.5, pot.hydrophobic(con.p, prot));
      }                                   // END of micro loop

      sys.update(pot.uself_popscmem(con.p, mem) 
                +pot.energy(con.p, mem, prot)
                +pot.extpot(con.p, prot)
                +pot.hydrophobic(con.p, prot));

      dist.write("dist.dat");
      dfpw.write("rdfw.dat");
      cout <<loop.timing(); 
      cout <<"#   Energy    (c, a, d)    = "<<sys.cur<<", "<<sys.uavg.avg()<<", "<<sys.cur-sys.sum<<endl;
      for(int i=mem.beg; i<=mem.end; i++)
        memz.add(con.p[i].z+con.zlen*0.5);
    }                                     // END of macro loop

    // Object 

    cout <<loop.info()<<z.info()<<r.info()<<mm.info()<<sys.info();

    aam.save("conf.aam",con.p);
    pqr.save("conf.pqr",con.p);
    popsend.write("end-dist.dat");
    memz.write("mem-densdist.dat");
                                       // END of simulation

}
