/*!\page test_pka Titration-GLU3
 * Proton titration of a GLU3 porphyrin dendrimer in explicit salt in the presence of
 * a phospholipid wall. 
 *
 * \author Bjorn Persson
 * 
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/potentials/pot_minimage.h"
#include "faunus/moves/membrane.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  cout << faunus_splash();             // Faunus spam
  string config = "gcglu3pka.conf";    // Default input (parameter) file
  if (argc==2) config = argv[1];       // ..also try to get it from the command line
  // SYSTEM AND INPUT
  inputfile in(config);                // Read input file
  mcloop loop(in);                     // Set Markov chain loop lengths
  slit con(in);                        // Use a cube with minimum image conv. in XY-direcitons
  grandcanonical nmt;                  // Use the canonical ensemble
  springinteraction<pot_r12minimageXY> pot(in);  // Specify pair potential
  // IO
  io io;
  ioaam aam;                           // Protein input file format is AAM
  iopqr pqr;                           // PQR coordinate output
  // PATRICLES
  glu3 glu3(con, in);                  // Group for the protein
  point pglu3;
  pglu3=-glu3.cm, pglu3.z-=(con.zlen*0.5-in.getflt("offset",35));
  glu3.move(con, pglu3);               // ..translate it
  glu3.accept(con);                    // ..accept translation
  popscmembrane mem;
  mem.load( in , con);
  membranemove memm(nmt, con, pot, mem);
  memm.dpmm=in.getflt("monomer_dp",3);
  memm.dpgp=in.getflt("graftpoint_dp",2);
  memm.dplatt=in.getflt("latticetranslation_dp",4);
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions
  // MARKOVMOVES
  saltmove sm(nmt, con, pot, in);      // Class for salt movements
  saltbath sb(nmt,con, pot ,in, salt); // GC salt reservoir
  GCglu3corechargereg tit(nmt,con,pot,in);  // GC titration routine (hardcoded for the core)
  monomermove mm(nmt,con,pot,in);      // ...rattle
  crankShaft cs(nmt,con,pot,in);       // Crankshaft
  macrorot r(nmt, con, pot);
  r.dp=in.getflt("glu3_rotdp", 1);
  zmove z(nmt, con, pot);
  z.zmax=0;
  z.dp=in.getflt("glu3_zdp",10);
  // RESTART?
  if(nmt.load(con, "gcgroup.conf")==true)
    aam.load(con, "confout.aam");        // Load old config (if present)
  for (int i=0; i<mem.pops.size(); i++)
    mem.pops[i].masscenter(con);
  for (int i=0; i<mem.popc.size(); i++)
    mem.popc[i].masscenter(con);


  widomSW wid2(30);                    // Class for single particle insertion w. charge scaling
  wid2.add( atom("NA"));
  wid2.add( atom("CA"));
  wid2.add( atom("LA"));
  wid2.add( atom("CL"));
  wid2.add( atom("SO"));
  wid2.runfraction=0.05;
  distributions dist(1.0, 0, con.zlen);
  histogram dfdw(1, 0, con.zlen);

  ioxtc xtc(con.len);

  vector<particle> psmear = con.p;

  systemenergy sys(pot.energy(con.p, glu3, salt)                           //The total system energy
                  +pot.uself_polymer(con.p, glu3.chains)                   //Note that it may have to redefined
                  +pot.internal(con.p, salt)                               //if the structure of GLU3 is altered
                  +pot.energy(con.p, glu3.core, glu3.chains)               //from that in GLU3coarse.mol2!!!
                  +pot.energy(con.p[4], con.p[13])
                  +pot.k*pow(con.p[51].dist(con.p[53])-pot.req,2)
                  +pot.k*pow(con.p[44].dist(con.p[96])-pot.req,2)
                  +pot.k*pow(con.p[30].dist(con.p[151])-pot.req,2)
                  +pot.k*pow(con.p[37].dist(con.p[124])-pot.req,2)
                  +pot.uself_popscmem(con.p, mem)
                  +pot.energy(con.p, mem, glu3)
                  +pot.energy(con.p, mem, salt));
             
//  cout << "#Self-int = "<<pot.energy(con.p[4], con.p[13])<<endl;

  cout << in.info() << con.info() << tit.info()                             // Some information
       << pot.info() << atom.info();

  while ( loop.macroCnt() ) {                                               // Markov chain 
    while ( loop.microCnt() ) {
      switch (rand()%8) {                                                 // Randomly chose move
        case 0:
          sys+=sm.move(salt);                                               // Displace salt particles
          break;
        case 1:
          sys+=tit.titrateall();                                            // Titrate protein sites
          glu3.charge(con.p);                                               // Re-calc. dendrimer(tot) charge
          glu3.dipole(con.p);                                               // Re-calc. dendrimer(tot) dipole moment
          break;
        case 3:
          sys+=tit.move(glu3);
          glu3.core.charge(con.p);                                          // Re-calc. porphoryn charge
          break;
        case 2:
          for (int i=0; i<glu3.chains.size(); i++)
            sys+=mm.move(glu3.chains);
          break;
        case 4:
          sys+=sb.move();
          break;
        case 5:
          sys+=memm.move(mem);
          //sys+=cs.move(glu3.chains,in.getint("crankshaft_num",1));
          break;
        case 6:
          sys+=r.move(glu3);
          break;
        case 7:
          sys+=z.move(glu3);
          break;
      }
      wid2.insert(con,pot);
      if(rand()>0.9) {
        glu3.masscenter(con);
        dfdw.add(glu3.cm.z+con.zlen/2);
        dist.add("Core charge", glu3.cm.z+con.zlen/2,con.p[4].charge+con.p[13].charge);
        dist.add("Glu3 charge", glu3.cm.z+con.zlen/2,glu3.getcharge(con.p));
        dist.add("Glu3 radius of gyration", glu3.cm.z+con.zlen/2,glu3.gradius(con.p));
        for(int i=0; i<10; i++) {
          con.randompos(pglu3);
          dist.add("Electric Potential", pglu3.z+con.zlen/2,pot.potential(con.p,pglu3));
        }
      }
    }                                                                       // END of micro loop
//  cout << "#Self-int = "<<pot.energy(con.p[4], con.p[13])<<endl;
    dfdw.write("glu3wall-dist.dat");
    dist.write("distributions.dat");

    sys.update(pot.energy(con.p, glu3, salt)
                  +pot.uself_polymer(con.p, glu3.chains) 
                  +pot.internal(con.p, salt) 
                  +pot.energy(con.p, glu3.core, glu3.chains)
                  +pot.energy(con.p[4], con.p[13])
                  +pot.k*pow(con.p[51].dist(con.p[53])-pot.req,2)
                  +pot.k*pow(con.p[44].dist(con.p[96])-pot.req,2)
                  +pot.k*pow(con.p[30].dist(con.p[151])-pot.req,2)
                  +pot.k*pow(con.p[37].dist(con.p[124])-pot.req,2)
                  // salt + glu3 OK
                  +pot.uself_popscmem(con.p, mem)
                  +pot.energy(con.p, mem, glu3)
                  +pot.energy(con.p, mem, salt)
                  //+pot.extpot(con.p, prot)
                  //+pot.hydrophobic(con.p, prot)
                   );

    aam.save("confout.aam", con.p);                                         // Save config. to disk
    pqr.save("confout.pqr", con.p);                                    // Save PQR file to disk - cool in VMD!
    io.writefile("gcgroup.conf", nmt.print());

    //tit.applycharges(psmear);
    aam.save("smeared.aam", psmear);
    cout << "# ENERGY "<<endl
         << "# Cur, Avg, Drift       = "<<sys.cur<<" , "<<sys.uavg.avg()<<" ("<<sys.uavg.stdev()<<") , "<<sys.cur-sys.sum<<endl
         << "# CHARGES (dendrimer):  Total = "<<glu3.Q.avg()<<"    Core = "<<glu3.core.Q.avg()<<endl;   
    cout << loop.timing();                                                  // Show progress
  }                                                                         // END of macro loop
  cout << sys.info() << con.info() <<loop.info() <<tit.info()<< sm.info() << sb.info() 
       << mm.info() <<z.info() <<r.info() <<memm.info() <<wid2.info()
       //<< cs.info()
       << salt.info(con)<<mem.info()
       << glu3.info();                                                      // Print final results
  xtc.close();
}

