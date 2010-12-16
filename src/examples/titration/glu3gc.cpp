/*!\page test_pka Titration-GLU3
 * Proton titration of a GLU3 porphyrin dendrimer in explicit salt.
 *
 * \author Bjorn Persson
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"
#include "faunus/potentials/pot_coulomb.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  cout << faunus_splash();             // Faunus spam
  string config = "gcglu3pka.conf";    // Default input (parameter) file
  if (argc==2) config = argv[1];       // ..also try to get it from the command line
  // SYSTEM AND INPUT
  inputfile in(config);                // Read input file
  mcloop loop(in);                     // Set Markov chain loop lengths
#ifdef MI
  box con(in);                         // Use a cube with minimum image conv.
#else
  cell con(in);                        // Use a spherical simulation container
#endif
  grandcanonical nmt;                  // Use the canonical ensemble
  springinteraction<pot_coulombr12> pot(in);  // Specify pair potential
  // IO
  io io;
  ioaam aam;                           // Protein input file format is AAM
  iopqr pqr;                           // PQR coordinate output
  // PATRICLES
  glu3 glu3(con, in);                  // Group for the protein
  glu3.move(con, -glu3.cm);            // ..translate it to origo (0,0,0)
  glu3.accept(con);                    // ..accept translation
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions
  // MARKOVMOVES
  saltmove sm(nmt, con, pot, in);      // Class for salt movements
  saltbath sb(nmt,con, pot ,in, salt); // GC salt reservoir
  GCglu3corechargereg tit(nmt,con,pot,in);  // GC titration routine (hardcoded for the core)
  monomermove mm(nmt,con,pot,in);      // ...rattle
  crankShaft cs(nmt,con,pot,in);       // Crankshaft
  // RESTART?
  if(nmt.load(con, "gcgroup.conf")==true)
    aam.load(con, "confout.aam");        // Load old config (if present)
  widomSW wid2(30);                    // Class for single particle insertion w. charge scaling
  wid2.add( atom("NA"));
  wid2.add( atom("CA"));
  wid2.add( atom("LA"));
  wid2.add( atom("CL"));
  wid2.add( atom("SO"));
  wid2.runfraction=0.05;

#ifdef MI
  ioxtc xtc(con.len);
#else
  ioxtc xtc(con.r*2);                  // Gromacs xtc output (if installed)
#endif

  vector<particle> psmear = con.p;

  systemenergy sys(pot.energy(con.p, glu3, salt)                           //The total system energy
                  +pot.uself_polymer(con.p, glu3.chains)                   //Note that it may have to redefined
                  +pot.internal(con.p, salt)                               //if the structure of GLU3 is altered
                  +pot.energy(con.p, glu3.core, glu3.chains)               //from that in GLU3coarse.mol2!!!
                  +pot.energy(con.p[4], con.p[13])
                  +pot.pair.k*pow(con.p[51].dist(con.p[53])-pot.pair.req,2)
                  +pot.pair.k*pow(con.p[44].dist(con.p[96])-pot.pair.req,2)
                  +pot.pair.k*pow(con.p[30].dist(con.p[151])-pot.pair.req,2)
                  +pot.pair.k*pow(con.p[37].dist(con.p[124])-pot.pair.req,2));
             
//  cout << "#Self-int = "<<pot.energy(con.p[4], con.p[13])<<endl;

  cout << in.info() << con.info() << tit.info()                             // Some information
       << pot.info() << atom.info();

  while ( loop.macroCnt() ) {                                               // Markov chain 
    while ( loop.microCnt() ) {
      switch (rand()%6) {                                                 // Randomly chose move
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
          sys+=cs.move(glu3.chains,in.getint("crankshaft_num",1));
          break;
      }
      wid2.insert(con,pot);
      if(rand()>0.9) {
    }}                                                                       // END of micro loop
//  cout << "#Self-int = "<<pot.energy(con.p[4], con.p[13])<<endl;
 
    sys.update(pot.energy(con.p, glu3, salt)
                  +pot.uself_polymer(con.p, glu3.chains) 
                  +pot.internal(con.p, salt) 
                  +pot.energy(con.p, glu3.core, glu3.chains)
                  +pot.energy(con.p[4], con.p[13])
                  +pot.pair.k*pow(con.p[51].dist(con.p[53])-pot.pair.req,2)
                  +pot.pair.k*pow(con.p[44].dist(con.p[96])-pot.pair.req,2)
                  +pot.pair.k*pow(con.p[30].dist(con.p[151])-pot.pair.req,2)
                  +pot.pair.k*pow(con.p[37].dist(con.p[124])-pot.pair.req,2));

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
  cout << sys.info() << con.info() <<loop.info() <<tit.info()<< sm.info() << sb.info() << mm.info() <<wid2.info()
       << cs.info()
       << salt.info(con)
       << glu3.info();                                                      // Print final results
  xtc.close();
}

