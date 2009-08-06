/*!\page test_pka Titration-GLU3
 * Proton titration of a GLU3 porphyrin dendrimer in explicit salt.
 *
 * \author Bjorn Persson
 * \include glu3pka.cpp
 */
#include "faunus/faunus.h"
#include "faunus/energy/springinteraction.h"

using namespace Faunus;
using namespace std;

int main(int argc, char* argv[]) {
  cout << faunus_splash();             // Faunus spam
//  random spl;                           // Random numbers
  string config = "pka.conf";          // Default input (parameter) file
  if (argc==2) config = argv[1];       // ..also try to get it from the command line
  inputfile in(config);                // Read input file
  mcloop loop(in);                     // Set Markov chain loop lengths
  cell con(in);                        // Use a spherical simulation container
  canonical nvt;                       // Use the canonical ensemble
  springinteraction<pot_coulombr12> pot(in);  // Specify pair potential
  glu3 glu3(con, in);                  // Group for the protein
  ioaam aam;                           // Protein input file format is AAM
  iopqr pqr;                           // PQR coordinate output
  glu3.move(con, -glu3.cm);            // ..translate it to origo (0,0,0)
  glu3.accept(con);                    // ..accept translation
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions
  saltmove sm(nvt, con, pot);          // Class for salt movements
  sm.dp=in.getflt("dp_salt", 20);
  monomermove mm(nvt,con,pot,in);      // ...rattle
  aam.load(con, "confout.aam");        // Load old config (if present)
  widom wid(30);
  widomSW wid2(30);                    // Class for single particle insertion w. charge scaling
  wid.add( atom("NA") );
  wid.add( atom("CL") );
  wid.runfraction=0.05;
  wid2.add( atom("NA"));
  wid2.add( atom("CA"));
  wid2.add( atom("LA"));
  wid2.add( atom("CL"));
  wid2.add( atom("SO"));
  wid2.runfraction=0.05;
//  average_vec<double> pot1(100), pot2(100), pot12(100);
//  double p1,p2,p3;

  ioxtc xtc(con.r*2);                  // Gromacs xtc output (if installed)
//#ifdef GCPKA // "Grand Canonical" titration
//  HAchargereg tit(nvt,con,pot,salt,in.getflt("pH", 7.),in.getflt("catpot"));
//#else        // "Normal" titration
  glu3corechargereg tit(nvt,con,pot,in,salt);
//#endif



  for (int i =0; i<con.p.size(); i++) {
    con.p[i].radius=1.0;
    con.trial[i].radius=1.0;
  }
  vector<particle> psmear = con.p;

  systemenergy sys(pot.energy(con.p, glu3, salt)                           //The total system energy
                  +pot.uself_polymer(con.p, glu3.chains)                   //Note that it may have to redefined
                  +pot.internal(con.p, salt)                               //if the structure of GLU3 is altered
                  +pot.energy(con.p, glu3.core, glu3.chains)               //from that in GLU3simp.mol2!!!
                  +pot.energy(con.p[4], con.p[13])
                  +pot.k*pow(con.p[51].dist(con.p[54])-pot.req,2)
                  +pot.k*pow(con.p[44].dist(con.p[131])-pot.req,2)
                  +pot.k*pow(con.p[30].dist(con.p[230])-pot.req,2)
                  +pot.k*pow(con.p[37].dist(con.p[181])-pot.req,2));
             
  cout << in.info() << con.info() << tit.info()                             // Some information
       << pot.info() << atom.info();

  while ( loop.macroCnt() ) {                                               // Markov chain 
    while ( loop.microCnt() ) {
      switch (rand()%4) {                                                 // Randomly chose move
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
          sys+=mm.move(glu3.chains);
          break;
      }
      wid.insert(con,pot);
      wid2.insert(con,pot);
      if(rand()>0.9) {
//        p1=pot.potential(con.p, 4 );
//        p1=pot.potential(con.p, 13);
//      if (slp.random_one()>.99)
//        xtc.save("ignored-name.xtc", con.p);                                // Save trajectory
    }}                                                                       // END of micro loop
 
    sys.update(pot.energy(con.p, glu3, salt)
                  +pot.uself_polymer(con.p, glu3.chains) 
                  +pot.internal(con.p, salt) 
                  +pot.energy(con.p, glu3.core, glu3.chains)
                  +pot.energy(con.p[4], con.p[13])
                  +pot.k*pow(con.p[51].dist(con.p[54])-pot.req,2)
                  +pot.k*pow(con.p[44].dist(con.p[131])-pot.req,2)
                  +pot.k*pow(con.p[30].dist(con.p[230])-pot.req,2)
                  +pot.k*pow(con.p[37].dist(con.p[181])-pot.req,2));

    aam.save("confout.aam", con.p);                                         // Save config. to disk
    pqr.save("confout.pqr", con.p, tit);                                    // Save PQR file to disk - cool in VMD!
    tit.applycharges(psmear);
    aam.save("smeared.aam", psmear);
    cout << "# ENERGY "<<endl
         << "# Cur, Avg, Drift       = "<<sys.cur<<" , "<<sys.uavg.avg()<<" ("<<sys.uavg.stdev()<<") , "<<sys.cur-sys.sum<<endl;
    cout << loop.timing();                                                  // Show progress
  }                                                                         // END of macro loop
  cout << sys.info() <<loop.info() <<tit.info()<< sm.info() << mm.info() << wid.info()<<wid2.info()
       << salt.info(con)
       << glu3.info();                                                      // Print final results
  xtc.close();
}

