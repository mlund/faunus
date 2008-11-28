/*!\page test_spce SPCE/E
 *
 * \author Mikael Lund and Bjorn Persson
 * \include spce.cpp
 */
#include "faunus/faunus.h"
#include "faunus/potentials/pot_rfield.h"

using namespace Faunus;
using namespace std;

class ionpmf {
  private:
    slump slump;
    particle a,b;
  public:
    double dr,rmax;
    vector< average<double> > pmf,u1,u2;
    ionpmf() {
      rmax=11.;
      dr=0.1;
      a.charge=0.1;
      b.charge=-a.charge;
      pmf.resize( int(rmax/dr+.5) );
      u1.resize( int(rmax/dr+.5) );
      u2.resize( int(rmax/dr+.5) );
    }
    void insert(container &c, energybase &pot, double r=5.0) {
      double u=0,ut1,ut2;
      a.x=2.*r*slump.random_half()+0.000001;
      b.x=-a.x;
      r=abs(a.x-b.x);
      ut1=a.charge*pot.potential(c.p, a);
      ut2=b.charge*pot.potential(c.p, b);
      u=ut1+ut2;
      pmf[ int( r/dr+.5) ] +=exp(-u);
      u1[ int( r/dr+.5) ] +=exp(-ut1);
      u2[ int( r/dr+.5) ] +=exp(-ut2);
    }
};

int main(int argc, char* argv[]) {
  slump slump;
  cout << faunus_splash();             // Faunus spam
  inputfile in("spce.conf");           // Read input file
  cell con(in);                        // Use a spherical simulation container
  ioaam aam(con.atom);                 // Protein input file format is AAM
  iogro gro(con.atom,in);
  mcloop loop(in);                     // Set Markov chain loop lengths
  canonical nvt;                       // Use the canonical ensemble
  interaction<pot_test> pot(in);       // Specify pair potential
  pot_rfield rfield(in);
  pot.pair.init(con.atom);
  ionpmf ip;

  macrorot mr(nvt, con, pot);          // Class for molecular rotation
  translate mt(nvt, con, pot);         // Class for molecular translation

  // Load central protein
  macromolecule protein;               // Group for the protein
  protein.add( con,
      aam.load(in.getstr("protein"))); // Load protein structure
  protein.move(con, -protein.cm);      // ..translate it to origo (0,0,0)
  protein.accept(con);                 // ..accept translation

  // Add salt
  salt salt;                           // Group for salt and counter ions
  salt.add( con, in );                 //   Insert sodium ions
  saltmove sm(nvt, con, pot);          // Class for salt movements

  // Add water
  molecules sol(3);                    // We want a three point water model
  sol.name="SPC/E Solvent";
  vector<particle>
    water = aam.load("water.aam");     // Load bulk water from disk (typically from MD)
  con.atom.reset_properties(water);    // Set particle parameters according to Faunus
  sol.add(con,water,sol.numatom);      // Inject water into the cell - avoid salt and protein overlap
  water.clear();                       // Free the (large) bulk water reservoir

  // Measure electric potential at all protein atoms
  pointpotential ppot;
  for (int i=protein.beg; i<=protein.end; i++)
    ppot.add(con.p[i], con.atom[con.p[i].id].name );

  // Distribution functions
  FAUrdf spccell(float(0.2), float(50.));
  FAUrdf nacell(float(0.2), float(50.));
  FAUrdf saltrdf(con.atom["NA"].id,con.atom["CL"].id,0.2,20.);
  FAUrdf catcat( con.atom["NA"].id,con.atom["NA"].id,0.2,20.);
  FAUrdf spcrdf( con.atom["OW"].id,con.atom["OW"].id,0.2,10.);

  histogram radorder(float(0.2), float(0.), float(con.r));

  mr.dp=0.6;
  mt.dp=0.6; //0.6
  sm.dp=0.4; //0.4

  aam.load(con, "confout.aam");        // Load old config (if present)

  systemenergy sys( 
      pot.internal(con.p, sol, sol.numatom) +
      //pot.internal(con.p, salt) +
      //pot.energy(con.p, protein, salt) +
      pot.energy(con.p, sol)
      );

  cout << in.info() << salt.info()
    << con.info() << con.atom.info()
    << pot.info() << sol.info() << rfield.info();

  group head;                            // Group of first three particles [0:2]
  head.set(0,2);                         // (Used to swap w. SPC for parallization reasons)
  macromolecule m;
  point origo;
  point trace = con.p[sol.beg];
  int cnt=0;
  while ( loop.macroCnt() ) {            // Markov chain 
    while ( loop.microCnt() ) {
      m=sol[ sol.random() ];
      m.cm=m.cm_trial=con.p[m.beg];
      switch (rand() % 3) {              // Randomly choose move
        case 0:
         // sys+=sm.move(salt,1);          // Displace a salt particle
          break;
        case 1:
          sys+=mr.move(m);             // Rotate solvent
          cnt++;
          break;
        case 2:
          sys+=mt.move(m);             // Translate solvent
          cnt++;
          break;
      }

      if (slump.random_one()<0.1) {
        ppot.sample(con,pot);
        ip.insert(con,pot);
      }
      if (slump.random_one()<0.05) {
        saltrdf.update(con, salt);
        catcat.update(con, salt);
        nacell.update(con, origo, "NA");
        int st=sol.size()/3;
        double l;
        //point mu;
        for (int i=0; i<st; i++) {
          m=sol[i];
          m.cm=m.masscenter(con.p);
          mulen=m.dipole(con.p);
          radorder.add(m.cm.len(), m.cm*m.mu/(m.cm.len()*m.mulen) );
        }
      }
      if (slump.random_one()<0.03) {
        spccell.update(con, origo, "OW");
      }
    }//end of micro-loop 

    spcrdf.update(con, sol);
    nacell.write("rdf-cell-NA.dat");
    spccell.write("rdf-cell-OW.dat");
    spcrdf.write("rdf-OW-OW.dat");
    catcat.write("rdf-Na-Na.dat");
    saltrdf.write("rdf-Na-Cl.dat");
    radorder.write("mu-OW.dat");

    sys.update( 
        pot.internal(con.p, sol, sol.numatom) +
        //pot.internal(con.p, salt) +
        //pot.energy(con.p, protein, salt) +
        pot.energy(con.p, sol)
        );

    gro.save("confout.gro", con.p);      // Save config. to disk
    aam.save("confout.aam", con.p);      // Save config. to disk
    cout << loop.timing()                // Show progress
         << "#   Energy (kT): sum average drift = "
         << sys.sum << " " << sys.uavg.avg() << " "
         << std::abs(sys.cur-sys.sum) << endl;
    cout << ppot.info();
    cout << "#   First water oxygen movement = " << trace.dist(con.p[sol.beg]) << endl;

  }//end of macro-loop
  cout << sys.info()
       << salt.info(con)
       << protein.info() << loop.info()  // Print final results
       << sm.info() << mr.info() << mt.info();

  // Print dipole intertion results
  particle a,b;
  a.charge=0.1;
  b.charge=-a.charge;
  for (float r=0; r<10.; r+=0.1) {
    a.x=.5*r;
    b.x=-a.x;
    cout
      << r << " "
      << -log(ip.pmf[ int(r/ip.dr+.5) ].avg()) << " "
      << -log(ip.u1[ int(r/ip.dr+.5) ].avg()) << " "
      << -log(ip.u2[ int(r/ip.dr+.5) ].avg()) << " ";
    double u = rfield.f * (
          rfield.selfenergy(a) + rfield.selfenergy(b)
          + rfield.pairpot(a,b) );
    cout << u << " " << u - rfield.f*a.charge*b.charge/a.dist(b) << endl;
  }

}

