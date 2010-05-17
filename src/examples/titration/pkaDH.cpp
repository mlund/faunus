#include "faunus/faunus.h"
#include "faunus/potentials/pot_debyehuckel.h"

using namespace Faunus;
using namespace std;

/*!
 * A class to measure volume of macromolecules
 * via hit and miss procedure.
 * \author Mikael Lund and Bjorn Persson
 * \date Feb 2010
 */
class measureVolume {
  private:
    unsigned long int hit,cnt;
    macromolecule g;
    hardsphere hs;
    cell con;
  public:
    measureVolume( container &c, macromolecule protein ) : con(protein.radius(c.p)+2.0) {
      hit=cnt=0;
      g.beg=0;
      for (int i=protein.beg; i<=protein.end; i++) {
        con.p.push_back( c.p[i] );
        con.trial.push_back( c.p[i] );
        g.end++;
      }
      g.masscenter(con.p);
      g.move(con, -g.cm); 
      g.accept(con);
    }
    void sample(unsigned int steps, double probeRadius=0.) {
      particle probe;
      probe.radius=probeRadius;
      for (unsigned int i=0; i<steps; i++) {
        cnt++;
        con.randompos(probe);
        if (hs.overlap(con.p, g, probe)==true)
          hit++;
      }
    }
    double getVolume() { return con.getvolume()*static_cast<double>(hit)/cnt; }
    double getRadius() { return pow( 3*getVolume() / (4*std::acos(-1)), 1./3.  ); }
};


int main(int argc, char* argv[]) {
  cout << faunus_splash();                // Faunus spam
#ifdef DHTEIXEIRA
  string config = "pkaAT.conf";           // Default input (parameter) file
#else
  string config = "pka.conf";             // Default input (parameter) file
#endif
  if (argc==2) config = argv[1];          // ..also try to get it from the command line
  inputfile in(config);                   // Read input file
  mcloop loop(in);                        // Set Markov chain loop lengths
  checkValue test(in);                    // Unit testing
  cell con(in);                           // Use a spherical simulation container
  canonical nvt;                          // Use the canonical ensemble
  interaction<pot_debyehuckel> pot(in);   // Specify pair potential
  vector <macromolecule> protein(1);      // Group for the protein
  ioaam aam;                              // Protein input file format is AAM
  protein[0].add( con,
      aam.load(in.getstr("protein")));    // Load protein structure
  protein[0].move(con, -protein[0].cm);   // ..translate it to origo (0,0,0)
  protein[0].accept(con);                 // ..accept translation
  aam.load(con, "confout.aam");           // Load old config (if present)

#ifdef DHTEIXEIRA
  ATchargereg tit(nvt,con,pot,in.getflt("pH", 7.),in,pot.pair);
  protein[0].conc = in.getflt("ProteinConc", 0.0001);

  measureVolume vol( con, protein[0] );
  cout << "# Encapsulating sphere      = " << protein[0].vradius(con.p) << endl;
  cout << "# Measuring protein radius... " << flush;
  vol.sample(1e6);
  protein[0].cm.radius = vol.getRadius() + 2.0;
  protein[0].cm_trial.radius = protein[0].cm.radius;
  cout << "done. (" << protein[0].cm.radius << " AA incl. salt)" << endl;
  systemenergy sys(tit.totalenergy(protein, con.p));
#else
  DHchargereg tit(nvt,con,pot,in.getflt("pH", 7.),in.getflt("mu_proton"));
  systemenergy sys(pot.energy(con.p));     // System energy analysis
#endif

  cout << con.info() << tit.info()         // Information
       << pot.info() << atom.info()
       << in.info();

  while ( loop.macroCnt() ) {              // Markov chain
    while ( loop.microCnt() ) {
#ifdef DHTEIXEIRA
      sys+=tit.titrateall( protein );
#else
      sys+=tit.titrateall();               // Titrate protein sites
#endif
      protein[0].charge(con.p);            // Re-calc. protein charge
      protein[0].dipole(con.p);            // Re-calc. dipole moment
    }                                      // END of micro loop
#ifdef DHTEIXEIRA
     sys.update(tit.totalenergy(protein, con.p));
#else
    sys.update(pot.energy(con.p));         // Update system energy
#endif
    aam.save("confout.aam", con.p);        // Save config. to disk
    cout << loop.timing();                 // Show progress
  }                                        // END of macro loop

  tit.check(test);
  sys.check(test);
  test.check("ProteinCharge", protein[0].Q.avg() );
  test.check("ProtienDipole", protein[0].dip.avg() );

  cout << loop.info() << sys.info()
    << tit.info() << protein[0].info(con)
    << test.report()
    << endl;                               // Print final results
}

