/*! \page example_bulk Example: Bulk
 This will simulate an arbitrary number of linear polymers in an an NPT simulation
 container with explicit salt particles and implicit solvent (dielectric continuum).
 We include the following Monte Carlo move:
 \li salt translation
 \li isobaric volume move (NPT ensemble)

 Information about the input file can be found in \c bulk.run in the \c src/examples
 directory.
*/
#include <faunus/faunus.h>
using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
typedef Potential::CoulombWolfLJ Tpairpot;// particle pair potential: primitive model
//typedef Potential::CoulombLJ Tpairpot;// particle pair potential: primitive model

namespace Faunus {
  namespace Energy {
    template<class T1, class T2>
      class NonbondedPairCut : public Energy::Nonbonded<T1,T2> {
        private:
          typedef Energy::Nonbonded<T1,T2> Tbase;
          double Rc2;
        public:
          NonbondedPairCut(InputMap &in) : Tbase(in) {
            Tbase::name+=" (cut)";
            Rc2=100;
          }
          double i2all(const p_vec &p, int i) {
            double u=0;
            int n=(int)p.size();
#pragma omp parallel for reduction (+:u)
            for (int j=0; j<n; ++j)
              if (i!=j) {
                double r2=Tbase::geometry.sqdist(p[i],p[j]);
                if (r2<Rc2)
                  u+=Tbase::pairpot(p[i],p[j],r2);
              }
            return u;
          }
          double g_internal(const p_vec &p, Group &g) { 
            double u=0;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
            for (int i=g.front(); i<g.back(); i++)
              for (int j=i+1; j<=g.back(); j++) {
                double r2=Tbase::geometry.sqdist(p[i],p[j]);
                if (r2<Rc2)
                  u+=Tbase::pairpot(p[i],p[j],r2);
              }
            return u;
          }
      };
  }//namespace
}//namespace

int main() {
  cout << textio::splash();           // show faunus banner and credits


  InputMap mcp("bulk.input");         // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  FormatPQR pqr;                      // PQR structure file I/O
  FormatAAM aam;                      // AAM structure file I/O
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);                 // class for unit testing

  // Energy functions and space
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  Space spc( pot.getGeometry() );

  // Markov moves and analysis
  Move::Isobaric iso(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Analysis::RadialDistribution<float,int> rdf_aa(0.1), rdf_bb(0.1), rdf_ab(0.1);

  // Add salt
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  mv.setGroup(salt);

  spc.load("state");                           // load old config. from disk (if any)

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global.randOne() < 0.5)
        sys+=mv.move( salt.size() );  // translate salt
      else 
        sys+=iso.move();              // isobaric volume move

      if (slp_global.randOne() < 0.1) {
        particle::Tid a=atom["Na"].id, b=atom["Cl"].id;
        for (auto i=salt.front(); i<salt.back(); i++) // salt rdf
          for (auto j=i+1; j<=salt.back(); j++) {
            if ( (spc.p[i].id==a && spc.p[j].id==b) || (spc.p[i].id==b && spc.p[j].id==a) )
              rdf_ab( spc.geo->dist(spc.p[i],spc.p[j]) )++;
          }
      }

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();

  } // end of macro loop

  // save to disk
  rdf_aa.save("bulk.rdf_aa.dat");
  rdf_bb.save("bulk.rdf_bb.dat");
  rdf_ab.save("bulk.rdf_ab.dat");
  pqr.save("confout.pqr", spc.p);
  spc.save("state");

  // perform unit tests
  iso.test(test);
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() << sys.info() << mv.info() << iso.info() << test.info();

  return test.numFailed();
}
