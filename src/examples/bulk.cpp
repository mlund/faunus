#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Geometry::Cuboid Tgeometry;   // geometry: cube w. periodic boundaries
typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot; // pair potential

/*
template<class Tpairpot, class Tgeometry>
class NonbondedEwald : public Energy::Energybase {
  private:
    Ewald<double> ew;
    Space* _spc;
  public:
    Tgeometry geo;
    Tpairpot pair;

    double p2p(const particle &a, const particle &b) FOVERRIDE {
      double r2=geo.sqdist(a,b);
      return ew.rSpaceEnergy(a.charge*b.charge, sqrt(r2)) + pair(a,b,r2);
    }

    double i2all(const p_vec &p, int i) FOVERRIDE {
      if (&p==&(spc->trial)) {
        ew.store();
        eq.kSpaceUpdate(p);
      }

      assert(i>=0 && i<int(p.size()) && "index i outside particle vector");
      double u=0;
      int n=(int)p.size();
      for (int j=0; j<i; ++j)
        u+=pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
      for (int j=i+1; j<n; ++j)
        u+=pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
      return u;
    }

};*/

int main() {
  cout << textio::splash();           // show faunus banner and credits


  InputMap mcp("bulk.input");         // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);                 // class for unit testing

  // Energy functions and space
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  Space spc( pot.getGeometry() );

  /*
  Ewald<double> ew(mcp);
  ew.store();
  ew.restore();
  ew.rSpaceEnergy(spc.p[0].charge*spc.p[1].charge, 2.0);
*/

  // Markov moves and analysis
  Move::Isobaric iso(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Analysis::RadialDistribution<> rdf_ab(0.1);      // 0.1 angstrom resolution

  // Add salt
  Group salt;
  salt.addParticles(spc, mcp);
  mv.setGroup(salt);

  spc.load("state");                               // load old config. from disk (if any)
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );// store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() < 0.5)
        sys+=mv.move( salt.size() );  // translate salt
      else 
        sys+=iso.move();              // isobaric volume move

      if (slp_global() < 0.05) {
        particle::Tid a=atom["Na"].id, b=atom["Cl"].id;
        for (auto i=salt.front(); i<salt.back(); i++) // salt radial distribution function
          for (auto j=i+1; j<=salt.back(); j++)
            if ( (spc.p[i].id==a && spc.p[j].id==b) || (spc.p[i].id==b && spc.p[j].id==a) )
              rdf_ab( spc.geo->dist(spc.p[i],spc.p[j]) )++;
      }
    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();

  } // end of macro loop

  // save to disk
  FormatPQR().save("confout.pqr", spc.p); // final PQR snapshot for VMD etc.
  rdf_ab.save("rdf.dat");         // g(r) - not normalized!
  spc.save("state");              // final simulation state

  // perform unit tests (irrelevant for the simulation)
  iso.test(test);
  mv.test(test);
  sys.test(test);
  nonbonded->pairpot.test(test);

  // print information
  cout << loop.info() << sys.info() << mv.info() << iso.info() << test.info();

  return test.numFailed();
}
/**
  @page example_bulk Example: Melted NaCl

  In this example we simulate melted NaCl in the NVT and NPT ensemble. We use a
  Lennard-Jones potential combined with a shifted Coulombic potential according to
  Wolf. This gives essentially identical results to the more elaborate Particle
  Mesh Ewald method - see figure below. In contrast, using the simple minimum image
  approach with a cubic cutoff, the system freezes.
  The `bulk.cpp` program can be used to simulate any atomic mixtures and the
  dielectric constant may also be varied (it is unity in this example).

  We have the following MC moves:
  - salt translation
  - isotropic volume move (NPT ensemble)

  Information about the input file can be found in `src/examples/bulk.run`.

  ![Na-Cl distribution function with various electrostatic potentials.](wolf.png)

  bulk.cpp
  ========
  \includelineno examples/bulk.cpp

*/
