#include <faunus/faunus.h>
#include <faunus/multipole.h>
using namespace Faunus;                          // use Faunus namespace
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace;

#ifdef POLARIZE
typedef Move::PolarizeMove<AtomicTranslation<Tspace> > TmoveTran;
typedef Move::PolarizeMove<AtomicRotation<Tspace> > TmoveRot;
#else
typedef Move::AtomicTranslation<Tspace> TmoveTran;   
typedef Move::AtomicRotation<Tspace> TmoveRot;
#endif

typedef CombinedPairPotential<LennardJones,DipoleDipole> Tpair;

int main() {
  ::atom.includefile("stockmayer.json");         // load atom properties
  InputMap in("stockmayer.input");               // open parameter file for user input
  Energy::NonbondedVector<Tspace,Tpair> pot(in); // non-bonded only
  EnergyDrift sys;                               // class for tracking system energy drifts
  Tspace spc(in);                // create simulation space, particles etc.
  Group sol;
  sol.addParticles(spc, in);                     // group for particles
  MCLoop loop(in);                               // class for mc loop counting
  Analysis::RadialDistribution<> rdf(0.1);       // particle-particle g(r)
  Analysis::Table2D<double,Average<double> > mucorr(0.2);       // particle-particle g(r)
  TmoveTran trans(in,pot,spc);
  TmoveRot rot(in,pot,spc);
  trans.setGroup(sol);                                // tells move class to act on sol group
  rot.setGroup(sol);                                  // tells move class to act on sol group
  spc.load("state_ST");
  spc.p = spc.trial;
  
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );   // initial energy
  while ( loop.macroCnt() ) {                         // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() > 0.5)
        sys+=trans.move( sol.size() );                // translate
      else
        sys+=rot.move( sol.size() );                  // rotate

      if (slp_global()<0.5)
        for (auto i=sol.front(); i<sol.back(); i++) { // salt rdf
          for (auto j=i+1; j<=sol.back(); j++) {
            double r =spc.geo.dist(spc.p[i],spc.p[j]); 
            rdf(r)++;
            mucorr(r) += spc.p[i].mu.dot(spc.p[j].mu)
              / (spc.p[i].muscalar*spc.p[j].muscalar);
          }
        }
    }    
    
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
  }

  FormatPQR().save("confout.pqr", spc.p);
  rdf.save("gofr.dat");                               // save g(r) to disk
  mucorr.save("mucorr.dat");                               // save g(r) to disk
  std::cout << spc.info() + pot.info() + trans.info()
    + rot.info() + sys.info(); // final info
  spc.save("state_ST");
}
