#include <faunus/faunus.h>
#include <faunus/multipole.h>
using namespace Faunus;                          // use Faunus namespace
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Geometry::Cuboid Tgeo;                   // select simulation geometry and pair potential
//typedef DipoleDipole Tpair;
//typedef CombinedPairPotential<CoulombHS,DipoleDipole> Tpair;
typedef CombinedPairPotential<HardSphere,IonDipole> Tpair;
//typedef CoulombLJ Tpair;
//typedef DipoleDipoleRF Tpair;




int main() {
  ::atom.includefile("stockmayer.json");         // load atom properties
  InputMap in("stockmayer.input");               // open parameter file for user input
  Energy::NonbondedVector<Tpair,Tgeo> pot(in);   // create Hamiltonian, non-bonded only
  EnergyDrift sys;                               // class for tracking system energy drifts
  Space spc( pot.getGeometry() );                // create simulation space, particles etc.
  GroupAtomic sol(spc, in);                      // group for particles
  MCLoop loop(in);                               // class for mc loop counting
  Analysis::RadialDistribution<> rdf(0.2);       // particle-particle g(r)
  Analysis::Table2D<double,Average<double> > mucorr(0.2);       // particle-particle g(r)
  Move::AtomicTranslation trans(in, pot, spc);   // particle move class
  Move::AtomicRotation rot(in, pot, spc);        // particle move class
  //PolarizeMove<AtomicTranslation> trans(in,pot,spc);
  //PolarizeMove<AtomicRotation> rot(in,pot,spc);
  trans.setGroup(sol);                                // tells move class to act on sol group
  rot.setGroup(sol);                                  // tells move class to act on sol group
  //spc.load("state_ST");
  //Analysis::getDielConst cdc(in.get<double>("dipdip_cutoff",pc::infty));   // 
  

  spc.p = spc.trial;
  
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy
  while ( loop.macroCnt() ) {                         // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() > 0.5)
        sys+=trans.move( sol.size() );                     // translate
      else
        sys+=rot.move( sol.size() );                       // rotate

      //cdc.sample(spc.p,*spc.geo);
      
      if (slp_global()<0.5)
        for (auto i=sol.front(); i<sol.back(); i++) { // salt radial distribution function
          for (auto j=i+1; j<=sol.back(); j++) {
            double r =pot.geometry.dist(spc.p[i],spc.p[j]); 
            rdf(r)++;
            mucorr(r) += spc.p[i].mu.dot(spc.p[j].mu)/(spc.p[i].muscalar*spc.p[j].muscalar);
          }
        }
    }    
    
    //cout << cdc.info();
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
  }
  //std::cout << " Mu_1: " << spc.p[0].mu.transpose()*spc.p[0].muscalar << " eÅ\n Mu_2: " << spc.p[1].mu.transpose()*spc.p[1].muscalar << " eÅ \n";

  FormatPQR().save("confout.pqr", spc.p);
  rdf.save("gofr.dat");                               // save g(r) to disk
  mucorr.save("mucorr.dat");                               // save g(r) to disk
  std::cout << spc.info() + pot.info() + trans.info() + rot.info() + sys.info(); // final info
  spc.save("state_ST");

}
