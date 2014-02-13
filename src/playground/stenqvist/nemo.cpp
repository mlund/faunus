#include <faunus/faunus.h>
#include <faunus/multipole.h>
using namespace Faunus;                          // use Faunus namespace
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace; 
typedef DipoleDipoleWolf TpairDDW;
typedef LennardJones TpairLJ;
typedef CombinedPairPotential<TpairLJ,TpairDDW> Tpair;

#ifdef POLARIZE
typedef Move::PolarizeMove<AtomicTranslation<Tspace> > TmoveTran; 
typedef Move::PolarizeMove<AtomicRotation<Tspace> > TmoveRot;
#else
//typedef Move::PolarizeMove<AtomicTranslation<Tspace> > TmoveTran; 
//typedef Move::PolarizeMove<AtomicRotation<Tspace> > TmoveRot;
typedef Move::AtomicTranslation<Tspace> TmoveTran;   
typedef Move::AtomicRotation<Tspace> TmoveRot;
#endif

template<class Tpairpot, class Tid>
bool savePotential(Tpairpot pot, TmoveRot rot, Tid ida, Tid idb, string file) {
  std::ofstream f(file.c_str());
  if (f) {
    double min=1.1 * (atom[ida].radius+atom[idb].radius);
    DipoleParticle a,b;
    a=atom[ida];
    b=atom[idb];
    a.mu = Point(1,0,0);
    b.mu = Point(1,0,0);
    /*f << "# Pair potential: " << pot.brief() << endl
      << "# Atoms: " << atom[ida].name << "<->" << atom[idb].name
      << endl;*/
    for (double r=min; r<=14; r+=0.05) {
      f << std::left << std::setw(10) << r << " "
        << exp(-pot(a,b,Point(r,0,0))) << endl;
    }
    return true;
  }
  return false;
}

template<class Tspace, class Tpairpot>
void getSystemEnergyExternal(Tspace &spc, Tpairpot pot, string file) {
    double E = 0;
    int N = spc.p.size();
    for(int n = 0; n < N-1; n++) {
      for(int m = n+1;m < N; m++) {
        Point r = spc.geo.vdist(spc.p[n],spc.p[m]);
        E = E + pot(spc.p[n],spc.p[m],r);
      }
    }
    cout << "System energy for " << file << "-potential is: " << E << endl;
}

template<class Tspace, class Tpairpot>
double getSystemEnergyExternalIn(Tspace &spc, Tpairpot pot) {
    double E = 0;
    int N = spc.p.size();
    for(int n = 0; n < N-1; n++) {
      for(int m = n+1;m < N; m++) {
        Point r = spc.geo.vdist(spc.p[n],spc.p[m]);
        E = E + pot(spc.p[n],spc.p[m],r);
      }
    }
    return E;
}

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
  //Analysis::LineDistribution<> rdf(0.1);       // particle-particle g(r)
  Analysis::Table2D<double,Average<double> > mucorr(0.1);       // particle-particle g(r)
  Analysis::Table2D<double,double> mucorr_distribution(0.1);
  Analysis::Table2D<double,double> mucorr_distribution1(0.1);
  Average<double> EnergyDDW;
  TmoveTran trans(in,pot,spc);
  TmoveRot rot(in,pot,spc);
  trans.setGroup(sol);                                // tells move class to act on sol group
  rot.setGroup(sol);                                  // tells move class to act on sol group
  spc.load("state");
  spc.trial = spc.p;
  UnitTest test(in);               // class for unit testing
  DipoleWRL sdp;
  FormatXTC xtc(spc.geo.len.norm());
  savePotential(TpairDDW(in),rot, atom["sol"].id, atom["sol"].id, "pot_dipdip.dat");
  Analysis::DielectricConstant diel(spc);
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );   // initial energy
  
  getSystemEnergyExternal(spc, TpairLJ(in),"InitLennardJones");
  getSystemEnergyExternal(spc, TpairDDW(in),"InitDipoleDipoleWolf");
  
  diel.sampleDP(spc);
  cout << "Diel before: " << diel.info() << endl;
  
  while ( loop.macroCnt() ) {                         // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() > 0.5)
        sys+=trans.move( sol.size() );                // translate
      else
        sys+=rot.move( sol.size() );                  // rotate

      if (slp_global()<1.2) {
        //EnergyDDW += getSystemEnergyExternalIn(spc,TpairDDW(in));
        double r, sca;
        for (auto i=sol.front(); i<=sol.back(); i++) { // salt rdf
          for (auto j=i+1.; j<=sol.back(); j++) {
            r = spc.geo.dist(spc.p[i],spc.p[j]); 
            rdf(r)++;
            sca = spc.p[i].mu.dot(spc.p[j].mu);
            mucorr(r) += sca;
            mucorr_distribution(sca) += 1.;
            mucorr_distribution1(r) += 1.5*(sca*sca -1.);
          }
        }
      }
      diel.sampleDP(spc);
      //if (slp_global()>0.99)
      //  xtc.save(textio::prefix+"out.xtc", spc.p);  
    }
    cout << "Std_eps: " << diel.getStdDiel() << endl;
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
  }
  diel.sampleKirkwood(spc);
  diel.saveKirkwood();

  //getSystemEnergyExternal(spc, TpairLJ(in),"EndLennardJones");
  getSystemEnergyExternal(spc, TpairDDW(in),"EndDipoleDipoleWolf");
  cout << "System average energy: " << EnergyDDW.avg() << endl;
  
  // perform unit tests
  trans.test(test);
  rot.test(test);
  sys.test(test);
  sdp.saveDipoleWRL("stockmayer.wrl",spc,sol);
  FormatPQR().save("confout.pqr", spc.p);
  rdf.save("gofr.dat");                               // save g(r) to disk
  mucorr.save("mucorr.dat");                               // save g(r) to disk
  mucorr_distribution.save("mucorr_distribution.dat");
  mucorr_distribution1.save("mucorr_distribution1.dat");
  std::cout << spc.info() + pot.info() + trans.info()
    + rot.info() + sys.info(); // final info
    cout << "Diel: " << diel.info() << endl;
  spc.save("state");

  return test.numFailed();
}
/**  @page example_stockmayer Example: Stockmayer potential
 *
 This will simulate a Stockmayer potential in a cubic box.

 Run this example from the `examples` directory:

 ~~~~~~~~~~~~~~~~~~~
 $ make
 $ cd src/examples
 $ ./stockmayer.run
 ~~~~~~~~~~~~~~~~~~~

 stockmayer.cpp
 ============

 @includelineno examples/stockmayer.cpp

*/
