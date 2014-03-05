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
    for (double r=min; r<=14; r+=0.05) {
      f << std::left << std::setw(10) << r << " "
        << pot(a,b,Point(r,0,0)) << endl;
    }
    return true;
  }
  return false;
}

template<class Tpairpot, class Tid>
void getKeesom(Tpairpot pot, TmoveRot rot, Tid ida, Tid idb, string file) {
  std::ofstream f(file.c_str());
  if (f) {
    double min=1.5 * (atom[ida].radius+atom[idb].radius);
    DipoleParticle a,b;
    a=atom[ida];
    b=atom[idb];
    a.mu = Point(1,0,0);
    b.mu = Point(1,0,0);
    
    double step = 1;
    for (double r=min; r<=4.5; r+=0.1) {
      double Ua = 0;
      double cnt = 0;
      for(double theta1 = 0; theta1 < pc::pi; theta1+=step) {
        for(double phi1 = 0; phi1 < 2*pc::pi; phi1+=step) {
          a.mu = Point(sin(theta1)*cos(phi1),sin(theta1)*sin(phi1),cos(theta1));
          for(double theta2 = 0; theta2 < pc::pi; theta2+=step) {
            for(double phi2 = 0; phi2 < 2*pc::pi; phi2+=step) {
              b.mu = Point(sin(theta2)*cos(phi2),sin(theta2)*sin(phi2),cos(theta2));
              Ua += exp(-pot(a,b,Point(r,0,0)));
              cnt++;
            }
          }
        }
      }
      //cout << "r: " << r << endl;
      f << std::left << std::setw(10) << r << " "
        << -log(Ua/cnt) << endl;
    }
  }
}

int runSim(string name) {
  ::atom.includefile("stockmayer.json");         // load atom properties
  InputMap in(name);               // open parameter file for user input
  Energy::NonbondedVector<Tspace,Tpair> pot(in); // non-bonded only
  Energy::NonbondedVector<Tspace,TpairDDW> potDDW(in);
  Energy::NonbondedVector<Tspace,TpairLJ> potLJ(in);
  EnergyDrift sys;                               // class for tracking system energy drifts
  Tspace spc(in);                // create simulation space, particles etc.
  Group sol;
  sol.addParticles(spc, in);                     // group for particles
  MCLoop loop(in);                               // class for mc loop counting
  TmoveTran trans(in,pot,spc);
  TmoveRot rot(in,pot,spc);
  trans.setGroup(sol);                                // tells move class to act on sol group
  rot.setGroup(sol);                                  // tells move class to act on sol group
  spc.load("state");
  spc.trial = spc.p;
  UnitTest test(in);               // class for unit testing
  
  Analysis::DipoleAnalysis dian(spc);
  Average<double> EnergyDDW;
  Average<double> EnergyLJ;
  FormatXTC xtc(spc.geo.len.norm());
  DipoleWRL sdp;
  
  savePotential(TpairDDW(in),rot, atom["sol"].id, atom["sol"].id, "pot_dipdip.dat");
  //getKeesom(TpairDDW(in),rot, atom["sol"].id, atom["sol"].id, "potKeesom.dat");
  
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );   // initial energy
  EnergyDDW += Energy::systemEnergy(spc,potDDW,spc.p);
  EnergyLJ += Energy::systemEnergy(spc,potLJ,spc.p);
  
  while ( loop.macroCnt() ) {                         // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() > 0.5)
        sys+=trans.move( sol.size() );                // translate
      else
        sys+=rot.move( sol.size() );                  // rotate

      if (slp_global() < 0.2) {
        EnergyDDW += Energy::systemEnergy(spc,potDDW,spc.p);
        EnergyLJ += Energy::systemEnergy(spc,potLJ,spc.p);
        dian.sampleMuCorrelationAndKirkwood(spc);
        xtc.save(textio::prefix+"xtcout.xtc", spc.p); 
      }
      dian.sampleDP(spc);
    }
    dian.save(std::to_string(loop.count()));
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
  }
  
  // Perform unit tests
  trans.test(test);
  rot.test(test);
  sys.test(test);
  
  // Show info
  std::cout << spc.info() + pot.info() + trans.info()
    + rot.info() + sys.info() + dian.info();
    
  cout << "\n  System DipoleDipoleWolf average energy: " << EnergyDDW.avg() << " kT" << endl;
  cout << "  System LennardJones average energy:     " << EnergyLJ.avg() << " kT" << "\n" << endl;
    
  // Save data
  dian.save();
  sdp.saveDipoleWRL("stockmayer.wrl",spc,sol);
  FormatPQR().save("confout.pqr", spc.p);
  spc.save("state");

  return test.numFailed();
}

int main() {
  cout << "Equilibration done ! Code: " << runSim("stockmayerEQ.input") << "!" << endl;
  int prod = runSim("stockmayer.input");
  cout << "Production done ! Code: " << prod << "!" << endl;
  return prod;
}