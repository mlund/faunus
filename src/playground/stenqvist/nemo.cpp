#include <faunus/faunus.h>
#include <faunus/ewald.h>

using namespace Faunus;
using namespace Faunus::Potential;

#define EWALD

typedef Space<Geometry::Cuboid> Tspace;
typedef LennardJonesLB_SF TpairpotLJ;
#ifdef EWALD
typedef TpairpotLJ Tpairpot;
typedef IonIonEwald TpairpotEl;
#else
typedef CombinedPairPotential<TpairpotEl,TpairpotLJ> Tpairpot;
typedef IonIonQ TpairpotEl; // IonIonQ IonIonSP3  MultipoleWolf<true,false,false,false,true>
#endif

int main() {
  InputMap mcp("nemo.json");   
  Tspace spc(mcp);
  int seed = mcp["seed"];
  cout << "Seed: " << seed << endl;
  spc.load("state"); 
  slump.seed(seed);

#ifdef EWALD
  auto pot = Energy::NonbondedEwald<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#else
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#endif
  auto potLJ = Energy::NonbondedVector<Tspace,TpairpotLJ>(mcp);
  auto potEl = Energy::NonbondedVector<Tspace,TpairpotEl>(mcp);
  
  auto waters = spc.findMolecules("water");
  Move::Propagator<Tspace> mv( mcp, pot, spc );
  vector<double> energy_T;
  vector<double> energy_LJ;
  //vector<double> energy_EP;
  Analysis::RadialDistribution<> rdfOO(0.02);
  Analysis::RadialDistribution<> rdfOH(0.02);
  Analysis::RadialDistribution<> rdfHH(0.02);
  Table2D<double,Average<double> > mucorr;
  mucorr.setResolution(0.02);
  Analysis::MultipoleAnalysis<Tspace> dian(spc,pot,mcp);
  FormatXTC xtc(1000);

  EnergyDrift sys; 
#ifdef EWALD
  pot.setSpace(spc);
#endif

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); 
  cout << atom.info() + spc.info() + pot.info();
  
  int cnt = 0;
  MCLoop loop(mcp);   
  while ( loop[0] ) {        
    while ( loop[1] ) {

      sys+=mv.move();

      double rnd = slump();
      if ( rnd > 0.8 ) {
        rdfOO.sample( spc, atom["OW"].id, atom["OW"].id );
        rdfOH.sample( spc, atom["OW"].id, atom["HW"].id );
        rdfHH.sample( spc, atom["HW"].id, atom["HW"].id );

        energy_T.push_back(Energy::systemEnergy(spc,pot,spc.p));
	energy_LJ.push_back(Energy::systemEnergy(spc,potLJ,spc.p));
	//energy_EP.push_back(Energy::systemEnergy(spc,potEl,spc.p));

	auto beg=spc.groupList().begin();
	auto end=spc.groupList().end();
	for (auto gi=beg; gi!=end; ++gi)
	  for (auto gj=gi; ++gj!=end;) {
	    Point mu_a = dipoleMoment(spc,*(*gi));
	    Point mu_b = dipoleMoment(spc,*(*gj));
	    double sca = mu_a.dot(mu_b);
	    double r = spc.geo.dist((*gi)->cm,(*gj)->cm); 
	    mucorr(r) += sca;
	  }

        if ( rnd > 0.99 ) {
          xtc.setbox( spc.geo.len );
          xtc.save( "traj.xtc", spc.p );
        }
        dian.sample(spc);
      }
    } 
    dian.save(std::to_string(cnt++));
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); 

    cout << loop.timing();
  } 

  // save to disk
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
  spc.save("state");
  rdfOO.save("rdfOO.dat");
  rdfOH.save("rdfOH.dat");
  rdfHH.save("rdfHH.dat");
  mucorr.save("mucorr.dat");
  dian.save();

  // print information
  cout << loop.info() + sys.info() + mv.info() + pot.info() + dian.info() + potEl.info();

  string file = "energy.dat";
  std::ofstream f(file.c_str());
  if (f)
    for (unsigned int i = 0; i < energy_T.size(); i++)
      f << std::left << std::setw(10) << energy_T.at(i)/spc.p.size()*3.0 << "     " << energy_LJ.at(i)/spc.p.size()*3.0  << endl;

  return 0;
}
