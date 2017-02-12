#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,CapParticle> Tspace;
//typedef DebyeHuckel TpairpotEl;
typedef Coulomb TpairpotEl;
//typedef CombinedPairPotential<TpairpotEl,HardSphereCap> Tpairpot;
typedef HardSphereCap Tpairpot;

template<class Tspace, class Tanalyze>
void analyzeDirection(Tspace &spc, Tanalyze &mucorrCC, Tanalyze &mucorrCS) {
  for(unsigned int i = 0; i < spc.p.size()-1; i++) {
    for(unsigned int j = i+1; j < spc.p.size(); j++) {
      if(!spc.p[i].is_sphere && !spc.p[j].is_sphere) {
	Point dirA = spc.p[i].cap_center_point/spc.p[i].cap_center;
	Point dirB = spc.p[j].cap_center_point/spc.p[j].cap_center;
	mucorrCC(spc.geo.dist(spc.p[i],spc.p[j])) += dirA.dot(dirB);
      }
      if(spc.p[i].is_sphere && !spc.p[j].is_sphere) {
	double r = spc.geo.dist(spc.p[i],spc.p[j]);
	Point dirA = spc.geo.vdist(spc.p[i],spc.p[j])/r;
	Point dirB = spc.p[j].cap_center_point/spc.p[j].cap_center;
	mucorrCS(r) += dirA.dot(dirB);
      }
      if(!spc.p[i].is_sphere && spc.p[j].is_sphere) {
	double r = spc.geo.dist(spc.p[i],spc.p[j]);
	Point dirA = spc.p[i].cap_center_point/spc.p[i].cap_center;
	Point dirB = -spc.geo.vdist(spc.p[i],spc.p[j])/r;
	mucorrCS(r) += dirA.dot(dirB);
      }
    }
  }
}

int main() {
  InputMap mcp("cap.json");   
  Tspace spc(mcp);
  int seed = mcp["seed"];
  cout << "Seed: " << seed << endl;
  spc.load("state"); 
  slump.seed(seed);
  double resolution = 0.2;

  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);
  auto caps = spc.findMolecules("cap");
  auto spheres = spc.findMolecules("sphere");
  Move::Propagator<Tspace> mv( mcp, pot, spc );
  Analysis::RadialDistribution<> rdfCC(resolution);
  Analysis::RadialDistribution<> rdfCS(resolution);
  Analysis::RadialDistribution<> rdfSS(resolution);
  Table2D<double,Average<double> > mucorrCC, mucorrCS;
  mucorrCC.setResolution(resolution);
  mucorrCS.setResolution(resolution);
  FormatXTC xtc(1000);
  EnergyDrift sys; 

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); 
  cout << atom.info() + spc.info() + pot.info();
  
  analyzeDirection(spc,mucorrCC,mucorrCS);
  
  MCLoop loop(mcp);   
  while ( loop[0] ) {        
    while ( loop[1] ) {
      sys+=mv.move();

      double rnd = slump();
      if ( rnd > 0.8 ) {
        rdfCC.sample( spc, atom["C"].id, atom["C"].id );
        rdfCS.sample( spc, atom["C"].id, atom["S"].id );
        rdfSS.sample( spc, atom["S"].id, atom["S"].id );

	analyzeDirection(spc,mucorrCC,mucorrCS);

        if ( rnd > 0.99 ) {
          xtc.setbox( spc.geo.len );
          xtc.save( "traj.xtc", spc.p );
        }
      }
    }
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); 
    cout << loop.timing();
  } 

  // save to disk
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);

  if(mcp["saveState"] == 1)
    spc.save("state");
  rdfCC.save("rdfCC.dat");
  rdfCS.save("rdfCS.dat");
  rdfSS.save("rdfSS.dat");
  mucorrCC.save("mucorrCC.dat");
  mucorrCS.save("mucorrCS.dat");

  // print information
  cout << loop.info() + sys.info() + mv.info() + pot.info();

  return 0;
}
