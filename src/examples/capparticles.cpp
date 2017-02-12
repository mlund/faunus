#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid, CapParticle> Tspace;
typedef HardSphereCap Tpairpot;

template<class Tspace, class Tanalyze>
void analyzeDirection(Tspace &spc, Tanalyze &mucorrCC, Tanalyze &mucorrCS) {
  
     for(unsigned int i = 0; i < spc.p.size()-1; i++) {
      for(unsigned int j = i+1; j < spc.p.size(); j++) {
       if(spc.p[i].is_sphere && !spc.p[j].is_sphere) {
	 Point dirA = spc.p[j].cap_center_point/spc.p[j].cap_center;
	 Point dirB = spc.geo.vdist(spc.p[i],spc.p[j]);
	 dirB = dirB/dirB.norm();
	 mucorrCS(spc.geo.dist(spc.p[i],spc.p[j])) += dirA.dot(dirB);
       }
       if(!spc.p[i].is_sphere && spc.p[j].is_sphere) {
	 Point dirA = spc.p[i].cap_center_point/spc.p[i].cap_center;
	 Point dirB = -spc.geo.vdist(spc.p[i],spc.p[j]);
	 dirB = dirB/dirB.norm();
	 mucorrCS(spc.geo.dist(spc.p[i],spc.p[j])) += dirA.dot(dirB);
       }
       if(!spc.p[i].is_sphere && !spc.p[j].is_sphere) {
	 Point dirA = spc.p[i].cap_center_point/spc.p[i].cap_center;
	 Point dirB = spc.p[j].cap_center_point/spc.p[j].cap_center;
	 mucorrCC(spc.geo.dist(spc.p[i],spc.p[j])) += dirA.dot(dirB);
       }
      }
     }
}

int main() {
  InputMap mcp("capparticles.json");   
  Tspace spc(mcp);
  spc.load("state");
  
  int seed = mcp["seed"];
  bool saveState = mcp["saveState"] | true;
  int traj_nbrs = mcp["traj_nbrs"] | 10;
  int steps = ( mcp["system"]["mcloop"]["micro"] | -10 )*(mcp["system"]["mcloop"]["macro"] | -1);
  double fraction = 1.0 - double(traj_nbrs)/double(steps);
  if(fraction > 1.0)
    fraction = 1.0;
  double resolution = mcp["resolution"] | 0.2;
  cout << "Seed: " << seed << ", saveState: " << saveState << ", traj_nbrs: " << traj_nbrs << ", steps: " << steps << ", fraction: " << fraction << ", resolution: " << resolution << endl;
  slump.seed(seed);
  
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp) + Energy::ExternalPressure<Tspace>(mcp);
  auto caps = spc.findMolecules("cap");
  auto spheres = spc.findMolecules("sphere");
  
  Move::Propagator<Tspace> mv( mcp, pot, spc );
  Analysis::RadialDistribution<> rdfCC(resolution);
  Analysis::RadialDistribution<> rdfCS(resolution);
  Analysis::RadialDistribution<> rdfSS(resolution);
  Table2D<double,Average<double> > mucorrCC, mucorrCS;
  mucorrCC.setResolution(resolution);
  mucorrCS.setResolution(resolution);
  FormatXTC xtc(2.0*traj_nbrs);
  EnergyDrift sys; 

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); 
  cout << atom.info() + spc.info() + pot.info();
  
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
      }
      if ( rnd > fraction ) {
        xtc.setbox( spc.geo.len );
        xtc.save( "traj.xtc", spc.p );
      }
    }
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); 
    cout << loop.timing();
  } 

  // save to disk
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);

  if(saveState)
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

