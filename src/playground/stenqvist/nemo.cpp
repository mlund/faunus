#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid> Tspace;
typedef DebyeHuckel TpairpotEl;
//typedef Coulomb TpairpotEl;
//typedef WeeksChandlerAndersen TpairpotSR;
typedef HardSphere TpairpotSR;
typedef CombinedPairPotential<TpairpotSR,TpairpotEl> Tpairpot;

/**
 * @brief Returns the direction of the cap.
 */
template<class Tspace, class Tgroup>
Point getDirectionOfCap(Tspace &spc, Tgroup& gi) {
  Point centerAtom(0,0,0);
  Point backAtom(0,0,0);
  int cnt1 = 0;
  int cnt2 = 0;
  for (auto i : *(*gi)) {
    Point t = spc.p[i] - (*gi)->cm;
    spc.geo.boundary(t);
    if(cnt1 == 0) {
      cnt1 = 1;
      centerAtom = t;
      continue;
    }
    if(cnt2 == 0) {
      cnt2 = 1;
      backAtom = t;
      break;
    }
  }
  Point temp = (centerAtom - backAtom);
  return temp/temp.norm();
}

template<class Tspace, class Tanalyze>
void analyzeDirection(Tspace &spc, Tanalyze &mucorrCC, Tanalyze &mucorrCS) {
	auto beg=spc.groupList().begin();
	auto end=spc.groupList().end();
	for (auto gi=beg; gi!=end; ++gi) {
	  for (auto gj=gi; ++gj!=end;) {
	   if ((*gi)->name.compare("cap") == 0 && (*gj)->name.compare("cap") == 0) {
	      Point dirA = getDirectionOfCap(spc,gi);
	      Point dirB = getDirectionOfCap(spc,gj);
	      mucorrCC(spc.geo.dist((*gi)->cm,(*gj)->cm)) += dirA.dot(dirB);
	    }
	   if ((*gi)->name.compare("cap") == 0 && (*gj)->name.compare("sphere") == 0) {
	      Point dirA = getDirectionOfCap(spc,gi);
	      Point dirB = -spc.geo.vdist((*gi)->cm,(*gj)->cm);
	      dirB = dirB/dirB.norm();
	      mucorrCS(spc.geo.dist((*gi)->cm,(*gj)->cm)) += dirA.dot(dirB);
	    }
	   if ((*gi)->name.compare("sphere") == 0 && (*gj)->name.compare("cap") == 0) {
	      Point dirA = spc.geo.vdist((*gj)->cm,(*gi)->cm);
	      dirA = dirA/dirA.norm();
	      Point dirB = getDirectionOfCap(spc,gj);
	      mucorrCS(spc.geo.dist((*gi)->cm,(*gj)->cm)) += dirA.dot(dirB);
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
  
  int traj_nbrs = mcp["traj_nbrs"];
  int micro_steps = 0;//mcp["system"]["micro"] | 10;
  int macro_steps = 0;//mcp["system"]["macro"] | 10;
  cout << "Microsteps: " << micro_steps << ", macrosteps: " << macro_steps << endl;
  
  double resolution = 0.2;

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);
  auto caps = spc.findMolecules("cap");
  auto spheres = spc.findMolecules("sphere");
  
  int sphereIndex = 0;
  auto beg=spc.groupList().begin();
  auto end=spc.groupList().end();
  for (auto gi=beg; gi!=end; ++gi) {
    if ((*gi)->name.compare("sphere") == 0 ) {
      for (auto i : *(*gi))
	sphereIndex = i;
      break;
    }
  }
  spc.p[sphereIndex] = Point(0,0,0);
  spc.trial = spc.p;
  cout << "Sphere: " << spc.p[sphereIndex].transpose() << endl;
  
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
          //xtc.setbox( spc.geo.len );
          xtc.save( "traj.xtc", spc.p );
        }
      }
    }
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); 
    cout << loop.timing();
  } 

  // save to disk
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);

  spc.save("state");
  rdfCC.save("rdfCC.dat");
  rdfCS.save("rdfCS.dat");
  rdfSS.save("rdfSS.dat");
  mucorrCC.save("mucorrCC.dat");
  mucorrCS.save("mucorrCS.dat");

  // print information
  cout << loop.info() + sys.info() + mv.info() + pot.info();
  
  cout << "Sphere: " << spc.p[sphereIndex].transpose() << endl;

  return 0;
}

