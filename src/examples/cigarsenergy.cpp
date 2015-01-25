#include <faunus/faunus.h>
using namespace Faunus;
typedef Space<Geometry::Cuboid,CigarParticle> Tspace;
using namespace Faunus::Potential;
typedef CombinedPairPotential<CosAttractMixed<>,WeeksChandlerAndersen> Tpair;
//the first interaction is a patchy, the second one is isotropic
typedef CigarSphereSplit<Tpair,Tpair,Tpair> Tpairpot;

int main() {
  InputMap mcp("cigarsenergy.json");   // open user input file

  Tspace spc(mcp);
  UnitTest test(mcp); 
  Energy::NonbondedVector<Tspace,Tpairpot> pot(mcp);

  FormatMXYZ::load("cigarsenergy.xyz", spc.p, spc.geo.len);

  test("energy", Energy::systemEnergy(spc,pot,spc.p), 1e-3 );

  std::cout << pot.info() << test.info();

  return test.numFailed();
}

