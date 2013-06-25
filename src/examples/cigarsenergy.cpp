#include <faunus/faunus.h>
using namespace Faunus;
typedef Space<Geometry::Cuboid,CigarParticle> Tspace;
using namespace Faunus::Potential;
typedef CombinedPairPotential<CosAttractMixed<>,WeeksChandlerAndersen> Tpair;
//the first interaction is a patchy, the second one is isotropic
typedef CigarSphereSplit<Tpair,Tpair,Tpair> Tpairpot;

int main() {
//  cout << textio::splash();           // show faunus banner and credits
  InputMap mcp("cigarsenergy.input");   // open user input file
  FormatMXYZ mxyz;                      // xyz structure file I/O
  UnitTest test(mcp); 

  // Energy functions and space
  Tspace spc(mcp);
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);

  // Load and add cigars to Space
  Group cigars;
  cigars.addParticles(spc, mcp);
  cigars.name="PSC";
  spc.enroll(cigars);

  cout << atom.info() << cigars.info() << spc.info() << pot.info();

  mxyz.load("cigarsenergy.xyz", spc.p, spc.geo.len);
  for (auto i : cigars) {
    Geometry::cigar_initialize(spc.geo, spc.p[i]);
    spc.trial[i]=spc.p[i];
  }

  test("energy", Energy::systemEnergy(spc,pot,spc.p), 1e-3 );

  cout << test.info();

  return test.numFailed();
}

