#include <faunus/faunus.h>
using namespace Faunus;
typedef Space<Geometry::Cuboid,CigarParticle> Tspace;
using namespace Faunus::Potential;
typedef CombinedPairPotential<CosAttractMixed<>,WeeksChandlerAndersen> Tpair;
//the first interaction is a patchy, the second one is isotropic
typedef CigarSphereSplit<Tpair,Tpair,Tpair> Tpairpot;

int main() {
  cout << textio::splash();                 // show faunus banner and credits
  //InputMap mcp("cigars2fibrils.input");     // open user input file
  InputMap mcp("psc-nvt.input");
  MCLoop loop(mcp);                         // class for handling mc loops
  FormatMXYZ mxyz;                          // MXYZ structure file I/O
  EnergyDrift sys;                         // class for tracking system energy drifts

    
  // Energy functions and space
  Tspace spc(mcp);
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);

  // Markov moves and analysis
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::AtomicRotation<Tspace> rot(mcp, pot, spc);

  cout <<"before adding cigars";

  // Add cigars
  Group cigars;
  cigars.addParticles(spc, mcp);
  cigars.name="PSC";
  for (auto i : cigars) {
    spc.p[i].dir.ranunit(slp_global);
    spc.p[i].patchdir.ranunit(slp_global);
    Geometry::cigar_initialize(spc.geo, spc.p[i]);
    spc.trial[i]=spc.p[i];
  }
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  mxyz.save("cigars2fibrils-mov", spc.p, spc.geo.len, loop.count());
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=slp_global.rand() % 2;
      switch (i) {
        case 0:
          mv.setGroup(cigars);
          sys+=mv.move( cigars.size() );  // translate cigars
          break;
        case 1:
          rot.setGroup(cigars);
          sys+=rot.move( cigars.size() ); // translate cigars
          break;
      }
    } // end of micro loop
  sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current

  cout << loop.timing();

  mxyz.save("cigars2fibrils-mov", spc.p, spc.geo.len, loop.count());
  } // end of macro loop

  // print information
  cout << loop.info() << sys.info() << mv.info() << rot.info();
}
