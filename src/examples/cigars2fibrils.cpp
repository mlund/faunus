#include <faunus/faunus.h>
using namespace Faunus;
typedef Space<Geometry::Cuboid,CigarParticle> Tspace;
using namespace Faunus::Potential;
typedef CombinedPairPotential<CosAttractMixed<>,WeeksChandlerAndersen> Tpair;
//the first interaction is a patchy, the second one is isotropic
typedef CigarSphereSplit<Tpair,Tpair,Tpair> Tpairpot;

int main() {
  cout << textio::splash();                 // show faunus banner and credits
  InputMap mcp("cigars2fibrils.input");     // open user input file
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

  // print information about simulationat the end
  cout << loop.info() << sys.info() << mv.info() << rot.info();
}

/**
 @page example_cigars2fibrils Example: Fibrils from patchy spherocylinders
 
 In this example we simulate self-assembly of patchy spherocylinders (cigars)
 in the NVT ensemble.
 Spherocylinder is a cylinder with hemispherical caps at both ends, thus 
 providing a smooth rod-like particle with a variable aspect ratio.
 An attractive angular was added to the side of a spherocylinder,
 creating a patchy spherocylinder. This patch  can either run along the whole
 axis including the ends (patchtype 1) as in this example, or it can be
 limited to the cylindrical part (patchtype 2).
 The attractive interaction is then determined from overlapping segments of
 the two spherocylinders withn their patches.
 We use a Weeks-Chandler-Anderson potential for the repulsive part of the
 particle combined  with a cosine squared potential for attractive
 interactions of patches.
 This is an implicit solvent model, where interactions are effective
 interactions and including solvent effects, hydrogen bonds, etc.
 
 In this short example we use a angular size of a patch 
 of 90 degrees. This results in self-assembly of fibrils with three
 spherocylinders in cross-section (see figure below). Note that the
 timescale of this short example does not  allow to reach equilibrium
 conditions.
 For more details see input files in `src/examples/cigars2fibrils.run`.
 
 The `cigars2fibrils.cpp` program can be used to simulate self-assembly of
 patchy spherocylinders with various aspect rations and patches at various
 temperatures.
 We have the following MC moves:
 - cigars translation
 - cigars rotations
 
 After simulation one can run a visualisation script `cigars-mov2pdb.py` to
 convert the output trajectory to input files, which can be visualised by
 VMD (Visual Molecular Dynamics) program freely available at
 `http://www.ks.uiuc.edu/Research/vmd/`. Note that it is not 
 straight forward to visualise the angular wedge on the spherocylinders
 in VMD. To guied an eye for direction and size of patch we display
 additional spherocylinder that represents the patch.
 
![Snapshot of short fibrils formed in cigars2fibrils example with patchye spherocylinders. Blue color represents the body of spherocylinders while red spherocylinders is used for visualisation of patch direction.](cigars2fibrils.png)

 cigars2fibrils.cpp
 ========
 \includelineno examples/cigars2fibrils.cpp
 
 */
