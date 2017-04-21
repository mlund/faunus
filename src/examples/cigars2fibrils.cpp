#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid, CigarParticle> Tspace;
typedef CombinedPairPotential<CosAttractMixed<>, WeeksChandlerAndersen> Tpair;
//the first interaction is a patchy, the second one is isotropic
typedef CigarSphereSplit<Tpair, Tpair, Tpair> Tpairpot;

int main()
{
    cout << textio::splash();
    InputMap mcp("cigars2fibrils.json");
    FormatMXYZ mxyz;

    // Energy functions and space
    Tspace spc(mcp);
    Energy::NonbondedVector<Tspace, Tpairpot> pot(mcp);

    // Markov moves and analysis
    Move::Propagator<Tspace> mv(mcp, pot, spc);

    cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");

    MCLoop loop(mcp);
    mxyz.save("cigars2fibrils-mov", spc.p, spc.geo.len, loop.innerCount());

    while ( loop[0] )
    {  // Markov chain
        while ( loop[1] )
            mv.move();

        cout << loop.timing();

        mxyz.save("cigars2fibrils-mov", spc.p, spc.geo.len, loop.innerCount());
    } // end of macro loop

    // print information about simulationat the end
    cout << loop.info() + mv.info();
}

/**
  @page example_cigars2fibrils Example: Fibrils from patchy spherocylinders

  In this example we simulate the self-assembly of patchy spherocylinders (cigars)
  in the NVT ensemble.
  A spherocylinder is a cylinder with hemispherical caps at both ends, thus 
  providing a smooth rod-like particle with a variable aspect ratio.
  An attractive angular patch was added to the spherocylinder side,
  creating a patchy spherocylinder. This patch  can run either along the whole
  axis including the ends (patchtype 1) as in this example, or it can be
  limited to the cylindrical part (patchtype 2).
  The attractive interaction is determined from overlapping segments of
  the two spherocylinders within their patches.
  We use a Weeks-Chandler-Anderson potential for the repulsive part of the
  particle combined with a cosine squared potential for attractive
  interactions of patches.
  This is an implicit solvent model, where interactions are effective
  interactions and include solvent effects, hydrogen bonds, etc.

  In this short example we use an angular patch size of
  of 90 degrees. This results in fibril formation with three
  spherocylinders in cross-section (see figure below). Note that the
  timescale of this short example does not allow to reach equilibrium
  conditions.
  For more details see input files in `src/examples/cigars2fibrils.run`.

  The `cigars2fibrils.cpp` program can be used to simulate self-assembly of
  patchy spherocylinders with various aspect rations and patches at various
  temperatures.
  We have the following MC moves:
  - cigars translation
  - cigars rotations

  After the simulation one can run a visualisation script `cigars-mov2pdb.py` to
  convert the output trajectory to input files, which can be visualised by
  VMD (Visual Molecular Dynamics) program freely available at
  `http://www.ks.uiuc.edu/Research/vmd`. Note that it is not 
  straight forward to visualise the angular wedge on the spherocylinders
  in VMD. To guide the eye for the direction and size of the patch we display
  a additional spherocylinder that represents the patch.

  ![Snapshot of short fibrils formed in cigars2fibrils example with patchye spherocylinders. Blue color represents the body of spherocylinders while red spherocylinders are used for visualisation of the patch direction.](cigars2fibrils.png)

  cigars2fibrils.cpp
  ==================
  @includelineno examples/cigars2fibrils.cpp

*/
