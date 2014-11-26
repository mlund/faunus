#include <faunus/faunus.h>
using namespace Faunus;

#ifdef CUBOID
typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
#else
typedef Geometry::Sphere Tgeometry;   // sphere with hard boundaries
#endif
typedef Potential::CoulombHS Tpairpot;// particle pair potential: primitive model

typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  
  InputMap mcp("gcmol.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);

  // Create Space and a Hamiltonian with three terms
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);

  string file = mcp.get<string>("polymer_file", "");

  //
  //    ADD CONFIGURATIONS FOR POOL INSERTS
  //
  p_vec conf;
  FormatAAM::load(file,conf);
  molecule.pushConfiguration("polymer", conf);   // p_vec is one molecule large
  molecule.pushConfiguration("polymer2",conf);

  // Markov moves and analysis
  Move::GCMolecular<Tspace> gc(mcp, pot, spc);
  Move::TranslateRotate<Tspace> mv(mcp,pot,spc);

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << molecule.info() << gc.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain
    while ( loop[1] ) {
      sys+=gc.move();
    }
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
  } // end of macro loop

  // save to disk

  spc.save("state");
  FormatPQR::save("confout.pqr", spc.p);

  sys.test(test);
  gc.test(test);

  // print information
  cout << loop.info() << sys.info() << gc.info() << spc.info() << test.info();

  //
  // clean allocated memory
  //
  spc.freeGroups();

  return test.numFailed();
}

/**
 * @page example_GCMolecular Example: Grand Canonical ensemble
 *
 * This is a example of Grand Canonical Monte Carlo simulation
 * with the following characteristics:
 * - Polymer insertion/deletion move
 * - Salt insertion/deletion move
 *
 * Run the example directly from the faunus directory:
 *
 *     $ make
 *     $ cd src/examples/
 *     $ ./gcmol.run
 *
*
* Input
* =====
* All listed files including the above C++ and python programs can be found in `src/examples/`
*
* gcmol.input    {#gcmol_input}
* -----------
* \include gcmol.input
*
*
* gcmol.aam    {#gcmol_aam}
* -----------
* configuration of rigid polymer
* \include gcmol.aam
*
*
* gcmol.json    {#gcmol_json}
* -----------
* State atom types
* \include gcmol.json
*
*
* gcmol.topo    {#gcmol_topo}
* -----------
* State molecular types and combinations for grand canonical move
*
* \include gcmol.topo
*
* * gcmol.cpp    {#gcmol_cpp}
* -----------
*
* \include gcmol.cpp
*/


