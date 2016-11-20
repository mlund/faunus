#include <faunus/faunus.h>

using namespace Faunus;

/**
 * @brief Setup interactions for coarse grained membrane
 *
 * More information: doi:10/chqzjk
 */
template<class Tbonded, class Tpotmap, class Tlipid, class Tinput>
void MakeDesernoMembrane(const Tlipid &lipid, Tbonded &b, Tpotmap &m, Tinput &in) {

  using namespace Potential;

  // non-bonded interactions
  auto hid=atom["HD"].id;
  auto tid=atom["TL"].id;
  double sigma   = in["system"]["lipid_sigma"] | 10.0;  // angstrom
  double epsilon = in["system"]["lipid_epsilon"] | 1.0; // kT

  auto js = in["energy"]["nonbonded"];
  WeeksChandlerAndersen WCA(js);
  WCA.customSigma(hid, hid, 0.95*sigma);            
  WCA.customSigma(hid, tid, 0.95*sigma);            
  WCA.customSigma(tid, tid, sigma);                 
  WCA.customEpsilon(hid, hid, epsilon);             
  WCA.customEpsilon(hid, tid, epsilon);             
  WCA.customEpsilon(tid, tid, epsilon);   

  // Add to main potential
  m.add(hid, hid, WCA + DebyeHuckel(js) );
  m.add(tid, tid, WCA + CosAttract(js) );
  m.add(hid, tid, WCA + ChargeNonpolar(js) );

  // bonded interactions
  double headtail_k=0.5*10*epsilon/(sigma*sigma);
  double headtail_req=4*sigma;
  double fene_k=30*epsilon/(sigma*sigma);
  double fene_rmax=1.5*sigma;

  assert(lipid.size() % 3 == 0);

  for (auto i=lipid.front(); i<lipid.back(); i+=3) {
    b.add(i,  i+1, FENE(fene_k,fene_rmax) );
    b.add(i+1,i+2, FENE(fene_k,fene_rmax) );
    b.add(i,  i+2, Harmonic(headtail_k,headtail_req) );
  }
}

typedef Space<Geometry::Cuboid> Tspace;
typedef Potential::PotentialMap<Potential::DebyeHuckel> Tpairpot;

int main() {

  cout << textio::splash();      // show faunus banner and credits
  InputMap mcp("membrane.json"); //read input file

  EnergyDrift sys;               // class for tracking system energy drifts

  // Energy functions and space
  auto pot = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::Bonded<Tspace>()
    + Energy::ExternalPressure<Tspace>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  auto nonbonded = std::get<0>( pot.tuple() );
  auto bonded    = std::get<1>( pot.tuple() );

  nonbonded->noPairPotentialCutoff=true;
  Tspace spc(mcp);

  auto lipids = spc.findMolecules("lipid");

  Group allLipids( lipids.front()->front(), lipids.back()->back() );
  allLipids.setMolSize(3);
  MakeDesernoMembrane(allLipids, *bonded, nonbonded->pairpot, mcp);

  // Place all lipids in xy plane (z=0);
  for ( auto g : lipids ) {
    double dz = spc.p[ g->back() ].z();
    for (auto i : *g) {
      spc.p[i].z() -= dz;
      spc.geo.boundary( spc.p[i] );
    }
    g->setMassCenter( spc );
  }
  spc.trial=spc.p;   // sync. particle trial vector
  spc.load("state"); // load old config. from disk (if any)

  // Markov moves and analysis
  Move::Propagator<Tspace> mv(mcp, pot, spc);
  Analysis::BilayerStructure lipidstruct;
  Analysis::VirialPressure<Tspace> virial(mcp, pot, spc);

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store total energy

  cout << atom.info() + spc.info() + pot.info();

  MCLoop loop(mcp);                      // class for handling mc loops
  while ( loop[0] ) {                    // Markov chain 
    while ( loop[1] ) {
      sys += mv.move();
      double ran = slump();
      if ( ran > 0.90 ) {
        virial.sample();
        lipidstruct.sample(spc.geo, spc.p, allLipids);
      }

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?

    // save to disk
    FormatPQR::save("confout.pqr", spc.p);
    spc.save("state");

    cout << loop.timing();
  } // end of macro loop

  // perform unit tests
  UnitTest test(mcp);
  mv.test(test);
  sys.test(test);
  lipidstruct.test(test);

  cout << loop.info() + mv.info()
    + lipidstruct.info() + sys.info() + virial.info() + test.info();

  return test.numFailed();
}

/** @page example_membrane Example: Membrane Bilayer

  This will simulate a 3-bead coarse grained membrane according to
  [Cooke and Deserno](http://dx.doi.org/10/chqzjk). Each bead interacts with a
  Weeks-Chandler-Andersen potential, while tail-tail interactions
  have an additional long range attractive potential. There is preliminary
  support for charged head groups (effect on elastic properties is unknown).

  The following moves are included:
  - Lipid rotation, translation and pivot
  - Monomer translation
  - Iso-tension move

  Run this example from the `examples` directory:

  ~~~~~~~~~~~~~~~~~~~
  $ make
  $ cd src/examples
  $ ./membrane.run
  ~~~~~~~~~~~~~~~~~~~

  ![Bilayer formed by 3-bead CG lipid model](membrane-3bead.jpg)

  membrane.cpp
  ============

  @includelineno examples/membrane.cpp

*/

