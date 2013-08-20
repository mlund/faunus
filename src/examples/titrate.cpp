#include <faunus/faunus.h>

using namespace Faunus;
typedef Space<Geometry::Sphere> Tspace;

int main() {
  cout << textio::splash();            // faunus splash info!
  ::atom.includefile("titrate.json");  // load atom properties
  InputMap mcp("titrate.input");
  EnergyDrift sys;                     // track system energy drift
  UnitTest test(mcp);

  // Space and hamiltonian
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Potential::DebyeHuckel>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  // Add molecule to middle of simulation container
  string file = mcp.get<string>("molecule", string());
  Tspace::ParticleVector v;
  FormatAAM::load(file,v);
  Geometry::cm2origo(spc.geo,v); // center molecule
  Group g = spc.insert(v);       // insert into space
  g.name="peptide"; // babtise
  spc.enroll(g);    // enroll group in space

  spc.load("state");

  Move::SwapMove<Tspace> tit(mcp,pot,spc,pot.second);
  Analysis::ChargeMultipole mp;

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  cout << atom.info() + spc.info() + pot.info() + tit.info()
    + textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);            // class for handling mc loops
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=tit.move();
      mp.sample(g,spc);
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );
    cout << loop.timing();

  } // end of macro loop

  cout << loop.info() + sys.info() + tit.info() + g.info() + mp.info();

  FormatPQR().save("confout.pqr", spc.p);
  spc.save("state");
}
