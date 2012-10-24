/*
 * We simulate a single polyelectrolyte in a rectangular slit where one wall
 * is charged. The following interactions are included,
 * 1) monomer-surface: Gouy-Chapman electrostatics
 * 2) monomer-monomer: stiff bonds - i.e. ideal chain
 *
 * This is of course an unphysical system to match exact polymer DFT theory.
 */

#include <faunus/faunus.h>

using namespace Faunus;

int main() {
  cout << textio::splash();

  InputMap mcp("slitpolymer.input");
  MCLoop loop(mcp);      // class for handling mc loops
  FormatPQR pqr;         // PQR structure file I/O
  FormatAAM aam;         // AAM structure file I/O
  EnergyDrift sys;       // class for tracking system energy drifts
  UnitTest test(mcp);

  // Set up energy field
  Geometry::Cuboidslit geo(mcp);
  Energy::GouyChapman pot(mcp);
  pot.setGeometry(geo);
  pot.setPosition( geo.len_half.z ); // Surface in xy plane at +z direction
  Space spc( pot.getGeometry() );

  // Add polymer
  string polyfile = mcp.get<string>("polymer_file", "");
  aam.load(polyfile);
  Geometry::FindSpace f;
  f.find(*spc.geo, spc.p, aam.particles()); // find empty spot
  GroupMolecular pol = spc.insert( aam.particles() );  // insert into space
  pol.name="polymer";
  spc.enroll(pol);

  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp,pot,spc);
  Move::CrankShaft crank(mcp,pot,spc);
  Move::Pivot pivot(mcp,pot,spc);
  Move::Reptation rep(mcp,pot,spc);

  Analysis::PolymerShape shape;
  Analysis::LineDistribution<float,int> surfmapall;

  spc.load("state");
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  cout << spc.info() << pot.info() << pol.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=rand() % 4;
      switch (i) {
        case 0: // translate and rotate polymer
          gmv.setGroup(pol);
          sys+=gmv.move(); 
          break;
        case 1: // pivot
          pivot.setGroup(pol);
          sys+=pivot.move();
          break;
        case 2: // crankshaft
          crank.setGroup(pol);
          sys+=crank.move();
          break;
        case 3: // reptation
          rep.setGroup(pol);
          sys+=rep.move();
          break;
      }

      pol.setMassCenter(spc);
      shape.sample(pol,spc);
      for (int i=pol.front(); i<=pol.back(); i++)
        surfmapall( pot.dist2surf( spc.p.at(i) ) )++;

    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );
    cout << loop.timing();

  } // end of macro loop

  pqr.save("confout.pqr", spc.p);
  spc.save("state");
  surfmapall.save("surfall.dat");

  gmv.test(test);
  mv.test(test);
  sys.test(test);
  pivot.test(test);
  crank.test(test);
  rep.test(test);
  shape.test(test);

  cout << loop.info() << sys.info() << gmv.info() << mv.info() << crank.info()
    << pivot.info() << rep.info() << shape.info() << spc.info() << test.info();

  return test.numFailed();
}
