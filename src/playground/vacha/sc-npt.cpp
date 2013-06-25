#include <faunus/faunus.h>
using namespace Faunus;

typedef Space<Geometry::Cuboid,CigarParticle> Tspace;

using namespace Faunus::Potential;

typedef CombinedPairPotential<CosAttract,WeeksChandlerAndersen> Tpair;
typedef CigarSphereSplit<Tpair,Tpair,Tpair> Tpairpot;

int main() {
  cout << textio::splash();           // show faunus banner and credits
  InputMap mcp("sc-npt.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  FormatPQR pqr;                      // PQR structure file I/O
  EnergyDrift sys; // class for tracking system energy drifts

  // Energy functions and space
  Tspace spc(mcp);
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);

  // Markov moves and analysis
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::AtomicRotation<Tspace> rot(mcp, pot, spc);
  Analysis::RadialDistribution<> rdf(0.2);

  // Add cigars
  Group cigars;
  cigars.addParticles(spc, mcp);
  cigars.name="cigars";
  for (auto i : cigars) {
    spc.p[i].halfl = 2.5;
    spc.p[i].dir.ranunit(slp_global);
    spc.trial[i]=spc.p[i];
  }

  spc.load("state");                                     // load old config. from disk (if any)

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  std::ofstream m("snapshot");
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=slp_global.rand() % 3;
      switch (i) {
        case 0:
          mv.setGroup(cigars);
          sys+=mv.move( cigars.size() );  // translate cigars
          break;
        case 1:
          rot.setGroup(cigars);
          sys+=rot.move( cigars.size() ); // translate cigars
          break;
        case 2:
          sys+=iso.move();                // isobaric volume move
          break;
      }

      // movie
      if (slp_global()<0.001)
      {
        m << spc.p.size() << "\n"
          << "sweep " << loop.count() << "; box "
          << spc.geo.len.x() << " "
          << spc.geo.len.y() << " "
          << spc.geo.len.z() << "\n";
        for (auto &i : spc.p) {
          m << i.x() << " " << i.y() << " " << i.z() << " "
            << i.dir.x() << " " << i.dir.y() << " " << i.dir.z() << " "
            << " 0 0 0 0\n";
        }

      }

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current

    cout << loop.timing();

  } // end of macro loop

  // save to disk
  rdf.save("rdf_p2p.dat");
  pqr.save("confout.pqr", spc.p);
  spc.save("state");

  // print information
  cout << loop.info() << sys.info() << mv.info() << rot.info() << iso.info();
}
