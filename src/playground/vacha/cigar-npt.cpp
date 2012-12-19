#include <faunus/faunus.h>
using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries

using namespace Faunus::Potential;

typedef CombinedPairPotential<CosAttract,WeeksChandlerAndersen> Tpair;
typedef CigarSphereSplit<Tpair,Tpair,Tpair> Tpairpot;

int main() {
  cout << textio::splash();           // show faunus banner and credits
  InputMap mcp("cigar-npt.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  FormatPQR pqr;                      // PQR structure file I/O
  EnergyDrift sys;                    // class for tracking system energy drifts

  // Energy functions and space
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::NonbondedVector<Tpairpot,Tgeometry>(mcp) );
  Space spc( pot.getGeometry() );

  // Markov moves and analysis
  Move::Isobaric iso(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::AtomicRotation rot(mcp, pot, spc);
  Analysis::RadialDistribution<float,int> rdf(0.2);

  // Add cigars
  GroupAtomic cigars(spc, mcp);
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
      int k,i=slp_global.rand() % 3;
      switch (i) {
        case 0:
          mv.setGroup(cigars);
          sys+=mv.move( cigars.size() );  // translate cigars
          break;
        case 1:
          rot.setGroup(cigars);
          sys+=rot.move( cigars.size() );  // translate cigars
          break;
        case 2:
          sys+=iso.move();              // isobaric volume move
          break;
      }

      // movie
      if (slp_global()<0.001)
      {
        m << spc.p.size() << "\n"
          << "sweep " << loop.count() << "; box "
          << nonbonded->geometry.len.x() << " "
          << nonbonded->geometry.len.y() << " "
          << nonbonded->geometry.len.z() << "\n";
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
