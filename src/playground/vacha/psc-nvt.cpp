#include <faunus/faunus.h>
using namespace Faunus;

typedef Space<Geometry::Cuboid,CigarParticle> Tspace;

using namespace Faunus::Potential;

typedef CombinedPairPotential<CosAttractMixed<>,WeeksChandlerAndersen> Tpair;
//the first interaction is a patchy, the second one is isotropic
typedef CigarSphereSplit<Tpair,Tpair,Tpair> Tpairpot;

int main() {
  cout << textio::splash();           // show faunus banner and credits
  InputMap mcp("psc-nvt.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  FormatMXYZ mxyz;                      // PQR structure file I/O
  EnergyDrift sys; // class for tracking system energy drifts

  // Energy functions and space
  Tspace spc(mcp);
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);

  // Markov moves and analysis
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::AtomicRotation<Tspace> rot(mcp, pot, spc);
  Analysis::RadialDistribution<> rdf(0.2);

  // Add cigars
  Group cigars;
  cigars.addParticles(spc, mcp);
    
    cout << atom[spc.p[0].id].name << endl;

  cigars.name="PSC1";
  for (auto i : cigars) {
    spc.p[i].dir.ranunit(slp_global);
    spc.p[i].patchdir.ranunit(slp_global);
    Geometry::cigar_initialize(spc.geo, spc.p[i]);
    spc.trial[i]=spc.p[i];
  }

  spc.load("state");                                     // load old config. from disk (if any)

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  std::ofstream m("snapshot");
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

      // movie
      if (slp_global()<0.001)
      {
        m << spc.p.size() << "\n"
          << "sweep " << loop.count() << "; box "
          << spc.geo.len.transpose() << "\n";
        for (auto &i : spc.p) {
          m << i.x() << " " << i.y() << " " << i.z() << "     "
            << i.dir.x() << " " << i.dir.y() << " " << i.dir.z() << "    "
            << i.patchdir.x() << " " << i.patchdir.y() << " " << i.patchdir.z() << " "
            << "\n";
        }

      }

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current

    cout << loop.timing();

  } // end of macro loop

  // save to disk
  rdf.save("rdf_p2p.dat");
  mxyz.save("confout.mxyz", spc.p, spc.geo.len, loop.count());
  spc.save("state");

  // print information
  cout << loop.info() << sys.info() << mv.info() << rot.info();
}
