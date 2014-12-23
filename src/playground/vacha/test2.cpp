#include <faunus/faunus.h>
using namespace Faunus;

typedef Space<Geometry::Cuboid,CigarParticle> Tspace;

using namespace Faunus::Potential;

typedef CombinedPairPotential<CosAttractMixed<>,WeeksChandlerAndersen> Tpair;
//the first interaction is a patchy, the second one is isotropic
typedef CigarSphereSplit<Tpair,Tpair,Tpair> Tpairpot;

int main() {
//  cout << textio::splash();           // show faunus banner and credits
  InputMap mcp("test2.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  FormatMXYZ mxyz;                      // PQR structure file I/O
  EnergyDrift sys; // class for tracking system energy drifts

  // Energy functions and space
  Tspace spc(mcp);
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);

  // Markov moves and analysis
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::AtomicRotation<Tspace> rot(mcp, pot, spc);

 // Load and add psc to Space
  Group cigars;
  cigars.addParticles(spc, mcp);
  mxyz.load("test2in.xyz", spc.p, spc.geo.len);
  //Tspace::ParticleVector v;                     // temporary, empty particle vector
  //mxyz.load("test2in.xyz",v ,spc.geo.len);      // loading data from MXYZ file
  //Geometry::FindSpace().find(spc.geo,spc.p,v);  // find empty spot in particle vector
  //Group cigars = spc.insert(v);                  // Insert into Space and return matching group
  cigars.name="PSC";                           // type of psc
  //spc.enroll(cigars);                            // Enroll psc in Space
  
  for (auto i : cigars) {
    Geometry::cigar_initialize(spc.geo, spc.p[i]);
    spc.trial[i]=spc.p[i];
  }

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

        for (auto &i : spc.p) {
          cout << i.x() << " " << i.y() << " " << i.z() << "     "
            << i.dir.x() << " " << i.dir.y() << " " << i.dir.z() << "    "
            << i.patchdir.x() << " " << i.patchdir.y() << " " << i.patchdir.z() << " "
            << "\n";
        }

  cout << sys.info();

//  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");
/*
  std::ofstream m("snapshot");
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=slump.rand() % 2;
        
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
      if (slump()<0.001)
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
  mxyz.save("confout.mxyz", spc.p, spc.geo.len,loop.count());
  spc.save("state");

   print information
  cout << loop.info() << sys.info() << mv.info() << rot.info();
*/
}
