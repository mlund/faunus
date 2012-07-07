#include <faunus/faunus.h>
#include <tclap/CmdLine.h>

using namespace Faunus;
using namespace TCLAP;

#ifdef CUBOID
typedef Geometry::Cuboid Tgeometry;
#else
typedef Geometry::Sphere Tgeometry;
#endif
typedef Potential::CombinedPairPotential<Potential::Coulomb, Potential::HardSphere> Tpairpot;

int main(int argc, char** argv) {
  string inputfile,istate,ostate;
  try {
    cout << textio::splash();
    CmdLine cmd("NPT Monte Carlo simulation of bead-chain polymers in explicit salt", ' ', "0.1");
    ValueArg<string> inputArg("i","inputfile","InputMap key/value file",true,"","inputfile");
    ValueArg<string> istateArg("c","instate","Name of input statefile",false,"state","instate");
    ValueArg<string> ostateArg("o","outstate","Name of output statefile",false,"state","outstate");
    cmd.add( inputArg );
    cmd.add( istateArg );
    cmd.add( ostateArg );
    cmd.parse( argc, argv );
    inputfile = inputArg.getValue();
    istate = istateArg.getValue();
    ostate = ostateArg.getValue();
  }
  catch (ArgException &e)  {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  InputMap mcp(inputfile);
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatTopology top;
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Analysis::PolymerShape shape;
  Analysis::RadialDistribution<float,int> rdf(0.2);

  // Add salt
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";

  // Add polymers
  vector<GroupMolecular> pol( mcp.get("polymer_N",0));
  string polyfile = mcp.get<string>("polymer_file", "");
  double req    = mcp.get<double>("polymer_eqdist", 0);
  double k      = mcp.get<double>("polymer_forceconst", 0);
  atom["MM"].dp = 10.;
  for (auto &g : pol) {                    // load polymers
    aam.load(polyfile);
    Geometry::FindSpace f;
    f.find(*spc.geo, spc.p, aam.p);        // find empty spot in particle vector
    g = spc.insert( aam.p );               // insert into space
    g.name="Polymer";
    spc.enroll(g);
    for (int i=g.front(); i<g.back(); i++)
      bonded->add(i, i+1, Potential::Harmonic(k,req)); // add bonds
  }
  Group allpol( pol.front().front(), pol.back().back() ); // make group w. all polymers

  spc.load(istate); // load old configuration from disk (if any)

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 4;
      switch (i) {
        case 0:
          mv.setGroup(salt);
          sys+=mv.move( salt.size() ); // translate salt
          break;
        case 1:
          mv.setGroup(allpol);
          sys+=mv.move( allpol.size() ); // translate monomers
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 2:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move(); // translate/rotate polymers
          }
          break;
        case 3:
          sys+=iso.move(); // isobaric volume move
          break;
      }
      for (auto i=pol.begin(); i!=pol.end()-1; i++)
        for (auto j=i+1; j!=pol.end(); j++)
          rdf( spc.geo->dist(i->cm,j->cm) )++;
      if ( slp_global.runtest(0.1) ) {
        //xtc.save("traj.xtc", spc);
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p)  ); // compare energy sum with current total

    cout << loop.timing();
  } // end of macro loop

  rdf.save("rdf_p2p.dat");
  pqr.save("confout.pqr", spc.p);
  top.save("mytopol.top", spc);
  spc.save(ostate);

  iso.test(test);
  gmv.test(test);
  mv.test(test);
  sys.test(test);

  cout << loop.info() << sys.info() << mv.info() << gmv.info() << iso.info() << shape.info() << test.info();
}
