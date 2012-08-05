#include <faunus/faunus.h>
#include <tclap/CmdLine.h>

using namespace Faunus;
using namespace TCLAP;

typedef Geometry::Cuboid Tgeometry;
typedef Potential::CombinedPairPotential<Potential::SquareWellHydrophobic, Potential::LennardJones> SRpot;
typedef Potential::CombinedPairPotential<Potential::DebyeHuckel, SRpot> Tpairpot;
//typedef Potential::CoulombSR<Tgeometry, Potential::DebyeHuckel, Potential::HardSphere> Tpairpot;

int main(int argc, char** argv) {
  string inputfile,istate,ostate;
  try {
    cout << textio::splash();
    CmdLine cmd("NPT Monte Carlo simulation of rigid bodies in continuum", ' ', "0.1");
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
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
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
  auto constrain = pot.create( Energy::MassCenterConstrain(pot.getGeometry()) );
  Space spc( pot.getGeometry() );

  //Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::SwapMove tit(mcp,pot,spc);
  //Move::SwapMoveMSR tit(mcp,pot,spc);
  Analysis::RadialDistribution<float,int> rdf(0.25);
  Analysis::ChargeMultipole poleTotal;
  Analysis::ChargeMultipole poleIonic;
  Analysis::ChargeMultipole poleProton;
  string X=mcp.get<string>("anion", "");
  poleProton.exclusionlist={X+"ALA",X+"ILE",X+"LEU", X+"MET", X+"PHE", X+"PRO", X+"TRP", X+"VAL", X+"Bbs", X+"Bbw"};
  poleIonic.exclusionlist={"ASP", "CTR", "GLU", "HHIS", "HNTR", "TYR", "HLYS", "CYS", "HARG"};
  
  //Rename buried titratable groups determined by their solvent accessible surface area 
  io fio;
  vector<string> v_sasa;

  int N1 = mcp.get("molecule_N1",0);                                                                          
  int N2 = mcp.get("molecule_N2",0);
  double threshold = mcp.get("sasa_threshold",25);

  string aamfile = mcp.get<string>("back_side1", "");
  string file_sasa = mcp.get<string>("sasa1", "");
  
  std::ostringstream of;
  of << threshold << "gammacryst-buried.aam";
  string outfile = of.str();
  
  std::istringstream stm;
  aam.load(aamfile);
  fio.readfile(file_sasa, v_sasa);
  assert (v_sasa.size()==aam.p.size() &&  "Sasa must be defined for each residue in aam file");
  int cnt_buried=0;

  for (unsigned int i=0; i < v_sasa.size(); i++){
    double sa = atof (v_sasa[i].c_str() );
    if (i < v_sasa.size()/2) {                         // Make sure that back bone particles come first in the aam file
      if (sa < threshold) {
        aam.p[i].id=atom["BBb"].id;                    // BBb = Buried Backbone 
        cnt_buried+=1;
      }
      else {
        if ((aam.p[i].id==atom["GLY"].id) || 
            (aam.p[i].id==atom["PRO"].id)){            // Anions show strong binding to PRO and GLY
          aam.p[i].id=atom["Bbs"].id;                  // Bbs = Backbone with strong binding
        }
        else { 
          aam.p[i].id=atom["Bbw"].id;                  // Bbw = Backbone with weak binding
        }
      }
    }
    else {
      if (sa < threshold){
        cnt_buried+=1;
        aam.p[i].id=atom["BSc"].id;                  // BSc = Buried Sidechain
      }
    }
  }
  aam.save(outfile, aam.p);

  // Add molecules
  string file;
  vector<GroupMolecular> pol(N1+N2);
  for (int i=0; i<N1+N2; i++) {
    GroupMolecular g;
    if (i>=N1)
      file = mcp.get<string>("molecule_file2", "");
    else
      file = mcp.get<string>("molecule_file1", "");
    aam.load(file);
    Geometry::FindSpace f;
    f.find(*spc.geo, spc.p, aam.p);        // find empty spot in particle vector
    pol[i] = spc.insert( aam.p );          // insert into space
    pol[i].name=file;
    spc.enroll( pol[i] );
  }
  if (pol.size()>1)
    constrain->addPair(pol[0], pol[1], 0, 90);

  Group allpol( pol.front().front(), pol.back().back() );
  
  for(unsigned int i=0; i < spc.p.size(); i++){
    spc.p[i].charge=atom[spc.p[i].id].charge;
    spc.trial[i].charge=spc.p[i].charge;
  }
  
  tit.findSites(spc.p);  // search for titratable sites
  spc.load(istate);

  double utot=pot.external() + pot.g_internal(spc.p, allpol);
  for (auto &g : pol)
    utot += pot.g_external(spc.p, g);
  sys.init( utot );

  pqr.save("initial.pqr",spc.p);
  cout << atom.info() << spc.info() << pot.info() << tit.info()
    << textio::header("MC Simulation Begins!");
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 2;
      switch (i) {
        case 0:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++)
              rdf( spc.geo->dist(i->cm,j->cm) )++;
          break;
        case 10:
          //sys+=iso.move();
          break;
        case 1:
          sys+=tit.move();
          poleProton.sample(pol, spc);
          poleTotal.sample(pol,spc);
          poleIonic.sample(pol,spc);
          break;
      }
      if (mcp.get<bool>("xtc_write","no"))
        if ( slp_global.runtest(0.0001) ) {
          //xtc.setbox( nonbonded->geo->len );
          xtc.save("traj.xtc", spc);
        }
    } // end of micro loop

    double utot=pot.external() + pot.g_internal(spc.p, allpol);
    for (auto &g : pol)
      utot += pot.g_external(spc.p, g);
    sys.checkDrift( utot );
    cout << loop.timing();

  } // end of macro loop

  /*
  iso.test(test);
  gmv.test(test);
  sys.test(test);
  */

  cout << loop.info() << sys.info() << gmv.info() << tit.info()<< endl 
       << textio::header("Total Charge Analysis") << poleTotal.info() 
       << textio::header("Charge Analysis w/only protons") << poleProton.info()
       << textio::header("Charge Analysis w/only anions") << poleIonic.info() << endl
       << "  Number of buried residues      "<< cnt_buried << endl;

  rdf.save("rdf_p2p.dat");
  pqr.save("confout.pqr", spc.p);
  top.save("mytopol.top", spc);
  spc.save(ostate);
}
