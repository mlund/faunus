#include <faunus/faunus.h>
#include <tclap/CmdLine.h>

using namespace Faunus;
using namespace std;

typedef Geometry::Cylinder Tgeometry;
typedef Potential::CoulombLJTS Tpairpot;

int main(int argc, char** argv) {
  Faunus::MPI::MPIController mpi;
  mpi.cout << textio::splash();
  InputMap mcp(textio::prefix+"temper.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  //auto constrain = pot.create( Energy::MassCenterConstrain(pot.getGeometry()) );
  Space spc( pot.getGeometry() );

  double pxtc = mcp.get("pxtc",0.99);

  // Add molecular species
  int cnt=0;
  int N1=mcp.get("polymer1_N",0);
  int N2=mcp.get("polymer2_N",0);
  double ppos[3];
  ppos[1] = mcp.get<double>("polymer1_pos",0);
  ppos[2] = mcp.get<double>("polymer2_pos",0);
  vector<GroupMolecular> pol( N1+N2);
  for (auto &g : pol) {
    cnt++;
    string polyfilekey = (cnt>N1) ? "polymer2_file" : "polymer1_file";
    aam.load( mcp.get<string>(polyfilekey, "") );
    Geometry::FindSpace f;
    //f.dir.x()=0; // put mass center
    //f.dir.y()=0; //   at [x,y,z] = [0,0,random]
    //f.dir.z()=1;
    //if (f.find(*spc.geo, spc.p, aam.p )) {
      Point v;
      v.x()=0.0;
      v.y()=0.0;
      v.z()=ppos[cnt];
      translate(*spc.geo, aam.p, -massCenter(*spc.geo, aam.p)-v);
      g = spc.insert( aam.p );
      g.name=mcp.get<string>(polyfilekey, "");
      spc.enroll(g);
    //} else
    //  return 1;
  }

  //constrain->addPair(pol[0], pol[1], 20, 150);

  // Add salt
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  spc.enroll(salt);

  spc.load(textio::prefix+"state");

  Move::TranslateRotateCluster gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::ParallelTempering pt(mcp,pot,spc,mpi);   //temper move
  
  Analysis::LineDistribution<float,unsigned long int> rdf(0.5);
  Analysis::LineDistributionNorm<float,unsigned long int> saltdistr(salt.size(), 0.2);
  
  Analysis::TwobodyForceDirect tfd(mcp, pot, spc, pol[0], pol[1], salt);
  Analysis::TwobodyForceMidp tfm(mcp, pot, spc, pol[0], pol[1], salt, &saltdistr);
  
  gmv.setMobile(salt); // specify where to look for clustered ions
  mv.setGroup(salt);   // specify atomic particles to be moved

  gmv.directions[ pol[0].name ].x()=0; // do not move in x
  gmv.directions[ pol[0].name ].y()=0; // do not move in y
  gmv.directions[ pol[0].name ].z()=0; // do move in z
  gmv.directions[ pol[1].name ].x()=0; // do not move in x
  gmv.directions[ pol[1].name ].y()=0; // do not move in y
  gmv.directions[ pol[1].name ].z()=0; // do move in z
  
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  mpi.cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");
  
  int macronum = 0;
  while ( loop.macroCnt() ) {  // Markov chain 
    macronum++;
    while ( loop.microCnt() ) {
      
      if (slp_global()>pxtc)
        xtc.save(textio::prefix+"out.xtc", spc.p);
      
      int k,i=rand() % 3;
      switch (i) {
        case 0:
          mv.setGroup(salt);
          sys+=mv.move( salt.size()/2+1 );
          break;
        case 1:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++) {
              double r=spc.geo->dist(i->cm,j->cm);
              if (r<rdf.maxdist)
                rdf(r)++;
            }
          break;
        case 2:
          sys+=pt.move();
          break;
      }
      
      for (auto j=salt.begin(); j!=salt.end()-1; j++) {
        saltdistr(spc.p[(*j)].z())++;
      }
      
      tfd.calc();
      tfm.calc();
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );

    rdf.save(textio::prefix+"rdf_p2p.dat");
    
    string strsaltdistr("");
    strsaltdistr += textio::prefix;
    strsaltdistr += "_";
    strsaltdistr += to_string(macronum);
    strsaltdistr += "saltdistr.dat";
    
    saltdistr.save(strsaltdistr);
    saltdistr.save(textio::prefix+"saltdistr.dat");
    
    string strtfd("");
    strtfd += textio::prefix;
    strtfd += "_";
    strtfd += to_string(macronum);
    strtfd += "dforce.dat";
    
    tfd.save(strtfd);
    tfd.save(textio::prefix+"dforce.dat");
    
    string strtfm("");
    strtfm += textio::prefix;
    strtfm += "_";
    strtfm += to_string(macronum);
    strtfm += "mpforce.dat";
    
    tfm.save(strtfm);
    tfm.save(textio::prefix+"mpforce.dat");
    
    string strcur("");
    strcur += textio::prefix;
    strcur += "_";
    strcur += to_string(macronum);
    strcur += "energy.dat";
    std::ofstream f(strcur.c_str());
    f.precision(10);
    if (f) {
      f << sys.current() << endl;
    }
    
    spc.save(textio::prefix+"state");
    mpi.cout << loop.info() << spc.info() << sys.info() << mv.info() << gmv.info() << pt.info();
    mpi.cout << loop.timing() << tfd.info() << tfm.info() ;
  } // end of macro loop

  pqr.save(textio::prefix+"confout.pqr", spc.p);

  mpi.cout << loop.info() << spc.info() << sys.info() << mv.info() << gmv.info() << pt.info() << tfd.info() << tfm.info();
  mpi.cout << pol[0].info();
  mpi.cout << pol[0].charge(spc.p) << endl;
}
