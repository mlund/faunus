#include <faunus/faunus.h>
#include <tclap/CmdLine.h>

using namespace Faunus;

typedef Geometry::PeriodicCylinder Tgeometry;
typedef Potential::CombinedPairPotential<Potential::DebyeHuckel, Potential::LennardJones> Tpairpot;

int main(int argc, char** argv) {
  Faunus::MPI::MPIController mpi;
  InputMap mcp(textio::prefix+"input");
  mpi.cout << textio::splash();
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatTopology top;
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  //auto constrain = pot.create( Energy::MassCenterConstrain(pot.getGeometry()) );
  Space spc( pot.getGeometry() );

  // Add molecular species
  int cnt=0;
  int N1=mcp.get("polymer1_N",0);
  int N2=mcp.get("polymer2_N",0);
  vector<GroupMolecular> pol( N1+N2);
  for (auto &g : pol) {
    cnt++;
    string polyfilekey = (cnt>N1) ? "polymer2_file" : "polymer1_file";
    aam.load( mcp.get<string>(polyfilekey, "") );
    Geometry::FindSpace f;
    f.dir.x()=0; // put mass center
    f.dir.y()=0; //   at [x,y,z] = [0,0,random]
    if (f.find(*spc.geo, spc.p, aam.p )) {
      g = spc.insert( aam.p );
      g.name=mcp.get<string>(polyfilekey, "");
      spc.enroll(g);
    } else
      return 1;
  }

  //constrain->addPair(pol[0], pol[1], 20, 150);

  // Add salt
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  //spc.enroll(salt);

  Move::SwapMove tit(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::ParallelTempering pt(mcp,pot,spc,mpi);
  mv.setGroup(salt);   // specify atomic particles to be moved

  tit.findSites(spc.p);  // search for titratable sites
  spc.load(textio::prefix+"state");

  typedef Analysis::Table2D<float, Average<double> > Tdiptable;
  Tdiptable dip1(0.25, Tdiptable::XYDATA);
  Tdiptable dip2(0.25, Tdiptable::XYDATA);
  Tdiptable dipdip(0.25, Tdiptable::XYDATA);

  gmv.directions[ pol[0].name ].x()=0; // do not move in x
  gmv.directions[ pol[0].name ].y()=0; // do not move in y
  gmv.directions[ pol[0].name ].z()=1; // do move in z
  gmv.directions[ pol[1].name ].x()=0; // do not move in x
  gmv.directions[ pol[1].name ].y()=0; // do not move in y
  gmv.directions[ pol[1].name ].z()=1; // do move in z
  Analysis::LineDistribution<float,unsigned long int> rdf(0.5);
  //rdf.maxdist=150;

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  mpi.cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 3;
      switch (i) {
        case 2:
          sys+=pt.move();
          //mv.setGroup(salt);
          //sys+=mv.move( salt.size()/2 );
          break;
        case 0:
          if (slp_global()<0.1)
            sys+=tit.move();
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
              if (r<rdf.maxdist) {
                rdf(r)++;
              }
            }
          double r=spc.geo->dist(pol[0].cm,pol[1].cm);
          Point mu1 = pol[0].dipolemoment(spc);
          Point mu2 = pol[1].dipolemoment(spc);
          dip1(r) += mu1.z()/mu1.len();
          dip2(r) += mu2.z()/mu2.len();
          dipdip(r) += mu1.z()*mu2.z() / mu1.len() / mu2.len();
          break;
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );
    rdf.save(textio::prefix+"rdf_p2p.dat");
    mpi.cout << loop.timing();
  } // end of macro loop

  pqr.save(textio::prefix+"confout.pqr", spc.p);
  spc.save(textio::prefix+"state");
  dip1.save(textio::prefix+"dip1.dat");
  dip2.save(textio::prefix+"dip2.dat");
  dipdip.save(textio::prefix+"dipdip.dat");

  mpi.cout << loop.info() << spc.info() << sys.info() << mv.info() << gmv.info()
    << tit.info() << pol[0].info() << pol[0].charge(spc.p) << pt.info();
}
