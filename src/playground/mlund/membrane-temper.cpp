#include <faunus/faunus.h>
#include <faunus/membrane.h>

using namespace Faunus;
using namespace Faunus::Potential;

namespace Faunus {
  namespace Potential {
    // custom pair potential that uses CosAttrac/WCA for lipid and Lennard-Jones for the rest
    class LJWCA : public WeeksChandlerAndersen {
      private:
        particle::Tid tailid, headid, hhisid;
        typedef pair_permutable<particle::Tid> Tpair;
        Tpair _tailPair, _headPair, _headtailPair;
      public:
        CosAttract cosattract;
        ChargeNonpolar qnop;
        LJWCA(InputMap &in) : WeeksChandlerAndersen(in), cosattract(in), qnop(in) {
          name="Lipids: WCA+CosAttract, HHIS-Tail: ChargeNonpolar, Rest: LJ";
          tailid=atom["TL"].id;
          headid=atom["HD"].id;
          _tailPair=Tpair(tailid,tailid);
          _headPair=Tpair(headid,headid);
          _headtailPair=Tpair(headid,tailid);
        }
        inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
          Tpair p(a.id,b.id);
          if (p==_tailPair) // two tails
            return cosattract(a,b,r2) + WeeksChandlerAndersen::operator()(a,b,r2);
          if (p==_headPair) // two heads
            return WeeksChandlerAndersen::operator()(a,b,r2);
          if (p==_headtailPair) { // head-tail (head could be charged)
            return WeeksChandlerAndersen::operator()(a,b,r2) + qnop(a,b,r2);
          }
          double u=WeeksChandlerAndersen::Tbase::operator()(a,b,r2); // default, base=LJ
          if (a.id==tailid || b.id==tailid) // if either a or b is tail, add possible non-polar repulsion
            u+=qnop(a,b,r2);
          return u;
        }
    };
  }
}

typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
typedef CombinedPairPotential<LJWCA, DebyeHuckel> Tpairpot;
//typedef LJWCA Tpairpot;

int main() {
  Faunus::MPI::MPIController mpi;                //init MPI
  InputMap mcp(textio::prefix+"membrane.input"); 
  mpi.cout << textio::splash();           // show faunus banner and credits

  MCLoop loop(mcp);                   // class for handling mc loops
  FormatPQR pqr;                      // PQR structure file I/O
  FormatAAM aam;
  FormatXTC xtc(1000);                // XTC gromacs trajectory format
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);                 // class for unit testing

  // Energy functions and space
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  double custeps = mcp.get<double>("custeps", 0);
  mpi.cout << "# Custom epsilon = " << custeps << "\n";
  nonbonded->pairpot.first.customEpsilon(atom["TL"].id, atom["ARG"].id, custeps);
  nonbonded->pairpot.first.customEpsilon(atom["TL"].id, atom["HARG"].id, custeps);
  nonbonded->pairpot.first.customEpsilon(atom["TL"].id, atom["NTR"].id, custeps);
  nonbonded->pairpot.first.customEpsilon(atom["TL"].id, atom["HNTR"].id, custeps);
  nonbonded->pairpot.first.customEpsilon(atom["TL"].id, atom["CTR"].id, custeps);
  nonbonded->pairpot.first.customEpsilon(atom["TL"].id, atom["HCTR"].id, custeps);
  //nonbonded->pairpot.customEpsilon(atom["TL"].id, atom["ARG"].id, custeps);

  // Markov moves and analysis
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::TranslateRotate gmv(mcp,pot,spc), gmvpol(mcp,pot,spc,"polymer");
  Move::Pivot pivot(mcp, pot, spc);
  Move::Reptation rep(mcp, pot, spc);
  Move::CrankShaft crank(mcp, pot, spc);
  Move::SwapMove tit(mcp,pot,spc);
  Move::ParallelTempering pt(mcp,pot,spc,mpi);   //temper move
  Analysis::PolymerShape shape;
  Analysis::RadialDistribution<> rdf(0.5);
  Analysis::LineDistribution<> hist_pepmem(0.5);
  Analysis::Table2D<float, Average<double> > z_pepmem(0.5), end2end(0.5);

  // Load membrane
  DesernoMembrane<Tgeometry> mem(mcp,pot,spc, nonbonded->pairpot.first, nonbonded->pairpot.first.cosattract);
  //DesernoMembrane<Tgeometry> mem(mcp,pot,spc, nonbonded->pairpot, nonbonded->pairpot.cosattract);

  // Load peptide
  string polyfile = mcp.get<string>("polymer_file", "");
  double req    = mcp.get<double>("polymer_eqdist", 0);
  double k      = mcp.get<double>("polymer_forceconst", 0);
  aam.load(polyfile);                                 // load polymer structure into aam class
  Geometry::FindSpace f;                              // class for finding empty space in container
  f.find(*spc.geo, spc.p, aam.particles());           // find empty spot
  GroupMolecular pol = spc.insert( aam.particles() ); // Insert particles into Space and return matching group
  pol.name="peptide";                                 // Give the polymer an arbitrary name
  spc.enroll(pol);                                    // All groups need to be enrolled in the Space
  //for (int i=pol.front(); i<pol.back(); i++)
  //  bonded->add(i, i+1, Potential::Harmonic(k,req));   // add bonds

  for (size_t i=0; i<spc.p.size(); i++){
    spc.p[i].charge=atom[spc.p[i].id].charge;
    spc.trial[i].charge=spc.p[i].charge;
  }
  tit.findSites(spc.p);  // search for titratable sites

  spc.load(textio::prefix+"state");                                     // load old config. from disk (if any)
  pqr.save(textio::prefix+"initial.pqr", spc.p);
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  mpi.cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  std::ofstream distfile(textio::prefix+"dist.dat", std::ios_base::app);

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k=mem.lipids.sizeMol(); //number of lipids
      int i=slp_global.rand() % 4;
      Group g;
      switch (i) {
        case 0:
          if (slp_global.randOne()>0.5) {
            mv.setGroup(mem.lipids);
            sys+=mv.move( mem.lipids.size() ); // translate lipid monomers
          } else {
            //mv.setGroup(pol);
            //sys+=mv.move( pol.size() )       ; // translate peptide monomers
            //pol.setMassCenter(spc);
          }
          break;
        case 1:
          while (k-->0) {
            i=mem.lipids.randomMol();   // pick random lipid molecule
            g=mem.lipids[i];
            g.name="subgrouplipid";
            g.setMassCenter(spc);       // mass center needed for rotation
            gmv.setGroup(g);            // tell what to move
            sys+=gmv.move();            // translate/rotate polymers
          }
          break;
        case 2:
          if (slp_global.randOne()>0.5) {
            gmvpol.setGroup(pol);
            sys+=gmvpol.move();
          } else {
            //crank.setGroup(pol);
            //sys+=crank.move();
          }
          break;
        case 3:
          sys+=pt.move();
          break;
        case 4:
          sys+=tit.move();
      }

      // peptide-membrane distribution
      double d = spc.geo->vdist(pol.cm, mem.lipids.massCenter(spc)).z();
      //d=abs(d);
      hist_pepmem(d)++;
      z_pepmem(d) += pol.charge(spc.p);
      end2end(d) += spc.geo->dist( spc.p[pol.front()], spc.p[pol.back()] );

      // gromacs trajectory
      if ( slp_global.randOne()<0.005 ) {
        distfile << d << "\n";
        xtc.save(textio::prefix+"traj.xtc", spc);
      }

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current

    // save to disk
    rdf.save(textio::prefix+"rdf.dat");
    pqr.save(textio::prefix+"confout.pqr", spc.p);
    spc.save(textio::prefix+"state");
    hist_pepmem.save(textio::prefix+"pepmem.dat");
    z_pepmem.save(textio::prefix+"zpepmem.dat");
    end2end.save(textio::prefix+"end2end.dat");

    mpi.cout << loop.timing();

  } // end of macro loop

  // print information
  mpi.cout << loop.info() << sys.info() << mv.info() << gmv.info() << gmvpol.info() << rep.info()
    << crank.info() << shape.info() << test.info() << spc.info() << tit.info() << pt.info();

  return test.numFailed();
}
