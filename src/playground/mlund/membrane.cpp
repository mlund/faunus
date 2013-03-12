#include <faunus/faunus.h>
#include <faunus/membrane.h>

using namespace Faunus;
using namespace Faunus::Potential;

class PenaltyEnergyBase : public Energy::Energybase {
  private:
    string _info() { return "penalty\n"; }
    double* _trialcoord, _coord;
    Space* spcPtr;
    Group *g1Ptr, *g2Ptr;
  public:
    Analysis::PenaltyFunction<double> f;
    void setCoord(Group *g1, Group *g2) {
      assert(g1!=nullptr && g2!=nullptr);
      g1Ptr=g1;
      g2Ptr=g2;
    }
    PenaltyEnergyBase(InputMap &in, Space &s, string pfx="penalty_") :
      f(in.get<double>(pfx+"du",0),1e5,4,0.5) { //!< Constructor - read input parameters
      name="Group-group mass center penalty";
      spcPtr=&s;
    }
    double i2all(const p_vec &p, int i) {
      if (g1Ptr->find(i) || g2Ptr->find(i)) {
        double r;
        if (&p==&spcPtr->trial) {   // p is trial
          g1Ptr->cm_trial = Geometry::massCenter(*spcPtr->geo,p,*g1Ptr); // if GroupArray, keep mass center up2date
          g2Ptr->cm_trial = Geometry::massCenter(*spcPtr->geo,p,*g2Ptr); // if GroupArray, keep mass center up2date
          r = abs(  spcPtr->geo->vdist(g1Ptr->cm_trial, g2Ptr->cm_trial).z()  );
        } else { 
          g1Ptr->cm = Geometry::massCenter(*spcPtr->geo,p,*g1Ptr); // if GroupArray, keep mass center up2date
          g2Ptr->cm = Geometry::massCenter(*spcPtr->geo,p,*g2Ptr); // if GroupArray, keep mass center up2date
          r = abs(  spcPtr->geo->vdist(g1Ptr->cm, g2Ptr->cm).z()  );
        }
        return f(r);
      }
      return 0;
    }
    double g2g(const p_vec &p, Group &g1, Group &g2) FOVERRIDE {
      // test if g1 AND g2 are part of the penalized groups
      if ( ( g1Ptr->find(g1.front()) && g2Ptr->find(g2.front()) )
          || ( g1Ptr->find(g2.front()) && g2Ptr->find(g1.front()) ) ) {
        double r;
        if (&p==&spcPtr->trial) {   // p is trial
          g1Ptr->cm_trial = Geometry::massCenter(*spcPtr->geo,p,*g1Ptr); // if GroupArray, keep mass center up2date
          g2Ptr->cm_trial = Geometry::massCenter(*spcPtr->geo,p,*g2Ptr); // if GroupArray, keep mass center up2date
          r = abs(  spcPtr->geo->vdist(g1Ptr->cm_trial, g2Ptr->cm_trial).z()  );
        } else {
          g1Ptr->cm = Geometry::massCenter(*spcPtr->geo,p,*g1Ptr);
          g2Ptr->cm = Geometry::massCenter(*spcPtr->geo,p,*g2Ptr);
          r = abs(  spcPtr->geo->vdist(g1Ptr->cm, g2Ptr->cm).z()  );
        }
        return f(r);
      }
      return 0;
    }
};

namespace Faunus {
  namespace Potential {
    // custom pair potential that uses CosAttrac/WCA for lipid and Lennard-Jones for the rest
    class LJWCA : public WeeksChandlerAndersen {
      private:
        particle::Tid tailid, headid, hhisid;
        typedef opair<particle::Tid> Tpair;
        Tpair _tailPair, _headPair, _headtailPair;
      public:
        CosAttract cosattract;
        ChargeNonpolar qnop;
        LJWCA(InputMap &in) : WeeksChandlerAndersen(in), cosattract(in), qnop(in) {
          name="Lipids: WCA+CosAttract, HHIS-Tail: ChargeNonpolar, Rest: LJ";
          tailid=atom["TL"].id;
          headid=atom["HD"].id;
          hhisid=atom["HHIS"].id;
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

  cout << textio::splash();           // show faunus banner and credits
  InputMap mcp("membrane.input");   //read input file

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
#ifdef PENALTY
  auto penalty   = pot.create( PenaltyEnergyBase(mcp,spc) );
  //Analysis::Table2D<float, Analysis::PenaltyFunction<double> > mik;
#endif

  // Markov moves and analysis
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::TranslateRotate gmv(mcp,pot,spc), gmvpol(mcp,pot,spc,"polymer");
  Move::Pivot pivot(mcp, pot, spc);
  Move::Reptation rep(mcp, pot, spc);
  Move::CrankShaft crank(mcp, pot, spc);
  Move::SwapMove tit(mcp,pot,spc);
  Analysis::PolymerShape shape;
  Analysis::RadialDistribution<> rdf(0.5);
  Analysis::LineDistribution<> hist_pepmem(0.5);
  Analysis::Table2D<float, Average<double> > z_pepmem(0.5), end2end(0.5);

  // Load membrane
  DesernoMembrane<Tgeometry> mem(mcp,pot,spc, nonbonded->pairpot.first, nonbonded->pairpot.first.cosattract);

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
  for (int i=pol.front(); i<pol.back(); i++)
    bonded->add(i, i+1, Potential::Harmonic(k,req));   // add bonds

#ifdef PENALTY
  penalty->setCoord(&mem.lipids, &pol);
#endif

  //for (size_t i=0; i<spc.p.size(); i++){
  //  spc.p[i].charge=atom[spc.p[i].id].charge;
  //  spc.trial[i].charge=spc.p[i].charge;
  //}
  tit.findSites(spc.p);  // search for titratable sites

  spc.load("state");                                     // load old config. from disk (if any)
  pqr.save("initial.pqr", spc.p);
#ifdef PENALTY
  penalty->f.load("penalty.dat");
  penalty->f.save("penalty_initial.dat");
#endif

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k=mem.lipids.sizeMol(); //number of lipids
      int i=slp_global.rand() % 4;
      Group g;
      switch (i) {
        case 0:
          if (slp_global()>0.5) {
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
          if (slp_global()>0.5) {
            gmvpol.setGroup(pol);
            sys+=gmvpol.move();
          } else {
            //crank.setGroup(pol);
            //sys+=crank.move();
          }
          break;
        case 30:
          if (slp_global()>0.9)
            sys+=tit.move();
          break;
      }

      // peptide-membrane distribution
      double d = spc.geo->vdist(pol.cm, mem.lipids.massCenter(spc)).z();
      d=abs(d);
      hist_pepmem(d)++;
      z_pepmem(d) += pol.charge(spc.p);
      end2end(d) += spc.geo->dist( spc.p[pol.front()], spc.p[pol.back()] );
#ifdef PENALTY
      sys+=penalty->f.update(d);
#endif

      // gromacs trajectory
      if ( slp_global()<0.01 )
        xtc.save("traj.xtc", spc);

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current

    // save to disk
    rdf.save("rdf.dat");
    pqr.save("confout.pqr", spc.p);
    spc.save("state");
    hist_pepmem.save("pepmem.dat");
    z_pepmem.save("zpepmem.dat");
    end2end.save("end2end.dat");

    cout << loop.timing();

#ifdef PENALTY
    penalty->f.save("penalty.dat");
    cout << penalty->f.info();
#endif
  } // end of macro loop

  // perform unit tests
  gmv.test(test);
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() << sys.info() << mv.info() << gmv.info() << gmvpol.info() << rep.info()
    << crank.info() << shape.info() << test.info() << spc.info() << tit.info();

  return test.numFailed();
}
