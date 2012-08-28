/* 
 * Parallel tempering example from the book of Frenkel and Smit, 2nd ed,
 * page 391 - Case Study 21 "Parallel Tempering of a Single Particle"
 */

#include <faunus/faunus.h>
using namespace Faunus;

class myenergy : public Energy::Energybase {
  private:
    string _info() { return "myenergy"; }
  public:
    double Tscale; // reduced temperature
    // external potential on i'th particle
    double i_external(const p_vec &p, int i)  {
      double s=( 1+std::sin(2*pc::pi*p[i].x) ) / Tscale;
      if (p[i].x>=-2 && p[i].x<=-1.25) return 1*s;
      if (p[i].x>=-1.25 && p[i].x<=-0.25) return 2*s;
      if (p[i].x>=-0.25 && p[i].x<=0.75) return 3*s;
      if (p[i].x>=0.75 && p[i].x<=1.75) return 4*s;
      if (p[i].x>=1.75 && p[i].x<=2) return 5*s;
      return pc::infty;
    }
    // external potential of group of particles
    double g_external(const p_vec &p, Group &g) {
      double u=0;
      for (auto i : g)
        u+=i_external(p, i);
      return u; // in kT
    }
};

int main() {
  Faunus::MPI::MPIController mpi;// init MPI
  InputMap mcp(textio::prefix+"temper.input"); // read input file
  MCLoop loop(mcp);              // handle mc loops
  EnergyDrift sys;               // track system energy drifts
  Geometry::Cuboid geo(mcp);     // set geometry
  myenergy pot;                  // our own potential! see above.
  Space spc(geo);                // create simulation space
  spc.insert( PointParticle() ); // insert a single particle
  Group mygroup(0,0);            // group with our single particle
  spc.enroll(mygroup);           // make sure space knows about this group

  Move::ParallelTempering temper(mcp,pot,spc,mpi); // temper move
  Move::AtomicTranslation trans(mcp,pot,spc);      // particle translation move
  Analysis::LineDistribution<float,int> dist(0.05);// distribution function

  trans.setGroup(mygroup);            // translate particles from mygroup
  pot.Tscale = mcp.get("Tscale",1.0); // get temperature from user input

  sys.init( pot.g_external(spc.p, mygroup) ); // store initial energy

  mpi.cout << spc.info() << loop.info();

  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      sys+=trans.move();      // and translate it
      dist(spc.p[0].x)++;     // add x-pos of first particle to histogram
    }
    sys+=temper.move();       // do temper moves
    sys.checkDrift( pot.g_external(spc.p, mygroup) ); // energy drift?
    mpi.cout << loop.timing();// print progress
  }

  dist.save(textio::prefix+"dist"); // save x-position distribution
  mpi.cout << trans.info() << temper.info() << sys.info();

  // save external potential to disk
  std::ofstream o(textio::prefix+"pot");
  for (spc.p[0].x=-2; spc.p[0].x<2; spc.p[0].x+=0.05)
    o << spc.p[0].x << " " << pot.i_external(spc.p, 0) << endl;
}
