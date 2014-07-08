#include <faunus/faunus.h>
using namespace Faunus;

typedef std::vector<Group*> Tvec;
typedef Space<Geometry::Cuboid> Tspace;
typedef Faunus::Analysis::PenaltyFunction2D<double> Tpenalty;

  struct myenergy : public Energy::Energybase<Tspace> { //custom energy class
    public:
      double i_external(const p_vec &p, int i)  { //pot. on particle
        double s= 1 + std::sin(2*pc::pi*p[i].x()) + std::cos(2*pc::pi*p[i].y());
        if (p[i].x() >=-2.00 && p[i].x() <=-1.25) return 1*s;
        if (p[i].x() >=-1.25 && p[i].x() <=-0.25) return 2*s;
        if (p[i].x() >=-0.25 && p[i].x() <= 0.75) return 3*s;
        if (p[i].x() >= 0.75 && p[i].x() <= 1.75) return 4*s;
        if (p[i].x() >= 1.75 && p[i].x() <= 2.00) return 5*s;
        return pc::infty;

      }
      double g_external(const p_vec &p, Group &g) { //pot. on group
        double u=0;
        for (auto i : g)
          u+=i_external(p, i);
        return u; // in kT
      }
      string _info() { return "myenergy"; }       //mandatory info 
  };

struct coordinates { //function that defines the reaction coordinates
  std::pair<double,double> operator()(Tspace *spc, const p_vec &p, Tvec &gv) {
    return std::make_pair(p[gv[0]->front()].x(),p[gv[0]->front()].y());
  }
};

int main() {
  Faunus::MPI::MPIController mpi;                  //init MPI
  InputMap mcp(textio::prefix+"penalty.input");    //read input file
  MCLoop loop(mcp);                                //handle mc loops
  auto pot = myenergy()                            //our custom potential!
    + Energy::PenaltyEnergy<Tspace,Tpenalty,coordinates>(mpi, mcp);
  auto penalty = &pot.second;                      //penalty function
  Tspace spc(mcp);                                 //create simulation space
  spc.insert( PointParticle() );                   //insert a single particle
  Group mygroup(0,0);                              //group with single particle
  spc.enroll(mygroup);                             //tell space about group
  EnergyDrift sys;                                 //class for tracking system energy drifts
  Move::AtomicTranslation<Tspace> trans(mcp,pot,spc); //translational move
  trans.setGroup(mygroup);                            //set translation group
  Tvec gvec;                                          //vector of pointers to penalized groups 
  gvec.push_back(&mygroup);
  penalty->setPenalized(&spc,gvec);
  penalty->pf.load(textio::prefix+"penalty");
  sys.init( Energy::systemEnergy(spc,pot,spc.p) ); // store total energy
  int seed_value = mcp.get<int>("seed_value",-13);
  slp_global.seed(seed_value);
  int sweeps = 0;
  while ( loop.macroCnt() ) {                 //start markov chain
    while ( loop.microCnt() ) {
      sys += trans.move();                    //translate particle
      ++sweeps;
      coordinates c;
      std::pair<double,double> coor = c(&spc,spc.p,gvec);
      sys += penalty->pf.update(coor.first,coor.second); //increment the histogram and/or update the penalty function
    }
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?
    // save to disk
    if (sweeps == 2.5e5) penalty->pf.save(textio::prefix+"penalty_init"); // initial histogram obtained without penalty
    else penalty->pf.save(textio::prefix+"penalty");
  }
  mpi.cout << trans.info() + loop.info() + penalty->pf.info() + sys.info();
  // perform unit 
  UnitTest test(mcp);
  trans.test(test);
  penalty->pf.test(test);
  sys.test(test);
  return test.numFailed();
  mpi.cout << trans.info() + loop.info() + penalty->pf.info() + sys.info();
}
/**
  @page example_penalty Example: Penalty Function

  This is a 2D version of the example used for parallel tempering in the book

  - _Understanding Molecular Simulation_ by Frenkel and Smit, 2nd edition.
  p.391 - Case Study 21, _Parallel Tempering of a Single Particle_.

  We simulate a single particle in an oscillating potential and use the 
  penalty function method to overcome energy barriers. 
  Below is the two-dimensional distribution functions sampled with and without
  the penalty function method, demonstrating how the penalty function 
  (at the bottom right) leads to flat distribution functions in rough energy landscapes.

  ![Particle distribution functions in an 2D oscillating field](penalty.png)

  To run this example, make sure faunus is compiled with MPI enabled:

  ~~~
  $ cmake . -DENABLE_MPI=on
  $ make
  $ cd src/examples
  $ ./penalty.run mpirun
  ~~~

  The penalty function routine in Faunus is general and implemented in
  MPI using a master-slave scheme. Each system has its own rank and random seed.
  Samplings of the configurational space from all processes are merged by 
  periodically summing up the two-dimensional distribution functions.
  The input files, `mpi$rank.penalty.input`, for this example look like this:

  ~~~
  loop_macrosteps         100       # 2 * number of updates of penalty function
  loop_microsteps         250000    # number of moves between printing histograms
  loop_update             500000    # number of moves between updates
  hist_size               2000      # must be >= max number of points in the histogram (i.e. 41x41=1681)
  res1                    0.1       # resolution of one reaction coordinate
  res2                    0.1       # resolution of the second reaction coordinate
  cuboid_len              4         # box side length Angstrom
  mv_particle_genericdp   0.5       # translational displacement [Angstrom]
  seed_value              $seed     # random seed
  ~~~

  penalty.cpp
  ==========

  @includelineno examples/penalty.cpp

*/
