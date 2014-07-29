#include <faunus/faunus.h>
using namespace Faunus;

typedef Space<Geometry::Cuboid> Tspace;

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
  static Group* pg; // penalized group
  std::pair<double,double> operator()(const p_vec &p) {
    return std::make_pair(p[pg->front()].x(),p[pg->front()].y());
  }
};
Group* coordinates::pg;

int main() {
  Faunus::MPI::MPIController mpi;                  //init MPI
  InputMap mcp(textio::prefix+"penalty.input");    //read input file
  MCLoop loop(mcp);                                //handle mc loops
  auto pot = myenergy()                            //our custom potential!
    + Energy::PenaltyEnergy<Tspace,coordinates>(mpi, mcp);
  auto penalty = &pot.second;                      //penalty function
  Tspace spc(mcp);                                 //create simulation space
  spc.insert( PointParticle() );                   //insert a single particle
  Group mygroup(0,0);                              //group with single particle
  spc.enroll(mygroup);                             //tell space about group
  EnergyDrift sys;                                 //class for tracking system energy drifts
  Move::AtomicTranslation<Tspace> trans(mcp,pot,spc); //translational move
  trans.setGroup(mygroup);                            //set translation group
  coordinates::pg = &mygroup;                         //set penalized group
  sys.init( Energy::systemEnergy(spc,pot,spc.p) ); // store total energy
  int seed_value = mcp.get<int>("seed_value",-13);
  slp_global.seed(seed_value);
  int sweeps = 0;
  while ( loop.macroCnt() ) {                 //start markov chain
    while ( loop.microCnt() ) {
      sys(trans.recycle());                    //translate particle
      ++sweeps;
      //increment the histogram and/or update the penalty function using waste recycling
      sys += penalty->update(penalty->coordpair,sys.weight,sys.rejection); 
    }
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?
    // save to disk
    if (sweeps == 2.5e5) penalty->save(textio::prefix+"init_"); // initial histogram obtained without penalty
    else penalty->save(textio::prefix);
  }
  // perform unit 
  UnitTest test(mcp);
  trans.test(test);
  sys.test(test);
  penalty->test(test);

  mpi.cout << trans.info() + loop.info() + penalty->info() + sys.info() + test.info();

  return test.numFailed();
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
  (at the bottom right) yields flat distribution functions in rough energy landscapes.

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
  penalty_update          500000    # number of moves between updates
  penalty_size            20000     # must be >= max number of points in the histogram (i.e. 41x41=1681)
  penalty_res1            0.1       # bin width of 1st coordinate
  penalty_res2            0.1       # bin width of 2nd coordinate
  penalty_lo1             -2.0      # lower limit of 1st coordinate
  penalty_hi1             2.0       # upper limit of 1st coordinate
  penalty_lo2             -2.0      # lower limit of 2nd coordinate
  penalty_hi2             2.0       # upper limit of 2nd coordinate
  cuboid_len              4         # box side length Angstrom
  mv_particle_genericdp   0.5       # translational displacement [Angstrom]
  seed_value              $seed     # random seed
  ~~~

  penalty.cpp
  ==========

  @includelineno examples/penalty.cpp

*/
