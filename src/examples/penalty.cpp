#include <faunus/faunus.h>
using namespace Faunus;

typedef Space<Geometry::Cuboid> Tspace;

struct myenergy : public Energy::Energybase<Tspace> { //custom energy class
  public:
    typedef typename Tspace::ParticleVector Tpvec;
    double i_external(const Tpvec &p, int i) FOVERRIDE { //pot. on particle
      double s= 1 + std::sin(2*pc::pi*p[i].x()) + std::cos(2*pc::pi*p[i].y());
      if (p[i].x() >=-2.00 && p[i].x() <=-1.25) return 1*s;
      if (p[i].x() >=-1.25 && p[i].x() <=-0.25) return 2*s;
      if (p[i].x() >=-0.25 && p[i].x() <= 0.75) return 3*s;
      if (p[i].x() >= 0.75 && p[i].x() <= 1.75) return 4*s;
      if (p[i].x() >= 1.75 && p[i].x() <= 2.00) return 5*s;
      return 1e10;
    }
    auto tuple() -> decltype(std::make_tuple(this)) {
      return std::make_tuple(this);
    }
    double g_external(const Tpvec &p, Group &g) FOVERRIDE { //pot. on group
      double u=0;
      for (auto i : g)
        u+=i_external(p, i);
      return u; // in kT
    }
    string _info() { return "myenergy"; }       //mandatory info 
};

struct coordinates { //function that defines the reaction coordinates
  static Group* pg; // penalized group
  std::pair<double,double> operator()(const typename Tspace::ParticleVector &p) {
    return std::make_pair(p[pg->front()].x(),p[pg->front()].y());
  }
};

Group* coordinates::pg;

int main() {
  // In MPI version:
  // Faunus::MPI::MPIController mpi; // init MPI

  InputMap mcp("penalty.json"); // read input file
  MCLoop loop(mcp); // class for handling mc loops

  auto pot
    = myenergy()                                      // our custom potential!
    + Energy::PenaltyEnergy<Tspace,coordinates>(mcp);

  // In MPI version:
  // + Energy::PenaltyEnergy<Tspace,coordinates>(mpi, mcp);

  auto penalty = std::get<1>( pot.tuple() );

  Tspace spc(mcp);                                    // create simulation space
  auto myparticle = spc.findMolecules("myparticle");
  Group mygroup(myparticle.front()->front(), myparticle.back()->back());
  mygroup.setMolSize(1);
  coordinates::pg = &mygroup;                         // set penalized group
  EnergyDrift sys;                                    // class for tracking system energy drifts
  
  // In MPI version:
  // slump.setDiscard(mpi.rank()+1);
  penalty->load("pf_");

  // Markov moves
  Move::Propagator<Tspace> mv(mcp, pot, spc);

  Table3D<double,double> histo(0.1, 0.1, Table3D<double,double>::HISTOGRAM);

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // store total energy

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      sys+=mv.move();  // translate particle
      ++histo(spc.p[mygroup.front()].x(),spc.p[mygroup.front()].y());
    }
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?
  }
  auto it_min = histo.min();
  histo.save("histo"+std::to_string(histo.getMap().size()),1./it_min->second);
  penalty->save("pf_");

  cout << loop.info() + mv.info() + penalty->info() + sys.info();

  // perform unit 
  if (penalty->penalty_update(true)) {
  UnitTest test(mcp);
  mv.test(test);
  sys.test(test);
  penalty->test(test);
  cout << test.info();
  return test.numFailed();
  }
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

  ![Particle distribution functions and penalty function in an 2D oscillating field](penalty.png)

  To run this example, make sure faunus is compiled with MPI enabled:

  ~~~
  $ cmake . -DENABLE_MPI=on
  $ make
  $ cd src/examples
  $ ./penalty.run mpirun
  ~~~

  A gnuplot script to generate the plots of the probability distributions and 
  penalty function is provided:

  ~~~
  $ gnuplot penalty.gnu
  ~~~

  The penalty function routine in Faunus is general and implemented in
  MPI using a master-slave scheme. Each system has its own rank and random seed.
  Samplings of the configurational space from all processes are merged by 
  periodically summing up the two-dimensional distribution functions.
  The parameters to be set in the input file are following:

  ~~~
  f0              1.0       # initial increment to the penalty function
  scale           0.8       # factor by which f0 is scaled
  update          1e5       # number of MC sweeps before scaling f0
  bw1             0.1       # bin width of 1st coordinate
  bw2             0.1       # bin width of 2nd coordinate
  lo1             -2.0      # lower limit of 1st coordinate
  hi1             2.0       # upper limit of 1st coordinate
  lo2             -2.0      # lower limit of 2nd coordinate
  hi2             2.0       # upper limit of 2nd coordinate
  ~~~

  penalty.cpp
  ==========

  @includelineno examples/penalty.cpp

*/
