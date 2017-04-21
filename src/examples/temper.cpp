#include <faunus/faunus.h>
using namespace Faunus;
typedef Space<Geometry::Cuboid> Tspace;

template<class Tspace>
struct myenergy : public Energy::Energybase<Tspace> {  //custom energy class
  typedef typename Tspace::ParticleVector Tpvec;
  double Tscale;                               //reduced temperature
  double i_external(const Tpvec &p, int i) override {  //pot. on particle
    double s=( 1+std::sin(2*pc::pi*p[i].x()) ) / Tscale;
    if (p[i].x() >=-2.00 && p[i].x() <=-1.25) return 1*s;
    if (p[i].x() >=-1.25 && p[i].x() <=-0.25) return 2*s;
    if (p[i].x() >=-0.25 && p[i].x() <= 0.75) return 3*s;
    if (p[i].x() >= 0.75 && p[i].x() <= 1.75) return 4*s;
    if (p[i].x() >= 1.75 && p[i].x() <= 2.00) return 5*s;
    return pc::infty;
  }
  double g_external(const Tpvec &p, Group &g) override {//pot. on group
    double u=0;
    for (auto i : g)
      u+=i_external(p, i);
    return u; // in kT
  }
  string _info() override { return "myenergy"; }         //mandatory info 
  auto tuple() -> decltype(std::make_tuple(this))
  {
      return std::make_tuple(this);
  }

};

int main() {
  Faunus::MPI::MPIController mpi;               //init MPI

  InputMap in( textio::prefix+"temper.json" );   // read input file
  Tspace spc( in );                             //create simulation space

  myenergy<Tspace> pot;                         //our custom potential!
  pot.Tscale = in["system"]["Tscale"] | 1.0;    //temperature from input (default: 1)

  Table2D<double,int> dst(0.05);                //distribution func.
  Move::Propagator<Tspace> mv(in,pot,spc,&mpi); //MC moves

  mpi.cout << spc.info();                       //print initial info
                                               
  MCLoop loop( in );                            //handle mc loops
  while ( loop[0] ) {                           //start markov chain
    while ( loop[1] ) {                        
      mv.move();                                //translate particle
      dst( spc.p[0].x() )++;                    //update histogram
    }                                          
    mpi.cout << loop.timing();                  //print progress
  }                                            

  UnitTest test( in );                          //unit testing
  mv.test( test );                          

  dst.save( textio::prefix+"dist" );            //save histogram

  mpi.cout << mv.info() << test.info();

  return test.numFailed();
}
/**
  @page example_temper Example: Parallel Tempering

  This is an example of parallel tempering taken from the book

  - _Understanding Molecular Simulation_ by Frenkel and Smit, 2nd edition.
  p.391 - Case Study 21, _Parallel Tempering of a Single Particle_.

  We simulate a single particle in an oscillating
  potential and use tempering in temperature (you may use any
  parameter in the Hamiltonian) to overcome energy barriers.
  Below is the one-dimensional distribution function sampled with and without
  tempering, demonstrating how the high T replicas help low T systems to
  pass steep barriers.

  ![Particle distribution function in an oscillating field](temper.png)

  To run this example, make sure faunus is compiled with MPI enabled:

  ~~~
  $ cmake . -DENABLE_MPI=on
  $ make example_temper
  $ cd src/examples
  $ ./temper.run mpirun
  ~~~

  The tempering routine in Faunus is general and implemented in
  MPI where each replica has its own rank. The temper parameter(s) is
  given simply by giving different input files for each replica.

  temper.cpp
  ==========

  @includelineno examples/temper.cpp
*/

