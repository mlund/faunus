#include <faunus/faunus.h>
using namespace Faunus;
typedef Space<Geometry::Cuboid> Tspace;

template<class Tspace>
struct myenergy : public Energy::Energybase<Tspace> {  //custom energy class
  typedef typename Tspace::ParticleVector Tpvec;
  double Tscale;                               //reduced temperature
  double i_external(const Tpvec &p, int i) FOVERRIDE {  //pot. on particle
    double s=( 1+std::sin(2*pc::pi*p[i].x()) ) / Tscale;
    if (p[i].x() >=-2.00 && p[i].x() <=-1.25) return 1*s;
    if (p[i].x() >=-1.25 && p[i].x() <=-0.25) return 2*s;
    if (p[i].x() >=-0.25 && p[i].x() <= 0.75) return 3*s;
    if (p[i].x() >= 0.75 && p[i].x() <= 1.75) return 4*s;
    if (p[i].x() >= 1.75 && p[i].x() <= 2.00) return 5*s;
    return pc::infty;
  }
  double g_external(const Tpvec &p, Group &g) FOVERRIDE {//pot. on group
    double u=0;
    for (auto i : g)
      u+=i_external(p, i);
    return u; // in kT
  }
  string _info() { return "myenergy"; }         //mandatory info 
};

int main() {
  Faunus::MPI::MPIController mpi;               //init MPI
  InputMap mcp(textio::prefix+"temper.json");   //read input file
  myenergy<Tspace> pot;                         //our custom potential!
  mcp.cd ("system");
  pot.Tscale = mcp.get("Tscale",1.0);           //temperature from input
  Tspace spc(mcp);                              //create simulation space

  Analysis::LineDistribution<> dst(.05);        //distribution func.
  Move::ParallelTempering<Tspace> pt(mcp,pot,spc,mpi);//temper move
  Move::AtomicTranslation<Tspace> trans(mcp,pot,spc); //translational move

  EnergyDrift sys;                              // class for tracking system energy drifts
  cout << "hej2" << endl;
  sys.init(Energy::systemEnergy(spc,pot,spc.p));// store initial total system energy
  cout << "hej3" << endl;

  mpi.cout << spc.info();                       //print initial info
                                               
  MCLoop loop(mcp);                             //handle mc loops
  while ( loop[0] ) {                           //start markov chain
    while ( loop[1] ) {                        
      sys+=trans.move();                        //translate particle
      dst(spc.p[0].x())++;                      //update histogram
    }                                          
    sys+=pt.move();                             //do temper move
    mpi.cout << loop.timing();                  //print progress
  }                                            

  UnitTest test( mcp );                         //unit testing
  trans.test( test );                          
  sys.test( test );
  pt.test( test );                             

  sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // calc. energy drift
  dst.save(textio::prefix+"dist");              //save histogram
  mpi.cout << trans.info() << pt.info()         //print final info
    << sys.info() << test.info();

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
  $ make
  $ cd src/examples
  $ ./temper.run mpirun
  ~~~

  The tempering routine in Faunus is general and implemented in
  MPI where each replica has its own rank. The temper parameter(s) is
  given simply by giving different input files for each replica.
  The input files, `mpi$rank.temper.input`, for this example look like this:

  ~~~~
  {
    "atomlist" : { "A" : { "dp":0.5 } },

    "moleculelist" : {
      "myparticles" : { "atoms":"A", "atomic":true, "Ninit":1 }
    },

    "moves" : {
      "temper"        : { "prob":1.0, "format":"XYZ" },
      "atomtranslate" : { "myparticles" : {} }
    },

    "system" : {
      "Tscale"       : 1.0,
      "cuboid"       : { "len" : 4. },
      "mcloop"       : { "macro":2000, "micro":10000 },
      "unittest"     : { "testfile":"mpi2.temper.test", "stable":false },
      "atomlist"     : "mpi2.temper.json",
      "moleculelist" : "mpi2.temper.json"
    }
  }
  ~~~~

  temper.cpp
  ==========

  @includelineno examples/temper.cpp
*/

