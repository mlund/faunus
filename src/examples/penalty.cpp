#include <faunus/faunus.h>

using namespace Faunus;

typedef Space<Geometry::Cuboid> Tspace;

struct myenergy : public Energy::Energybase<Tspace>
{ //custom energy class
public:
    typedef typename Tspace::ParticleVector Tpvec;

    double i_external( const Tpvec &p, int i ) override
    { //pot. on particle
        double s = 1 + std::sin(2 * pc::pi * p[i].x()) + std::cos(2 * pc::pi * p[i].y());
        if ( p[i].x() >= -2.00 && p[i].x() <= -1.25 )
            return 1 * s;
        if ( p[i].x() >= -1.25 && p[i].x() <= -0.25 )
            return 2 * s;
        if ( p[i].x() >= -0.25 && p[i].x() <= 0.75 )
            return 3 * s;
        if ( p[i].x() >= 0.75 && p[i].x() <= 1.75 )
            return 4 * s;
        if ( p[i].x() >= 1.75 && p[i].x() <= 2.00 )
            return 5 * s;
        return 1e10;
    }

    auto tuple() -> decltype(std::make_tuple(this))
    {
        return std::make_tuple(this);
    }

    double g_external( const Tpvec &p, Group &g ) override
    { //pot. on group
        double u = 0;
        for ( auto i : g )
            u += i_external(p, i);
        return u; // in kT
    }

    string _info() override { return "myenergy"; }       //mandatory info
};

int main()
{
    // For the MPI version add the following line:
    // Faunus::MPI::MPIController mpi; // init MPI

    InputMap mcp("penalty.json"); // read input file
    MCLoop loop(mcp); // class for handling mc loops
    Tspace spc(mcp);
    auto pot
        = myenergy()                                      // our custom potential
            + Energy::PenaltyEnergy<Tspace>(mcp, spc); // to be subsituted with
    // + Energy::PenaltyEnergy<Tspace>(mpi,mcp,spc); // in the MPI version
    auto penalty = std::get<1>(pot.tuple());

    auto myparticle = spc.findMolecules("myparticle");
    Group mygroup(myparticle.front()->front(), myparticle.back()->back());
    mygroup.setMolSize(1);
    // In MPI version add:
    // slump.setDiscard(mpi.rank()+1);
    penalty->load("pf_");
    // Markov moves
    Move::Propagator<Tspace> mv(mcp, pot, spc);

    Table3D<double, double> histo(0.1, 0.1, Table3D<double, double>::HISTOGRAM);

    while ( loop[0] )
    {  // Markov chain
        while ( loop[1] )
        {
            mv.move();  // translate particle
            ++histo(spc.p[mygroup.front()].x(), spc.p[mygroup.front()].y());
        }
    }
    auto it_min = histo.min();
    histo.save("histo" + std::to_string(histo.getMap().size()), 1. / it_min->second);
    penalty->save("pf_", 1, {1.3, 1.3, -2, -2});
    penalty->saveRow("row_", {1.6}, 1, {1.3, 1.3, -2, -2});

    cout << loop.info() + mv.info() + penalty->info();

    // perform unit
    if ( penalty->update(true) != 0 )
    { // this ensures that the test is run only in the 2nd simulation
        UnitTest test(mcp);
        mv.test(test);
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
  penalty function method to overcome the energy barriers. 
  Below is the two-dimensional distribution functions sampled with and without
  the penalty function method, demonstrating how the penalty function 
  (at the bottom right) yields flat distribution functions in rough energy landscapes.

  ![Particle distribution functions and penalty function in an 2D oscillating field](penalty.png)

  ~~~
  $ cmake . -DENABLE_MPI=OFF
  $ make
  $ cd src/examples
  $ python ./penalty.py
  ~~~

  A gnuplot script to generate the plots of the probability distributions and 
  penalty function is provided:

  ~~~
  $ gnuplot penalty.gnu
  ~~~

  The penalty function routine in Faunus is also implemented in
  MPI using a master-slave scheme. Each system has its own rank and random seed.
  Samplings of the configurational space from all processes are merged by 
  periodically averaging over the two-dimensional distribution functions.
  The lines of code that need to be modified to run the program in parallel
  are indicated in the comments in penalty.cpp.
  
  The parameters to be set in the input file are the following:

  ~~~
  f0              0.5       # initial increment to the penalty function
  scale           0.2       # factor by which f0 is scaled
  update          1e4       # number of MC sweeps before scaling f0
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
