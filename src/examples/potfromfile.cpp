#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Move;
using namespace Faunus::Potential;


typedef Space<Geometry::Cuboid> Tspace;
typedef CombinedPairPotential<HardSphere,Potfromfile> Tpair;
typedef CombinedPairPotential<Tpair,DebyeHuckel> Tpair2;


int main()
{
    InputMap in("potfromfile.json");                // open parameter file for user input
    Tspace spc(in);                                // sim.space, particles etc.


    auto pot = Energy::Nonbonded<Tspace,Tpair2>(in);


    spc.load("state");


    Move::Propagator<Tspace, false> mv(in, pot, spc);

    Analysis::CombinedAnalysis analyzer(in,pot,spc);

    MCLoop loop(in);                               // class for mc loop counting
    while ( loop[0] )
    {                            // Markov chain
        while ( loop[1] )
        {
            mv.move();
	    analyzer.sample();
        }
        cout << loop.timing() << std::flush;
    }

    UnitTest test(in);
    mv.test(test);
    analyzer.test(test);
   

    std::cout << spc.info() + pot.info() + mv.info() + test.info() + analyzer.info();

    return test.numFailed();
}
/**  @page example_potfromfile Example: Preform a calculation using a Potential saven in an external file. 

  This will simulate a combined potential of Hard Shpere a Potential 
  From an external file and a Debye Hyckel potential in a cubic box.

  Run this example from the `examples` directory:

  ~~~~~~~~~~~~~~~~~~~
 
  ~~~~~~~~~~~~~~~~~~~

  stockmayer.cpp
  ==============

  @includelineno examples/stockmayer.cpp

 */
