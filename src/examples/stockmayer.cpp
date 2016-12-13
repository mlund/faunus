#include <faunus/faunus.h>
#include <faunus/multipole.h>

using namespace Faunus;
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid, DipoleParticle> Tspace;
typedef CombinedPairPotential<LennardJonesLB, DipoleDipole> Tpair;

int main()
{
    InputMap in("stockmayer.json");                // open parameter file for user input
    Tspace spc(in);                                // sim.space, particles etc.

    Energy::NonbondedVector<Tspace, Tpair> pot(in); // non-bonded only

    spc.load("state");

#ifdef __POLARIZE
    Move::Propagator<Tspace,true> mv(in,pot,spc);
#else
    Move::Propagator<Tspace, false> mv(in, pot, spc);
#endif

    Analysis::DipoleAnalysis dian(spc, in);
    DipoleWRL sdp;

    MCLoop loop(in);                               // class for mc loop counting
    while ( loop[0] )
    {                            // Markov chain
        while ( loop[1] )
        {
            mv.move();
            if ( slump() > 0.5 )
                dian.sampleMuCorrelationAndKirkwood(spc);
            dian.sampleDP(spc);
        }
        cout << loop.timing() << std::flush;
    }

    UnitTest test(in);
    mv.test(test);
    sdp.saveDipoleWRL("stockmayer.wrl", spc, Group(0, spc.p.size() - 1));
    FormatPQR().save("confout.pqr", spc.p);
    dian.save();
    spc.save("state");

    std::cout << spc.info() + pot.info() + mv.info() + test.info() + dian.info();

    return test.numFailed();
}
/**  @page example_stockmayer Example: Stockmayer potential

  This will simulate a Stockmayer potential in a cubic box.

  Run this example from the `examples` directory:

  ~~~~~~~~~~~~~~~~~~~
  $ make
  $ cd src/examples
  $ ./stockmayer.run
  ~~~~~~~~~~~~~~~~~~~

  stockmayer.cpp
  ==============

  @includelineno examples/stockmayer.cpp

*/
