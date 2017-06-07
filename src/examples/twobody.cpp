#include <faunus/faunus.h>
#include <docopt.h>

using namespace Faunus;
using namespace Faunus::Potential;

#ifdef CYLINDER
    typedef Space<Geometry::Cylinder> Tspace;
#else
    typedef Space<Geometry::Cuboid> Tspace;
#endif

typedef CombinedPairPotential<CoulombGalore, LennardJonesLB> Tpairpot;

static const char USAGE[] =
R"(Faunus - A Framework for Molecular Simulation.

    http://github.com/mlund/faunus

    Usage:
      faunus [-q] [--rerun=TRAJ] [--state=FILE] (CONFIGFILE)
      faunus (-h | --help)
      faunus --version

    Options:
      CONFIGFILE                 JSON configuration file.
      --rerun=TRAJ               Rerun w. trajectory file (.xtc).
      -s FILE, --state=FILE      Initialize using state file.
      -q, --quiet                Less verbose output.
      -h --help                  Show this screen.
      --version                  Show version.
)";

int main( int argc, char **argv )
{
    auto args = docopt::docopt( USAGE,
            { argv + 1, argv + argc }, true, "Faunus 1.0.1");

    if ( args["--quiet"].asBool() )
        cout.setstate( std::ios_base::failbit ); // hold kæft!

    cout << textio::splash();

    Tmjson mcp = openjson( args["CONFIGFILE"].asString() );

    Tspace spc(mcp);

    auto pot = Energy::Nonbonded<Tspace, Tpairpot>(mcp)
        + Energy::EquilibriumEnergy<Tspace>(mcp)
        + Energy::MassCenterConstrain<Tspace>(mcp, spc);

    if ( args["--state"] ) // load old state
        spc.load( args["--state"].asString(), Tspace::RESIZE );

    Analysis::CombinedAnalysis analysis(mcp, pot, spc);
    Move::Propagator<Tspace> mv(mcp, pot, spc);

    cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

    MCLoop loop(mcp);
    while ( loop[0] )
    {
        while ( loop[1] )
        {
            mv.move();
            analysis.sample();
        } // end of micro loop
        cout << loop.timing();
    } // end of macro loop

    cout << loop.info() + mv.info() << endl;
}

/**
  @page example_twobody Example: Osmotic Second Virial Coefficient

  @todo Under serious construction!

  In this example we calculate the potential of mean force between
  two rigid bodies with a fluctuating charge distribution. Salt and
  solvent are implicitly treated using a Debye-Huckel potential and
  the rigid bodies (here proteins) are kept on a line and allowed to
  translate, rotate, and protons are kept in equilibrium with bulk
  using particle swap moves.

  During simulation we simply sample the probability of observing
  the two bodies at a specific separation.
  As the proteins are confined onto a line, each volume element along
  `r` is constant and there's thus no need to normalize with a spherical
  shell as done for freely moving particles.
  The minimum and maximum allowed distance between the two bodies can
  specified via the `Faunus::Energy::MassCenterConstrain` term in the
  Hamiltonian.

  ![Two rigid molecules kept on a line inside a cylinder](twobody.png)

  Run this example from the `examples` directory:
  The python script `twobody.py` generates the simulation input files
  and all changes to the simulation setup should be done from inside
  the script.
  This example can run with the following (the first, optional
  command checks out a specific commit with which this tutorial was created),

  ~~~~~~~~~~~~~~~~~~~
  $ (git checkout 3cbc360e7af0a8ef57e5abe687fed61517eee516)
  $ make example_twobody
  $ cd src/examples
  $ python twobody.py
  ~~~~~~~~~~~~~~~~~~~

  For calculating the virial coefficient, @f$B_2@f$, the angularly averaged
  radial distribution function, g(r), must be integrated according to,

  @f[ B_2 = -2\pi \int_0^{\infty} [ g(r) -1 ] r^2 dr @f]

  where it should be noted that noise at large separations is strongly
  amplified. To fix this, the tail of the g(r) can be fitted to
  a Yukawa potential of the form (see [doi:10/xqw](http://dx.doi.org/10/xqw)),

  @f[ w(r)/k_BT = -\ln g(r) = \lambda_B s^2 z_iz_j / r \exp{(-\kappa r )}  @f]

  where @f$ s=\sinh{(a\kappa)}/a\kappa @f$ scales the charges due to spreading
  them over a sphere with effective radius, `a`.
  The fitting can be done using the supplied python script `twobody.virial.py` which
  require information about protein charges, fitting range, and Debye length.
  Edit and run with,

  ~~~~
  $ python twobody.virial.py
  ~~~~
  (requires `scipy` and `matplotlib`)

  The output could look like below, and it is important to ensure
  that the fitted parameters (Debye length or effective radius or both)
  are in the vicinity of the simulated system.
  Note that in the
  limit @f$ a\kappa \rightarrow 0 @f$, i.e. when the Debye length
  is much bigger then the particle size, @f$ s \rightarrow 1 @f$
  and fitting via the effective radius, @f$ a @f$, is unfeasible
  and should be disabled in the script.

  ~~~~
  Excecuting python script to analyse the virial coeffient...
  Loaded g(r) file      =  rdf.dat
  Saved w(r) file       =  wofr.dat
  Particle charges      =  [7.1, 7.1]
Particle weights      =  [6520.0, 6520.0] g/mol
Fit range [rmin,rmax] =  40.0 90.0 A
Fitted Debye length   =  12.1634707424 A
Fitted ionic strength =  62.4643374077 mM
Fitted radius         =  14.4531689553 A
Virial coefficient (cubic angstrom):
    Hard sphere  [    0:   19] =  14365.4560073
    Loaded data  [   19:   40] =  72438.2668732
    Debye-Huckel [   40:  500] =  74969.9046117
    TOTAL        [    0:  500] =  161773.627492 A^3
    =  0.0229167635392 ml*mol/g^2
    Reduced, B2/B2_HS          =  11.2612942749
    ~~~~

    If `matplotlib` is installed, the python script will offer to
    visualize the fit.

![Fitting of the PMF tail to a Yukawa potential](twobody-plot.png)

    twobody.cpp
    ===========
    This is the underlying C++ program where one should note that the
    pair potential is hard coded but can naturally be customised to
    any potential from namespace `Faunus::Potential`.

    @includelineno examples/twobody.cpp
    */
