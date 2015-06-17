#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::PeriodicCylinder> Tspace;
typedef DebyeHuckelLJ Tpairpot;

int main(int argc, char** argv) {
  InputMap mcp("twobody.json");
  EnergyDrift sys;                     // class for tracking system energy drifts

  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp)
    + Energy::MassCenterConstrain<Tspace>(mcp, spc);

  Analysis::LineDistribution<> rdf(0.5);
  Analysis::ChargeMultipole mpol;
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  spc.load("state"); // load previous state, if any

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);
  while ( loop[0] ) {
    while ( loop[1] ) {
      sys += mv.move();

      rdf( spc.geo.dist( spc.groupList()[0]->cm, spc.groupList()[1]->cm ))++;
      for (auto i : spc.groupList())
        mpol.sample(*i, spc);

    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // detect energy drift

    cout << loop.timing();
    rdf.save("rdf.dat");
    FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
    spc.save("state");

  } // end of macro loop

  cout << loop.info() + sys.info() + mv.info() + mpol.info() << endl;

}

/**
 * @page example_twobody Example: Osmotic Second Virial Coefficient
 *
 * In this example we calculate the potential of mean force between
 * two rigid bodies with a fluctuating charge distribution. Salt and
 * solvent are implicitly treated using a Debye-Huckel potential and
 * the rigid bodies (here proteins) are kept on a line and allowed to
 * translate, rotate, and protons are kept in equilibrium with bulk
 * using particle swap moves.
 *
 * During simulation we simply sample the probability of observing
 * the two bodies at a specific separation (using `Analysis::LineDistribution`).
 * As the proteins are confined onto a line, each volume element along
 * `r` is constant and there's thus no need to normalize with a spherical
 * shell as done for freely moving particles.
 *
 * ![Two rigid molecules kept on a line inside a periodic cylinder](twobody.png)
 *
 * Run this example from the `examples` directory:
 *
 * ~~~~~~~~~~~~~~~~~~~
 * $ make example_twobody
 * $ cd src/examples
 * $ python twobody.py
 * ~~~~~~~~~~~~~~~~~~~
 *
 * For calculating the virial coefficient, @f$B_2@f$, the sampled radial distribution
 * function, g(r), needs to be integrated according to,
 *
 * @f[ B_2 = -2\pi \int_0^{\infty} [ g(r) -1 ] r^2 dr @f]
 *
 * where it should be noted that noise at large separations is strongly
 * amplified. To fix this, the tail of the g(r) can be replaced by
 * a Yukawa potential of the form (see [doi:10/xqw](http://dx.doi.org/10/xqw)),
 *
 * @f[ w(r)/k_BT = -\ln g(r) = \lambda_B s^2 z_iz_j / r \exp{(-\kappa r )}  @f]
 *
 * where @f$ s=\sinh{(a\kappa)}/a\kappa @f$ scales the charges due to spreading
 * them over a sphere with effective radius, `a`.
 * The fitting can be done using the supplied python script `twobody.virial.py` which
 * require information about protein charges, fitting range, and Debye length.
 * Edit and run with,
 *
 * ~~~~
 * $ python twobody.virial.py
 * ~~~~
 * (requires `scipy` and `matplotlib`)
 *
 * The output could look like below, and it is important to ensure
 * that the fitted parameters (Debye length or effective radius or both)
 * are in the vicinity of the simulated system.
 * ~~~~
 * Excecuting python script to analyse the virial coeffient...
 * Loaded g(r) file      =  rdf.dat
 * Saved w(r) file       =  wofr.dat
 * Particle charges      =  [7.1, 7.1]
 * Particle weights      =  [6520.0, 6520.0] g/mol
 * Fit range [rmin,rmax] =  40.0 90.0 A
 * Fitted Debye length   =  12.1634707424 A
 * Fitted ionic strength =  62.4643374077 mM
 * Fitted shift          =  8.31069879399 kT
 * Fitted radius         =  14.4531689553 A
 * Virial coefficient (cubic angstrom):
 * Hard sphere  [    0:   19] =  14365.4560073
 * Loaded data  [   19:   40] =  72438.2668732
 * Debye-Huckel [   40:  500] =  74969.9046117
 * TOTAL        [    0:  500] =  161773.627492 A^3
 *                            =  0.0229167635392 ml*mol/g^2
 * Reduced, B2/B2_HS          =  11.2612942749
 * ~~~~
 *
 * If `matplotlib` is installed, the python script will offer to
 * visualize the fit.
 *
 * ![Fitting of the PMF tail to a Yukawa potential](twobody-plot.png)
 *
 * twobody.cpp
 * ===========
 * @includelineno examples/twobody.cpp
 */
