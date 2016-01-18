#include <faunus/faunus.h>

/*
 * This will simulate:
 * - an arbitrary number of rigid molecules (two different kinds)
 * - any number of atomic species
 * - equilibrium swap moves
 * - external pressure (NPT)
 * - arbitrary pair potential
 */

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid> Tspace;
typedef CombinedPairPotential<DebyeHuckelShift,CutShift<LennardJones> > Tpairpot;

int main(int argc, char** argv) {
  Faunus::MPI::MPIController mpi;

  InputMap mcp("manybody.json");
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts

  Tspace spc(mcp);
  auto pot
    = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::HydrophobicSASA<Tspace>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  auto eqenergy = &pot.second;

  Analysis::RadialDistribution<> rdf(0.5);
  Analysis::ChargeMultipole mpol;
  Scatter::DebyeFormula<Scatter::FormFactorUnity<>> debye(mcp);
  Scatter::DebyeFormula<Scatter::FormFactorUnity<>> debye2(mcp);

  spc.load("state"); // load previous state, if any

  Move::Propagator<Tspace> mv(mcp,pot,spc);
  mv.setMPI( &mpi );

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  std::ofstream cmfile, energyfile;

  if ( mpi.isMaster() ) {
    cmfile.open("cm.xyz");
    energyfile.open("energy.dat");
    cout << eqenergy->info();
    cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");
  }

  vector<Point> cm_vec; // vector of mass centers
  auto pol = spc.findMolecules("protein");

  MCLoop loop(mcp);
  while ( loop[0] ) {  // Markov chain
    while ( loop[1] ) {
      sys += mv.move();

      double rnd = slump();
      if ( rnd > 0.995 ) {
        // sample molecule-molecule g(r)
        for (auto i=pol.begin(); i!=pol.end()-1; i++)
          for (auto j=i+1; j!=pol.end(); j++)
            rdf( spc.geo.dist( (*i)->cm, (*j)->cm) )++;

        cm_vec.clear();
        for (auto &i : pol)
          cm_vec.push_back( i->cm );
        debye.sample(cm_vec,spc.geo.getVolume());

        if (mpi.isMaster())
          if (cmfile) {
            cmfile << cm_vec.size() << "\n" << "cm" << "\n";
            for (auto &m : cm_vec)
              cmfile << "H " << ((m+spc.geo.len_half)/10).transpose() << "\n";
          }
        if (energyfile)
          energyfile << loop.innerCount() << " " << sys.current()
            << " " << std::cbrt(spc.geo.getVolume()) << "\n";

        mpol.sample(pol,spc);
      }

      if ( rnd > 0.99 ) {
        xtc.setbox( spc.geo.len );
        xtc.save("traj.xtc", spc.p);
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // detect energy drift

    if (mpi.isMaster()) {
      cout << loop.timing();
      rdf.save("rdf_p2p.dat");
      FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
      spc.save("state");
      mcp.save("mdout.mdp");
      debye.save("debye.dat");
      debye2.save("debye.g2g.dat");
    }

  } // end of macro loop

  if (mpi.isMaster()) {
    cout << loop.info() + sys.info() + mv.info() + mpol.info() << endl;

    // save first molecule with average charges (as opposed to instantaneous)
    eqenergy->eq.copyAvgCharge(spc.p);
    spc.p.erase( spc.p.begin() + spc.groupList()[0]->size(), spc.p.end() );
    Geometry::cm2origo(spc.geo, spc.p);
    FormatPQR::save("avgcharge.pqr", spc.p, spc.geo.len);
    FormatAAM::save("avgcharge.aam", spc.p);
  }
}
