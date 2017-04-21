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

  Tspace spc(mcp);
  auto pot
    = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::HydrophobicSASA<Tspace>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  auto eqenergy = &pot.second;

  spc.load("state"); // load previous state, if any

  Analysis::CombinedAnalysis analysis(mcp,pot,spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);
  mv.setMPI( &mpi );

  std::ofstream cmfile;

  if ( mpi.isMaster() ) {
    cmfile.open("cm.xyz");
    cout << eqenergy->info();
    cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");
  }

  vector<Point> cm_vec; // vector of mass centers
  auto pol = spc.findMolecules("protein");

  MCLoop loop(mcp);
  while ( loop[0] ) {  // Markov chain
    while ( loop[1] ) {
      mv.move();

      if (mpi.isMaster()) {
          analysis.sample();

          if (cmfile) {
              cm_vec.clear();
              for (auto &i : pol)
                  cm_vec.push_back( i->cm );
              cmfile << cm_vec.size() << "\n" << "cm" << "\n";
              for (auto &m : cm_vec)
                  cmfile << "H " << ((m+spc.geo.len_half)/10).transpose() << "\n";
          }
      }

    } // end of micro loop

    if (mpi.isMaster())
      cout << loop.timing();

  } // end of macro loop

  if (mpi.isMaster()) {
    cout << loop.info() + mv.info() + analysis.info() << endl;

    // save first molecule with average charges (as opposed to instantaneous)
    eqenergy->eq.copyAvgCharge(spc.p);
    spc.p.erase( spc.p.begin() + spc.groupList()[0]->size(), spc.p.end() );
    Geometry::cm2origo(spc.geo, spc.p);
    FormatPQR::save("avgcharge.pqr", spc.p, spc.geo.len);
    FormatAAM::save("avgcharge.aam", spc.p);
  }
}
