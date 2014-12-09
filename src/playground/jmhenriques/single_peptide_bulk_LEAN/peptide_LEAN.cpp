#include <faunus/faunus.h>
/*
 * @brief  Lean version of the original peptide.cpp (bulk only!, i.e no surface) with extra comments 
 *         for teaching purposes.
 *
 * @author Joao Henriques
 * @date   2014/12/04
 */
using namespace Faunus;
//
typedef Space<Geometry::Cuboid> Tspace; // Define the simulation box (space) type (cuboid in this case).
//
// Sum three non-bonded pair potentials (Debye-Huckel electrostatics, Pauli repulsion [aka. Lennar-Jones R12]
// and Hydrophobic interactions [aka. Square Well]). Note: It can only sum 2 at a time, therefore we sum twice.
//
typedef Potential::CombinedPairPotential<Potential::DebyeHuckel, 
                                         Potential::LennardJonesR12> 
                                         Tpairpot1;
//
typedef Potential::CombinedPairPotential<Tpairpot1,
                                         Potential::SquareWellHydrophobic>
                                         Tpairpot2;
//
// Main function, where everything happens.
//
int main() {
  cout << textio::splash();         // Print credits.
  ::atom.includefile("param.json"); // Load amino acid properties (name, radius, mass, charge, pKa, etc).
  InputMap mcp("simulation.input"); // Load simulation input file (box size, how many cycles, temperature, etc).
  EnergyDrift sys;
  Tspace spc(mcp);                  // Name space and get its properties from the input file (mcp), loaded before.
  //
  // Define the system Hamiltonian.
  //
  auto pot = Energy::Nonbonded<Tspace, Tpairpot2>(mcp) // Non-bonded terms (using Tpairpot2 defined above).  
           + Energy::EquilibriumEnergy<Tspace>(mcp)    // Titration.
           + Energy::Bonded<Tspace>();                 // Bonded terms (not defined yet).
  auto bonded    = &pot.second;                        // Give a name to the bonded terms.
  //
  string file = mcp.get<string>("molecule", "");           // Find the file with the molecule geometry (inside the sim. input, mcp).  
  double req  = mcp.get<double>("harmonic_eqdist", 0);     // Find the equilibrium distance between bonds.
  double k    = mcp.get<double>("harmonic_forceconst", 0); // Find the bond force constant.
  //
  Tspace::ParticleVector v;                      // Create a container for our protein.
  FormatAAM::load(file, v);                      // Load protein structure (found in the 'file' we created above) into 'v'.
  Geometry::FindSpace().find(spc.geo, spc.p, v); // Finds space to put it inside the simulation space.
  Group pol = spc.insert(v);                     // Insert it.
  pol.name  = "peptide";                         // Name it.
  spc.enroll(pol);                               // Enroll it.
  //
  // Here we are going to add the harmonic bonds in between 'i' and 'i+1' (we also remove the Pauli repulsion in between them).
  //
  for (int i = pol.front(); i < pol.back(); i++) 
    bonded->add(i, i + 1, Potential::Harmonic(k, req)
                        + Potential::LennardJonesR12(mcp, "r12repex"));
  //
  FormatGRO gro;           // Create a GRO file type to save the simulation snapshots.
  gro.len=spc.geo.len.x(); // Box size (for visualization purposes).
  //
  // Types of moves our MC program will attempt.
  //
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);           // Pick a single particle (amino acid in this case) and displace it.
  Move::TranslateRotate<Tspace> gmv(mcp, pot, spc);            // Translate and rotate the protein (will be necessary further ahead in the project).
  Move::CrankShaft<Tspace> crank(mcp, pot, spc);               // Crankshaft move.
  Move::Pivot<Tspace> pivot(mcp, pot, spc);                    // Pivot move.
  Move::SwapMove<Tspace> tit(mcp, pot, spc, pot.first.second); // Titration move, ie. attempt to protonate/deprotonate amino acids.
  //
  Analysis::PolymerShape shape; // Create an analysis routine to sample the protein shape. 
  Analysis::ChargeMultipole mp; // Create an analysis routine to sample the protein charge/capacitance. 
  //
  spc.load("simulation.state");                    // Load checkpoint file if it exists (useful to do sim. restarts).
  sys.init(Energy::systemEnergy(spc, pot, spc.p)); // Compute and store initial system energy.
  //
  // Print information about many different things.
  //
  cout << atom.info()
       << "\nNumber of hydrophobic sites = " << ::numHydrophobic(spc.p, pol) << endl
       << pol.info()        
       << spc.info() 
       << pot.info() 
       << sys.info()
       << textio::header("- - - - - S T A R T I N G - - - T H E - - - S I M U L A T I O N - - - - -");
  //
  // Simulation loop.
  //
  MCLoop loop(mcp);
  //
  // Start the macro loop.
  //
  while (loop[0]) {
    //
    // Start the micro loop.
    //
    while (loop[1]) {
      int i = slp_global.rand() % 5; // Pick a random number between 0 and 4.
      switch (i) {
      // Particle translation.
      case 0:
        mv.setGroup(pol);
        sys += mv.move();
        break;
      // Group translation/rotation.
      case 1:
        gmv.setGroup(pol);
        sys += gmv.move();
        break;
      // Crankshaft move.
      case 2:
        crank.setGroup(pol);
        sys += crank.move();
        break;
      // Pivot move. 
      case 3:
        pivot.setGroup(pol);
        sys += pivot.move();
        break;
      // Titration move.
      case 4:
        sys += tit.move();
        mp.sample(pol,spc);
      }
      //
      // Sample the shape analysis.
      //
      double rnd = slp_global();
      if (rnd < 0.05) {
        shape.sample(pol, spc);
      }
    }
    //
    // End of micro loop.
    //
    sys.checkDrift(Energy::systemEnergy(spc, pot, spc.p)); // Check energy drift.
    cout << loop.timing();                                 // Print estimated time to completion.
    //
    gro.save("simulation.gro", spc.p, "append"); // Save snapshot to the GRO file.
    spc.save("simulation.state");                // Save the current simulation state/checkpoint.
  }
  //
  // End of macro loop
  //
  FormatPQR::save("simulation.pqr", spc.p, spc.geo.len); // Save a PQR file with last conformation, masses, charges, etc.
  //
  // Print simulation/analysis statistics.
  //
  cout << sys.info() 
       << mv.info() 
       << gmv.info() 
       << pivot.info() 
       << crank.info() 
       << tit.info()
       << mp.info()
       << shape.info() 
       << loop.info();
}
