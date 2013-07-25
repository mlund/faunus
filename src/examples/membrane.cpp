#include <faunus/faunus.h>

using namespace Faunus;

/**
 * @brief Setup interactions for coarse grained membrane
 *
 * More information: doi:10/chqzjk
 */
template<class Tbonded, class Tnonbonded, class Tlipid, class Tinput>
void MakeDesernoMembrane(const Tlipid &lipid, Tbonded &bond, Tnonbonded &nb, Tinput &in) {

  using namespace Potential;

  // non-bonded interactions
  auto hid=atom["HD"].id;
  auto tid=atom["TL"].id;
  double sigma   = in("lipid_sigma", 10);  // angstrom
  double epsilon = in("lipid_epsilon", 1); // kT

  auto headhead = WeeksChandlerAndersen(in) + DebyeHuckel(in);
  auto tailtail = WeeksChandlerAndersen(in) + CosAttract(in);
  auto headtail = WeeksChandlerAndersen(in) + ChargeNonpolar(in);

  // initialize WCA for head-head
  headhead.first.customSigma(hid, hid, 0.95*sigma);            
  headhead.first.customSigma(hid, tid, 0.95*sigma);            
  headhead.first.customSigma(tid, tid, sigma);                 
  headhead.first.customEpsilon(hid, hid, epsilon);             
  headhead.first.customEpsilon(hid, tid, epsilon);             
  headhead.first.customEpsilon(tid, tid, epsilon);   

  tailtail.first=headhead.first;      // copy initialized WCA to tail-tail and
  headtail.first=headhead.first;      // head-tail pair potential

  nb.pairpot.add(hid, hid, headhead); // Add to main pair potential
  nb.pairpot.add(tid, tid, tailtail);
  nb.pairpot.add(hid, tid, headtail);

  // bonded interactions
  double headtail_k=0.5*10*epsilon/(sigma*sigma);
  double headtail_req=4*sigma;
  double fene_k=30*epsilon/(sigma*sigma);
  double fene_rmax=1.5*sigma;

  assert(lipid.size() % 3 == 0);

  for (int i=lipid.front(); i<lipid.back(); i+=3) {
    bond.add(i,  i+1, FENE(fene_k,fene_rmax) );
    bond.add(i+1,i+2, FENE(fene_k,fene_rmax) );
    bond.add(i,  i+2, Harmonic(headtail_k,headtail_req) );
  }
}

typedef Space<Geometry::Cuboid> Tspace;
typedef Potential::PotentialMap<Potential::DebyeHuckelLJ> Tpairpot;

int main() {

  cout << textio::splash();      // show faunus banner and credits
  InputMap mcp("membrane.input");//read input file

  FormatXTC xtc(1000);
  EnergyDrift sys;               // class for tracking system energy drifts

  // Energy functions and space
  auto pot = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::Bonded<Tspace>()
    + Energy::ExternalPressure<Tspace>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);
  auto nonbonded = &pot.first.first.first;
  auto bonded = &pot.first.first.second;
  auto eqenergy = &pot.second;
  Tspace spc(mcp);

  // Load and add polymer to Space
  string file = mcp.get<string>("lipid_file","");
  int Nlipid=mcp("lipid_N",1);
  vector<Group> lipids(Nlipid);
  for (int i=0; i<Nlipid; i++) {
    Tspace::ParticleVector v;                   // temporary, empty particle vector
    FormatAAM::load(file,v);                    // load AAM structure into v
    Geometry::FindSpace().find(spc.geo,spc.p,v);// find empty spot in particle vector
    Group pol = spc.insert(v);                  // Insert into Space
    if (slp_global()>0.5)
      std::swap( spc.p[pol.front()].z(), spc.p[pol.back()].z() );
    if (slp_global()<mcp("lipid_chargefraction", 0.0))
      spc.p[ pol.front() ].charge = -1;
    pol.name="lipid";
    lipids[i]=pol;
    spc.enroll(lipids[i]);
  }

  // Set up bonded and non-bonded interactions
  Group allLipids(lipids.front().front(), lipids.back().back());
  allLipids.setMolSize(3);
  MakeDesernoMembrane(allLipids, *bonded, *nonbonded, mcp);

  // Place all lipids in xy plane (z=0);
  for (auto &g : lipids) {
    double dz=spc.p[ g.back() ].z();
    for (auto i : g) {
      spc.p[i].z() -= dz;
      spc.geo.boundary(spc.p[i]);
      g.setMassCenter(spc);
    }
  }
  spc.trial=spc.p;   // sync. particle trial vector
  spc.load("state"); // load old config. from disk (if any)

  // Markov moves and analysis
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::Pivot<Tspace> piv(mcp,pot,spc);
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::SwapMove<Tspace> swap(mcp,pot,spc,*eqenergy);
  Analysis::BilayerStructure lipidstruct;

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store total energy

  cout << atom.info() + spc.info() + pot.info();

  MCLoop loop(mcp);            // class for handling mc loops
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=slp_global.rand() % 3;
      int k=lipids.size(), j;
      Group g;
      switch (i) {
        case 0:
          mv.setGroup(allLipids);
          sys+=mv.move(allLipids.size()); // translate lipid monomers
          break;
        case 1:
          while (k-->0) {
            j=slp_global.rand() % (lipids.size());
            if (slp_global()>0.5) {
              gmv.setGroup(lipids[j]);          // tell what to move
              sys+=gmv.move();          // translate/rotate polymers
            } else {
              piv.setGroup(lipids[j]);          // tell what to move
              sys+=piv.move();          // translate/rotate polymers
            }
          }
          break;
        case 2:
          sys+=iso.move();
          break;
      }
      double ran = slp_global();
      if (ran>0.99) {
        xtc.setbox( spc.geo.len );
        xtc.save("traj.xtc", spc.p);
      }
      if (ran>0.90)
        lipidstruct.sample(spc.geo, spc.p, allLipids);

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?

    // save to disk
    FormatPQR::save("confout.pqr", spc.p);
    spc.save("state");

    cout << loop.timing();
  } // end of macro loop

  // perform unit tests
  UnitTest test(mcp);
  iso.test(test);
  mv.test(test);
  gmv.test(test);
  piv.test(test);
  sys.test(test);
  lipidstruct.test(test);

  // print information
  cout << loop.info() + mv.info() + gmv.info() + iso.info() + piv.info()
    + lipidstruct.info() + test.info() + sys.info();

  return test.numFailed();
}

/** @page example_membrane Example: Membrane Bilayer

  This will simulate a 3-bead coarse grained membrane according to
  Cooke and Deserno (doi:10/chqzjk). Each bead interacts with a
  Weeks-Chandler-Andersen potential, while tail-tail interactions
  have an additional long range attractive potential. There is preliminary
  support for charged head groups (effect on elastic properties is unknown).

  The following moves are included:
  - Lipid rotation, translation and pivot
  - Monomer translation
  - Iso-tension move

  Run this example from the `examples` directory:

  ~~~~~~~~~~~~~~~~~~~
  $ make
  $ cd src/examples
  $ ./membrane.run
  ~~~~~~~~~~~~~~~~~~~

  ![Bilayer formed by 3-bead CG lipid model](membrane-3bead.jpg)

  membrane.cpp
  ============

  \includelineno examples/membrane.cpp
  */

