#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace;
typedef CombinedPairPotential<LennardJonesLB,MultipoleWolf<true,true,true,false>> Tpairpot;
//typedef CombinedPairPotential<LennardJonesLB,Coulomb> Tpairpot1;
//typedef CombinedPairPotential<Tpairpot1,IonDipole> Tpairpot2;
//typedef CombinedPairPotential<Tpairpot2,DipoleDipole> Tpairpot;

//typedef CombinedPairPotential<LennardJonesLB,IonIonGaussianDamping> Tpairpot1;
//typedef CombinedPairPotential<Tpairpot1,IonDipoleGaussianDamping> Tpairpot2;
//typedef CombinedPairPotential<Tpairpot2,DipoleDipoleGaussianDamping> Tpairpot;

int main() {
  InputMap mcp("water2.input");//read input file
  slp_global.seed(mcp.get<int>("seed", -7));

  // Energy functions and space
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
  Tspace spc(mcp);

  // Load and add polymer to Space
  auto N    = mcp.get<int>("mol_N",1);
  auto file = mcp.get<string>("mol_file");
  vector<Group> water(N);
  Tspace::ParticleVector v;                   // temporary, empty particle vector
  FormatAAM::load(file,v);                    // load AAM structure into v
  for (auto &i : water) {
    Geometry::FindSpace f;
    f.allowMatterOverlap=true;
    f.find(spc.geo,spc.p,v);// find empty spot in particle vector
    i = spc.insert(v);                          // Insert into Space
    i.name="h2o";
    spc.enroll(i);
  }

  // Markov moves and analysis
  Move::PolarizeMove<TranslateRotate<Tspace> >  gmv(mcp,pot,spc);
  //Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  
  Analysis::RadialDistribution<> rdfOO(0.05);
  Analysis::RadialDistribution<> rdfOH(0.05);
  Analysis::RadialDistribution<> rdfHH(0.05);

  spc.load("state"); // load old config. from disk (if any)
  Analysis::DipoleAnalysis dian(spc,mcp);
  dian.load(mcp.get<string>("dipole_data_ext", ""));

  FormatPQR::save("initial.pqr", spc.p);

  EnergyDrift sys;   // class for tracking system energy drifts
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store total energy

  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);            // class for handling mc loops
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=slp_global.rand() % 2;
      int j,k=water.size();
      Group g;
      switch (i) {
        case 0:
          while (k-->0) {
            j=slp_global.rand() % (water.size());
            gmv.setGroup(water[j]);
            sys+=gmv.move();          // translate/rotate polymers
          }
          break;
        case 1:
          sys+=iso.move();
          break;
      }
      dian.sampleDP(spc);
      // sample oxygen-oxygen rdf
      if (slp_global()>0.5) {
        dian.sampleMuCorrelationAndKirkwood(spc);
        auto idO = atom["OW"].id;
        auto idH = atom["HW"].id;
        rdfOO.sample(spc,idO,idO);
        rdfOH.sample(spc,idO,idH);
        rdfHH.sample(spc,idH,idH);
      }
    } // end of micro loop
    //rdfOO.save("rdfOO.dat"+std::to_string(loop.count()));
    //rdfOH.save("rdfOH.dat"+std::to_string(loop.count()));
    //rdfHH.save("rdfHH.dat"+std::to_string(loop.count()));
    //dian.save(std::to_string(loop.count()));
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?

    cout << loop.timing();
  } // end of macro loop

  // save to disk
  spc.save("state");
  rdfOO.save("rdfOO.dat");
  rdfOH.save("rdfOH.dat");
  rdfHH.save("rdfHH.dat");
  spc.save("state");
  dian.save();
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);

  // print information
  cout << loop.info() + sys.info() + gmv.info() + iso.info() + dian.info();
}
