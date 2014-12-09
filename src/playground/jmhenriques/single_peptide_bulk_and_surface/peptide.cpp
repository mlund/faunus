#include <faunus/faunus.h>
/**
 * @brief  Single polymer/peptide/protein in bulk OR at a charged (and/or sticky/hydrophobic) 
 *         surface. Coarse-grained constant-pH, NVT, MC simulation with implicit solvent and 
 *         salt.
 *        
 * @author Joao Henriques
 * @date   2014/12/04
 */
using namespace Faunus;
//
#ifdef SURFACE
typedef Space<Geometry::Cuboidslit> Tspace;
typedef Potential::GouyChapman<>    Textpot1;
typedef Potential::StickyWall<>     Textpot2;
#else
typedef Space<Geometry::Cuboid>     Tspace;
#endif
//
typedef Potential::CombinedPairPotential<Potential::DebyeHuckel, 
                                         Potential::LennardJonesR12> 
                                         Tpairpot1;
//
typedef Potential::CombinedPairPotential<Tpairpot1,
                                         Potential::SquareWellHydrophobic>
                                         Tpairpot2;
//
int main() {
  cout << textio::splash();
  ::atom.includefile("param.json");
  InputMap mcp("simulation.input");
  EnergyDrift sys;
  Tspace spc(mcp);
  //  
#ifdef SURFACE
  auto pot = Energy::Nonbonded<Tspace, Tpairpot2>(mcp) 
           + Energy::ExternalPotential<Tspace, Textpot1>(mcp)
           + Energy::ExternalPotential<Tspace, Textpot2>(mcp)
           + Energy::MassCenterRestrain<Tspace>(mcp)           
           + Energy::EquilibriumEnergy<Tspace>(mcp)        
           + Energy::Bonded<Tspace>();
  pot.first.first.first.first.second.expot.setSurfPositionZ(&spc.geo.len_half.z());
  pot.first.first.first.second.expot.setSurfPositionZ(&spc.geo.len_half.z());
#else
  auto pot = Energy::Nonbonded<Tspace, Tpairpot2>(mcp) 
           + Energy::MassCenterRestrain<Tspace>(mcp)
           + Energy::EquilibriumEnergy<Tspace>(mcp)
           + Energy::Bonded<Tspace>();
#endif 
  auto bonded    = &pot.second;
  auto restraint = &pot.first.first.second;
  //  
  string file = mcp.get<string>("molecule", "");  
  double req  = mcp.get<double>("harmonic_eqdist", 0);
  double k    = mcp.get<double>("harmonic_forceconst", 0);
  //
  // TO DO: modify things to use the new addFastaSequence() 
  //
  Tspace::ParticleVector v;
  FormatAAM::load(file, v);
  Geometry::FindSpace().find(spc.geo, spc.p, v);
  Group pol = spc.insert(v);
  pol.name  = "peptide";
  spc.enroll(pol);
  for (int i = pol.front(); i < pol.back(); i++)
    bonded->add(i, i + 1, Potential::Harmonic(k, req)
                        + Potential::LennardJonesR12(mcp, "r12repex"));
  //
  restraint->add(pol);
  //
  double midpoint;
  double min = mcp.get<double>("bin_min", 0);
  double max = mcp.get<double>("bin_max", pc::infty);
  if (max == pc::infty)
    midpoint = 0;
  else
    midpoint = (min + max - spc.geo.len.z()) / 2;
  Point pt;
  pt.x() = 0;
  pt.y() = 0;
  pt.z() = midpoint;
  Point mc = Geometry::massCenter(spc.geo, spc.p, pol);
  pol.translate(spc, - mc + pt);
  pol.accept(spc);
  //
  FormatGRO gro;
  gro.len=spc.geo.len.x();
  //
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::TranslateRotate<Tspace> gmv(mcp, pot, spc);
  Move::CrankShaft<Tspace> crank(mcp, pot, spc);
  Move::Pivot<Tspace> pivot(mcp, pot, spc);
  Move::SwapMove<Tspace> tit(mcp, pot, spc, pot.first.second);
  //
  Analysis::PolymerShape shape;
  Average<double> avrg2;
  Analysis::ChargeMultipole mp;
#ifdef SURFACE
  Analysis::LineDistribution<> surfmcdist;
  int scnt = 0;
  double nmax;              // avoid 
  if (max != pc::infty)     // problems 
    nmax = max;             // with 
  else                      // infinite
    nmax = spc.geo.len.z(); // loops
  std::map<int, Analysis::LineDistribution<> > surfresdist;
  Table2D<double, Average<double> > netqtable;
  Table2D<double, Average<double> > rg2table;
#endif
  //
  spc.load("simulation.state");
  sys.init(Energy::systemEnergy(spc, pot, spc.p));
  //
  cout << "\nNumber of hydrophobic sites = " << ::numHydrophobic(spc.p, pol) << endl; // NEW
  //
  cout << atom.info()
       << pol.info()        
       << spc.info() 
       << pot.info() 
       << sys.info()
       << textio::header("The right man in the wrong place can make all the difference in the world.\n"
                         "  So, wake up, Mister Freeman. Wake up and smell the ashes.");
  //
  std::ofstream f1("rg_step.dat");
#ifdef SURFACE
  std::ofstream f2("surf_res_dist.dat");
#endif  
  //
  MCLoop loop(mcp);
  while (loop[0]) {
    while (loop[1]) {
      int i = slp_global.rand() % 5;
      switch (i) {
      case 0:
        mv.setGroup(pol);
        sys += mv.move();
        break;
      case 1:
        gmv.setGroup(pol);
        sys += gmv.move();
        break;
      case 2:
        crank.setGroup(pol);
        sys += crank.move();
        break;
      case 3:
        pivot.setGroup(pol);
        sys += pivot.move();
        break;
      case 4:
        sys += tit.move();
        mp.sample(pol,spc);
	break;
      }
      //   
      double rnd = slp_global();
      if (rnd < 0.05) {
        shape.sample(pol, spc);
        Point r2 = shape.vectorgyrationRadiusSquared(pol, spc);
        double rg2 = r2.x()+r2.y()+r2.z();
        avrg2 += rg2;
      }
      //
#ifdef SURFACE
      if (rnd < 0.05) {
        Point mc = Geometry::massCenter(spc.geo, spc.p, pol);                // mass center coord
        double dist = pot.first.first.first.first.second.expot.surfDist(mc); // mc dist to surf
        surfmcdist(dist)++;                                                  // mass center prob distr along the z axis
        netqtable(dist) += ::netCharge(spc.p, pol);                          // NEW
        Point rg2 = shape.vectorgyrationRadiusSquared(pol, spc);             // Rg2 coord
        rg2table(dist) += rg2.x() + rg2.y() + rg2.z();                       // Rg2 vs. dist to surf
        for (int i = pol.front(); i <= pol.back(); i++) {
          // res dist to surf
          double resdist = spc.geo.len.z()/2 - spc.p[i].z();
          // res prob distr along the z axis
          surfresdist[i](resdist)++;
        }
        scnt += 1;
      }
#endif
      //
    } // End of micro loop
    //
    sys.checkDrift(Energy::systemEnergy(spc, pot, spc.p));
    cout << loop.timing();
    //
    f1 << loop.innerCount() << " " << sqrt(avrg2.avg()) << "\n";
    //    
#ifdef SURFACE  
    netqtable.save("netq_dist.dat");
    rg2table.save("rg2_dist.dat");
#endif
    //
    gro.save("simulation.gro", spc.p, "append");
    spc.save("simulation.state"); // NEW
  //
  } // End of macro loop
  //
  FormatPQR::save("simulation.pqr", spc.p, spc.geo.len);
  //
#ifdef SURFACE
  surfmcdist.save("surf_mc_dist.dat");
  for (double d = spc.geo.len.z() - nmax; d <= spc.geo.len.z() - min; d += 0.25) {
    for (int i = pol.front(); i <= pol.back(); i++)
      f2 << d << " " << i + 1 << " " << -log(double(surfresdist[i](d))/double(scnt)) << endl; // Note: This will write 'inf' if probability is 0
    f2 << endl;
  }
#endif
  // 
  f1.close();
#ifdef SURFACE
  f2.close();
#endif
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
