#include <faunus/faunus.h>
#include <faunus/multipole.h>
using namespace std;
using namespace Faunus;                          // use Faunus namespace
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Geometry::Cuboid Tgeo;                   // select simulation geometry and pair potential
typedef CombinedPairPotential<Hertz,DipoleDipoleRF> Tpair;
//typedef Potential::DipoleDipoleRF Tpair;
int main() {
  int counter = 0;
  cout << textio::splash();

  ::atom.includefile("stockmayer.json");         // load atom properties
  InputMap in("stockmayer.input");               // open parameter file for user input
  Energy::NonbondedVector<Tpair,Tgeo> pot(in);   // create Hamiltonian, non-bonded only
  EnergyDrift sys;                               // class for tracking system energy drifts
  Space spc( pot.getGeometry() );                // create simulation space, particles etc.
  GroupAtomic sol(spc, in);                      // group for particles
  sol.name = "sol";
  MCLoop loop(in);                               // class for mc loop counting
  Analysis::RadialDistribution<> rdf_ab(1e2);       // particle-particle g(r)

  Move::AtomicTranslation trans(in, pot, spc);   // particle move class
  Move::AtomicRotation rot(in, pot, spc);        // particle move class
  
  

  //PolarizeMove<AtomicTranslation> trans(in,pot,spc);
  //PolarizeMove<AtomicRotation> rot(in,pot,spc);
  trans.setGroup(sol);                                // tells move class to act on sol group
  rot.setGroup(sol);                                  // tells move class to act on sol group

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy
  cout << atom.info() << spc.info() << pot.info() << textio::header("MC begins");

  /*
    ofstream myfile;
    myfile.open ("potential.dat");
  
    

   particle a = spc.p[0];

   particle b = spc.p[1];

   a.x()=0;

   a.y()=0;

   a.z()=0;



   b.x()=0;

   b.y()=0;

   b.z()=0.0;

   

   // myfile << "a " << a << endl;

   //   myfile << "b " << b << endl;

   

   double dz =  50;

   for (int i = 0; i <  1420; i++) {

     b.z()+=dz;

     double p = pot.p2p(a, b);

     myfile << b.z() << " " << p << endl;

   }
   myfile.close();
   return 0;

  */

 
  while ( loop.macroCnt() ) {                         // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() > 0.5)
        sys+=trans.move( sol.size() );                     // translate
      else
        sys+=rot.move( sol.size() );                       // rotate

       if (counter == 5) {
        counter = 0;
        particle::Tid a=atom["sol"].id, b=atom["sol"].id;
        rdf_ab.sample(spc,sol,a,b);
        
      }
      counter++;
    }

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
  }

  FormatPQR().save("confout.pqr", spc.p);
  rdf_ab.save("Stock.dat");                               // save g(r) to disk
  std::cout << spc.info() + pot.info() + trans.info() + rot.info() + sys.info(); // final info
}
