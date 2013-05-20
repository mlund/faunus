#include <faunus/faunus.h>
#include <fstream>
using namespace std;
using namespace Faunus;
using namespace Potential;
typedef Geometry::Cuboid Tgeometry;
typedef Hertz Tpairpot;
//typedef CombinedPairPotential<YukawaGel,Hertz> Tpairpot;
int main(){
  int counter = 0;
  cout << textio::splash();
  
  
  InputMap mcp("HS.input");
  MCLoop loop(mcp);
  
  EnergyDrift sys;
  UnitTest test(mcp);
  
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp));
  Space spc(pot.getGeometry());
  
  
  Move::AtomicTranslation mv(mcp,pot,spc);
  
  
  GroupAtomic HS(spc, mcp);
  HS.name="HS";
  mv.setGroup(HS);
  Analysis::RadialDistribution<double, unsigned long long int> rdf_ab(4e2);
   spc.load("state");
  
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );
  
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
  while(loop.macroCnt()){
    while(loop.microCnt()){
      
      
      sys+=mv.move(HS.size());
      
      
      if (counter == 5) {
        counter = 0;
        particle::Tid a=atom["HS"].id, b=atom["HS"].id;
        rdf_ab.sample(spc,HS,a,b);
        
      }
      counter++;
    } // end of micro loop
    
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
    
  } // end of macro loop
  
  // save to disk
  FormatPQR().save("HS.pqr", spc.p); // final PQR snapshot for VMD etc.
  rdf_ab.save("HS.dat");         // g(r) - not normalized!
  spc.save("state");              // final simulation state
  
  // perform unit tests (irrelevant for the simulation)
  
  
  mv.test(test);
  sys.test(test);
  nonbonded->pairpot.test(test);
  
  // print information
  cout << loop.info() << sys.info() << mv.info()<< test.info();
  
  return test.numFailed();
}

