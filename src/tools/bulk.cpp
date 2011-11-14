#include <faunus/faunus.h>
//#include "range.hpp"

using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;                // select simulation geometry
typedef Potential::CoulombSR<Tgeometry, Potential::Coulomb, Potential::HardSphere> Tpairpot;

namespace Faunus {
  namespace Analysis {
    class hist {
      protected:
        typedef double Tx;
        typedef unsigned int Ty;
        typedef std::map<Tx,Ty> Tmap;
        Ty count() {
          Ty cnt=0;
          for (auto &m : map)
            cnt+=m.second;
          return cnt;
        }
        Tx dx;
      private:
        Tmap map;
        Tx round(Tx x) {
          return floor(x);
        }
        virtual Ty normalize(Tx x) {
          return map[ round(x) ];
        }
      public:
        hist(Tx resolution=0.2) { dx=resolution; }
        Ty& operator[] (Tx x) {
          x=round(x);
          if (map.find(x)==map.end())
            map[x]=0;
          return map[x];
        }
        void save(string filename) {
          std::ofstream f(filename.c_str());
          if (f)
            for (auto &m : map)
              f << m.first << " " << normalize( m.second ) << endl;
        }
    };

    class rdf : public hist {
      private:
        Ty normalize(Tx x) {
          return 4./3.*pc::pi*( pow(x+0.5*dx,3) - pow(x-0.5*dx,3) );
        }
      public:
        rdf(Tx res=0.2) : hist(res) {}
    };
  }//namespace
}//namespace

int main() {
  cout << textio::splash();
  atom.includefile("atomlist.inp");    // load atom properties
  InputMap mcp("bulk.inp");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  EnergyDrift sys;                     // class for tracking system energy drifts

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot>(mcp) );
  Space spc( pot.getGeometry() );

  // Handle particles
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  spc.load("space.state", Space::RESIZE);

  /*
  salt.resize( 0 );
  cout << salt.beg << " " << salt.last << endl;
  cout << salt.front() << " " << salt.back() << endl;
  return 0;
  for (auto i=salt.begin(); i!=salt.end()-1; i++)
    for (auto j=i+1; j!=salt.end(); j++)
        cout << *i << " " << *j << endl;
  cout << endl;

  return 0;
*/

  Move::GrandCanonicalSalt gc(mcp,pot,spc,salt);
  Move::AtomicTranslation mv(mcp, pot, spc);  // Particle move class
  mv.setGroup(salt);

  // Widom particle insertion
  Analysis::rdf rdf(0.2);
  particle a;
  Analysis::Widom widom(spc, pot);
  a = atom["Cl"];
  widom.addGhost( spc);
  widom.addGhost( a );

#define UTOTAL \
  + pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt)\
  + pot.external()
  sys.init( UTOTAL );

  cout << atom.info() << spc.info() << pot.info() << mv.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move( salt.size() );
      sys+=gc.move( salt.size()/2 );
      if (slp_global.randOne()>0.9) {
        widom.sample();
        /*
        short ida=atom["Mg"].
        for (auto i=spc.p.begin(); i<spc.p.end()-1; i++)
          for (int j=i+1; j<spc.p.end(); j++)
            if (atom[].name==)
            */
      }
    }
    sys.checkDrift( UTOTAL );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");

  cout << sys.info() << loop.info() << mv.info() << gc.info() << widom.info();
}
