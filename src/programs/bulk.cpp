#include <faunus/faunus.h>

using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;                // select simulation geometry
typedef Potential::CoulombSR<Tgeometry, Potential::Coulomb, Potential::LennardJones> Tpairpot;

//typedef Potential::CoulombSR<Tgeometry, Potential::Coulomb, Potential::HardSphere> Tpairpot;

template<class T>
class distributions {
  private:
    typedef xytable<double,T> Ttable;
    typedef std::map<string,Ttable> Tmap;
    double xmin, xmax, dx;
    Tmap map;
  public:
    distributions(double min, double max, double delta) {
    }
    void add(string name, double val) {
      map[name]+=val;
    }
};

/*
 * \brief Auxillary class for tracking atomic species
 * \author Mikael Lund
 * \date Malmo 2011
 *
 * This class keeps track of individual particle positions based on particle type (id).
 * It contains functions to insert and erase particles while automatically moving particles
 * above the deletion or insertion point in sync.
 *
 * Example:
 * \code
 * AtomTracker track(myspace);
 * track.insert( myparticle );
 * ...
 * int i=track[ myparticle.id ].random();
 * \endcode
 */
class AtomTracker {
  private:
    typedef short Tid;  // particle id type
    typedef int Tindex; // particle index type
    Space* spc;
    class data {
      public:
        Tid id;
        vector<Tindex> index;
        Tindex random();   //!< Pick random particle index
    };
    std::map<Tid,data> map; 
  public:
    AtomTracker(Space&);
    Tid randomAtomType() const;   //!< Select a random atomtype from the list
    bool insert(const particle&); //!< Insert particle into Space and track position
    bool erase(Tindex);           //!< Delete particle from Space at specific particle index
    data& operator[] (Tid);       //!< Access operator to atomtype data
};

AtomTracker::Tid AtomTracker::randomAtomType() const {
  assert(!map.empty() && "No atom types have been added yet");
  vector<Tid> vid;
  for (auto &m : map)
    vid.push_back(m.second.id);
  std::random_shuffle(vid.begin(), vid.end());
  return vid.front();
}

AtomTracker::AtomTracker(Space &s) { spc=&s; }

AtomTracker::Tindex AtomTracker::data::random() {
  std::random_shuffle( index.begin(), index.end() );
  return index.back();
}

AtomTracker::data& AtomTracker::operator[](AtomTracker::Tid id) { return map[id]; }

/*!
 * \warning Assumes Space inserts in the back of vector
 */
bool AtomTracker::insert(const particle &a) {
  spc->insert(a);
  Tindex index=(Tindex)spc->p.size()-1; // !!!!
  map[a.id].index.push_back(index);
  for (auto &m : map)
    for (auto &i : m.second.index)
      if (i>index) i++;
   return true;
}

bool AtomTracker::erase(AtomTracker::Tindex index) {
  spc->remove(index);
  bool deleted=false;
  for (auto &m : map) {
    auto f=std::find(m.second.index.begin(), m.second.index.end(), index);
    if (f!=m.second.index.end()) {
      assert( spc->p[index].id==m.second.id && "Atom id mismatch");
      m.second.index.erase(f);
      deleted=true;
      break;
    }
  }
  if (deleted)
    for (auto &m : map)
      for (auto &i : m.second.index)
        if (i>index) i--;
  assert(deleted && "Could not delete specified index");
  return deleted;
}

int main() {
  cout << textio::splash();
  distributions<double> dst(0,100,0.5);
  atom.includefile("atomlist.inp");    // load atom properties
  InputMap mcp("bulk.inp");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                           // PQR structure file I/O
  EnergyDrift sys;                     // class for tracking system energy drifts

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  // Handle particles
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  Move::ParticleTranslation mv(mcp, pot, spc);  // Particle move class
  mv.setGroup(salt);

  GroupMolecular test;
  test.beg=spc.p.size();
  spc.insert("Na", 2);
  test.end=spc.p.size()-1;
  test.name="Test";
  spc.p[test.beg].x=0;
  spc.p[test.beg].y=0;
  spc.p[test.beg].z=0;
  spc.p[test.end].x=2;
  spc.p[test.end].y=2;
  spc.p[test.end].z=2;
  spc.trial=spc.p;
  spc.enroll(test);

  spc.load("space.state");
  salt.setMassCenter(spc);
  test.setMassCenter(spc);

  //bonded->bonds.add(0,1, Potential::Harmonic(0.2, 10.0));
  //bonded->bonds.add(1,2, Potential::Harmonic(0.3,  5.0));

  // Particle titration
  Move::SwapMove tit(mcp,pot,spc);
  Move::Isobaric iso(mcp,pot,spc);
  Move::RotateGroup gmv(mcp,pot,spc);

  // Widom particle insertion
  Analysis::Widom widom(spc, pot);
  widom.addGhost(spc);

  //FastaSequence fasta;
  //Group protein1 = fasta.insert( "AAK", spc, bonded->bonds );

#define UTOTAL \
  pot.g_external(spc.p, test)\
  + pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt)\
  + pot.g2g(spc.p,salt,test) + pot.external()

  //#define UTOTAL pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt) + pot.external()

  sys.init( UTOTAL );

  cout << atom.info() << spc.info() << pot.info() << mv.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move();
      if (slp_global.randOne()>0.9)
        widom.sample();
      sys+=tit.move();
      gmv.setGroup(test);
      sys+=gmv.move();
      sys+=iso.move();
    }
    sys.checkDrift( UTOTAL );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");

  cout << mv.info() << sys.info() << loop.info() << widom.info()
    << tit.info() << gmv.info() << iso.info();
}
