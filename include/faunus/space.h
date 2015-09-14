#ifndef FAUNUS_SPACE_H
#define FAUNUS_SPACE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/slump.h>
#include <faunus/point.h>

#include <faunus/inputfile.h>
#include <faunus/species.h>
#include <faunus/physconst.h>
#include <faunus/group.h>
#include <faunus/point.h>
#include <faunus/space.h>
#include <faunus/textio.h>
#include <faunus/io.h>
#include <faunus/molecule.h>
#endif

namespace Faunus {

  /**
   * @brief Particle and molecule tracker
   *
   * This will track data (`T`) associated with an
   * id (`Tid`) and is typically used to track atoms and
   * molecules in grand canonical moves. Only one instance
   * of a data element can exist at any given time.
   *
   * Internally, the structure is a map of `Tid` key and
   * a vector of `T`.
   *
   * Example:
   *
   * ~~~~
   *   Tracker<Group*> molTrack;
   *   Tracker<int> atomTrack;
   *
   *   // build list of atoms and molecules
   *   for ( auto g : spc.groupList() )
   *     if ( g->isAtomic() )
   *       for ( auto i : *g )
   *         atomTrack.insert( spc.p[i].id, i );
   *     else
   *       molTrack.insert( g->molId, g );
   *
   *   // find three, random sodium atoms
   *   vector<int> result;
   *   auto id = name["Na"].id
   *   atomTrack.find( id, 3, result );
   *
   *   // erase one of the found atoms from tracker
   *   atomTrack.erase( id, result[0] );
   * ~~~~
   */
  template<class T, class Tid=int>
    class Tracker {
      private:
        std::map<Tid, std::vector<T> > _map;
        std::map<Tid, Average<double> > Navg; // average vector length. TO BE REMOVED
      public:
        /** @brief Number of elements of type id */
        size_t size(Tid id) const {
          auto i=_map.find(id);
          if (i!=_map.end())
            return i->second.size();
          return 0;
        }

        /** @brief Sum of all elements **/
        size_t size() const {
          size_t sum=0;
          for (auto &m : _map)
            sum+=m.second.size();
          return sum;
        }

        const std::vector<T>& operator[](Tid i) const { return _map.at(i); };
        std::vector<T>& operator[](Tid i) { return _map[i]; };

        /** @brief Clear map -- preserve averages */
        void clear() { _map.clear(); }

        /** @brief Update average number of particles. TO BE REMOVED. */
        void updateAvg() {
          for (auto &m : _map)
            Navg[m.first] += m.second.size();
        }

        /** @brief Get average number of particles */
        Average<double> getAvg(Tid id) { return Navg[id]; }

        /** @brief Update atom tracker. TO BE REMOVED. */
        template<class Tpvec>
          void update(const Tpvec &p) {
            _map.clear();
            for (size_t i=0; i<p.size(); i++)
              if ( atom[ p[i].id ].activity > 1e-6 )
                _map[ p[i].id ].push_back( i );
          }

        /**
         * @brief Find random data based on id
         * @param id Id for data
         * @param N Number of unique elements (data) to return
         * @param dst Destination vector -- new elements added to end
         * @return True if found; false otherwise
         */
        bool find(Tid id, size_t N, std::vector<T> &dst) {
          size_t Ninit=dst.size();
          auto i=_map.find(id);
          if (i!=_map.end()) // id found
            if (i->second.size()>=N) {  // enough elements?
              do {
                auto j= slump.element(i->second.begin(), i->second.end());
                if ( std::find(dst.begin(), dst.end(), *j)==dst.end() )
                  dst.push_back(*j);
              } while (dst.size()!=Ninit+N);
              return true;
            }
          return false;
        }

        /**
         * @brief Insert data (if not already present)
         * @param id Id of data
         * @param data Data to track
         */
        void insert(Tid id, T data) {
          auto i = _map.find(id);
          if (i!=_map.end()) {
            auto pos=std::find(i->second.begin(), i->second.end(), data);
            if (pos==i->second.end())
              i->second.push_back(data);
          }
          else
            _map[id].push_back(data);
        }

        /**
         * @brief Erase data
         * @param id Id of data
         * @param data Data to erase
         */
        bool erase(Tid id, T data) {
          auto i = _map.find(id);
          if (i!=_map.end()) {
            auto pos=std::find(i->second.begin(), i->second.end(), data);
            if (pos!=i->second.end()) {
              std::swap( *pos, i->second.back() );
              i->second.pop_back();
              return true;
            }
          }
          return false;
        }

        /**
         * @brief Erase data using search (slower than using id)
         */
        bool erase(Tid id) {
          for (auto &m : _map) {
            auto it = std::find(m.second.begin(), m.second.end(), id);
            if (it!=m.second.end()) {
              m.second.erase(it);
              return true;
            }
          }
          return false;
        }
    };

  /**
   * @brief Placeholder for particles and groups
   *
   * Every simulation must have a `Space` instance as this contains
   * the particles and information about groups (particle ranges).
   * Space also contains the geometry as given by the template
   * parameter.
   *
   * @todo insert/erase functions are a mess.
   */
  template<typename Tgeometry, typename Tparticle=PointParticle>
    class Space {
      protected:
        std::ifstream fin;
      private:
        bool checkSanity();                    //!< Check group length and vector sync
        std::vector<Group*> g;                 //!< Pointers to ALL groups in the system

      public:
        typedef std::vector<Tparticle, Eigen::aligned_allocator<Tparticle> > p_vec;
        typedef p_vec ParticleVector;          //!< Particle vector type
        typedef Tparticle ParticleType;        //!< Particle type
        typedef Tgeometry GeometryType;        //!< Geometry type
        typedef MoleculeData<ParticleVector> MoleculeType;

        enum keys {OVERLAP_CHECK,NOOVERLAP_CHECK,RESIZE,NORESIZE};

        Tgeometry geo;                         //!< System geometry
        ParticleVector p;                      //!< Main particle vector
        ParticleVector trial;                  //!< Trial particle vector.
        MoleculeMap<ParticleVector> molecule;  //!< Map of molecules

        Tracker<int> atomTrack;                //!< Track atom index based on atom type
        Tracker<Group*> molTrack;              //!< Track groups pointers based on molecule type

        /**
         * @brief Constructor
         *
         * This will pass the json object to `MoleculeMap`
         * in order to load all molecule types. The molecule map
         * is searched for molecules with non-zero `Ninit` and
         * will insert accordingly.
         */
        Space( Tmjson &j ) : geo( j ) {
          pc::setT( j["system"]["temperature"] | 298.15 );
          atom.include( j );
          molecule.include( j );
          for (auto mol : molecule)
            while (mol.Ninit-- > 0)
              insert( mol.id, mol.getRandomConformation(geo, p) );
        }

        std::vector<Group*>& groupList() { return g; };   //!< Vector with pointers to all groups

        AtomMap& atomList() { return atom; } //!< Vector of atoms

        MoleculeMap<ParticleVector>& molList() { return molecule; } //!< Vector of molecules

        bool save(string);                  //!< Save container state to disk
        bool load(string, keys=NORESIZE);   //!< Load container state from disk

        /** @brief insert p_vec of MolID to end of p and trial */
        Group* insert(PropertyBase::Tid, const p_vec&); // inserts to trial and p

        /*
         * @brief insert p_vec of MolID
         * @param enlarge - sets whether to enlarge group of MolID, or to add new isMolecular()==true group
         */
        //Group* insert(const p_vec&, PropertyBase::Tid, bool enlarge=true);
        //Group insert(const p_vec&, int=-1);
        bool insert(const Tparticle&, int=-1); //!< Insert particle at pos n (old n will be pushed forward).
        //bool insert(string, int, keys=OVERLAP_CHECK);
        bool erase(int);             //!< Remove n'th particle and downshift/remove groups
        //bool eraseGroup(Group& group); ///< find and remove group in g and its particles from p and trial
        bool eraseGroup(int);        //!< Remove n'th group as well as its particles
        //int enroll(Group&);          //!< Store group pointer
        void reserve(int);           //!< Reserve space for particles for better memory efficiency
        string info();               //!< Information string

        /** @brief Reset and refill atom- and molecular trackers*/
        inline void initTracker() {
          atomTrack.clear();
          molTrack.clear();
          for (auto g : groupList()) {
            molTrack.insert(g->molId, g);
            for (auto i : *g)
              atomTrack.insert(p[i].id, i);
          }
        }

        /**
         * @brief Find which group given particle index belongs to.
         *
         * If not found, `nullptr` is returned.
         */
        inline Group* findGroup(int i) {
          for (auto g : groupList())
            if (g->find(i))
              return g;
          return nullptr;
        }

        /**
         * @brief Find group index for given group pointer.
         *
         * If not found, `-1` is returned.
         */
        inline int findIndex(Group* group) {
          auto it = std::find(g.begin(), g.end(), group); 
          return (it!=g.end()) ? it-g.begin() : -1;
        }

        /**
         * @brief Finds index of particle using addresses
         *
         * Given a particle reference, this function will try
         * to locate the particle position in the particle
         * vector. Both `p` and `trial` and scanned and in
         * that order. If the particle is not from `Space`
         * -1 will be returned.
         */
        inline int findIndex(const Tparticle &a) const {
          int i=&a-&p.at(0); // calc. index from addresses
          if (i>=0)
            if (i<(int)p.size())
              return i;
          i=&a-&trial.at(0); // calc. index from addresses
          if (i>=0)
            if (i<(int)trial.size())
              return i;
          return -1;
        }

        /** @brief Returns pointer to random molecule of type `molid` */
        inline Group* randomMol( int molId ) {
          auto v = findMolecules( molId );
          if ( !v.empty() )
            return *slump.element( v.begin(), v.end() );
          return nullptr;
        }

        /** @brief Returns vector of molecules with matching molid */
        inline std::vector<Group*> findMolecules( int molId, bool sort=false ) {
          auto v = molTrack[molId];
          if (sort)
            std::sort( v.begin(), v.end(),
                       [](Group* a, Group* b) { return a->front() < b->front(); });
          return v;
/*
          std::vector<Group*> v;
          v.reserve( g.size() );
          for ( auto i : g )
            if (i->molId == molId)
              v.push_back( i );

          // now lets make sure the molecules are sorted according
          // to their position in the particle vector:
          if (sort)
            std::sort( v.begin(), v.end(), 
                [](Group* a, Group* b) { return a->front() < b->front(); }); 
          return v;*/
        }

        /** @brief Count number of molecules with matching molid */
        inline int numMolecules( int molId ) const {
          return molTrack.size(molId);
          //int cnt=0;
          //for (auto i : g)
          //  if (i->molId == molId)
          //    cnt++;
          //return cnt;
        }
 
        /**
         * @brief Returns vector of molecules with matching name
         *
         * @todo Add regex support? Inspiration,
         *
         * ~~~~
         * bool ExpandWildCard(vector<string>& names, vector<string>& result, string& wildcard) {
         *   auto oldsize = result.size();
         *   std::copy_if(std::begin(names), std::end(names),
         *       std::back_inserter(result),
         *       [&](string const& name) {
         *       return std::regex_match(name, make_regex(wildcard)); } );
         *   return (result.size() > oldsize);
         * }
         * ~~~~
         */
        inline std::vector<Group*> findMolecules( const std::string &name, bool sort=false ) {
          auto it = molecule.find(name);
          if ( it != molecule.end() )
            return findMolecules( it->id, sort );
          else return std::vector<Group*>();
        }

    };

  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::checkSanity() {
      assert(p.size()==trial.size() && "Trial/P vector-size mismatch!");
      assert(p==trial && "Trial/P vector mismatch!");
      bool rc=true;
      if (p.size()!=trial.size())
        rc=false;
      else {
        for (size_t i=0; i<p.size(); i++) {
          if ( geo.collision(p[i], p[i].radius) ) {
            rc=false;
            assert(!"Particle/container collision!");
            break;
          }
          if (p[i].x()!=trial[i].x() ||
              p[i].y()!=trial[i].y() ||
              p[i].z()!=trial[i].z() ||
              p[i].charge!=trial[i].charge)
          {
            rc=false;
            std::cerr << "Trial/P vector mismatch!\n";
            break;
          };
        }
        size_t cnt=0;
        for (auto gi : g)
          if (!gi->empty())
            cnt+=gi->size();

        if (cnt!=p.size()) {
          std::cerr << "Sum of enrolled group sizes does not match particle vector";
          rc=false;
        }
      }

      if (atomTrack.size() != p.size())
        rc=false;

      if (molTrack.size() != g.size())
        rc=false;

      if ( rc==false )
        throw std::runtime_error("Space sanity check failed.");

      return rc;
    }

  /**
   * This will insert a particle vector into the end current space (p and trial).
   *
   * @param pin Particle vector to insert
   * @param i Insert position (PRESENTLY IGNORED). Default = -1 which means end of current vector
   *
  template<class Tgeometry, class Tparticle>
    Group Space<Tgeometry,Tparticle>::insert(const p_vec &pin, int i) {
      assert(1==2 && "Function obsolete");
      assert(i==-1 && "Vector insertion at random position unimplemented.");
      Group g;
      if ( !pin.empty() ) {
        g.setrange( p.size(), -1);
        assert(g.size()==0 && "Group range broken!");

        for (size_t i=0; i<pin.size(); i++)
          atomTrack.insert(pin[i].id, p.size()+i);

        p.insert( p.end(), pin.begin(), pin.end() );
        trial.insert( trial.end(), pin.begin(), pin.end() );

        g.resize( pin.size() );
        g.setMassCenter( *this );
        g.setMolSize( pin.size() );
      }
      return g;
    }*/

  /**
   * @param a Particle to insert
   * @param i Insert position in particle vector. Old i will be pushed forward.
   *          Default is -1 = end of vector.
   *
   * This will insert a particle in both `p` and `trial` vectors and
   * expand or push forward groups.
   *
   * @todo REMOVE, still used by `Move::AtomTracker`.
   */
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::insert(const Tparticle &a, int i) {
      if ( i==-1 || i>(int)p.size() ) {
        i=p.size();
        p.push_back(a);
        trial.push_back(a);
      } else {
        p.insert(p.begin()+i, a);
        trial.insert(trial.begin()+i, a);
      }
      atomTrack.insert(a.id, i);

      for (auto gj : g) {
        if ( gj->front() > i ) gj->setfront( gj->front()+1  ); // gj->beg++;
        if ( gj->back() >= i ) gj->setback( gj->back()+1 );    //gj->last++; // +1 is a special case for adding to the end of p-vector
      }
      return true;
    }

  /**
   * @param atomname Name if atom to insert
   * @param n Number of atoms to insert
   * @param key Specify `NOOVERLAPCHECK` if overlap is allowed [default: `OVERLAPCHECK`]
   *
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::insert(string atomname, int n, keys key) {
      Tparticle a;
      a=atom[atomname];
      while (n>0) {
        geo.randompos(a);
        if ( !Geometry::overlap(geo,p,a) || key==NOOVERLAP_CHECK) {
          insert(a,-1);
          n--;
        }
      }
      return true;
    }*/

  /** @brief Erase i'th particle */
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::erase(int i) {
      assert( i < (int)p.size() );
      assert( checkSanity() && "erase(int) fails at sanity check.");

      if (i<(int)p.size()) {
        atomTrack.erase( p[i].id, i );
        p.erase( p.begin()+i );
        trial.erase( trial.begin()+i );

        Group* is_empty=nullptr;
        for (auto gj : g) {               // downshift groups
          if ( i < gj->front() )
            gj->setfront( gj->front()-1);
          if ( i<=gj->back() )
            gj->setback( gj->back()-1);
          if (gj->size()==0)
            is_empty=gj;
        }

        if (is_empty!=nullptr) {          // remove empty group
          molTrack.erase( is_empty->molId, is_empty );
          g.erase( g.begin()+findIndex(is_empty) );
          delete(is_empty);
        }
        assert( checkSanity() );
        return true;
      }
      return false;
    }

  /**
   * This will remove the specified group (given as index in `groupList()`)
   * from the space. Later groups will be shuffled down.
   * The group pointer in `groupList` will be destructed.
   *
   * @todo Optimize by swapping groups of similar size to delete always at
   * the end of `g`.
   */
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::eraseGroup(int i) {

      assert( !groupList().empty() );
      assert( i>=0 && i<(int)g.size() );

      if ( !groupList().empty() ) {

        assert( checkSanity() );
        assert( !g.empty() );
        assert( i < (int)g.size() );

        int n   = g.at(i)->size(); // number of particles in group
        int beg = g[i]->front();   // first particle
        int end = g[i]->back();    // last particle

        // erase from trackers
        molTrack.erase( g[i]->molId, g[i] );
        for (auto i : *g[i])
          atomTrack.erase( p[i].id, i);

        // erase data
        delete( g[i] );            // destruct group
        g.erase( g.begin()+i );    // remove group pointer
        p.erase( p.begin()+beg, p.begin()+end+1); // remove particles
        trial.erase( trial.begin()+beg, trial.begin()+end+1 );

        // move later groups down to reflect new particle index
        size_t cnt=0;
        for (auto l : groupList()) { // loop over groups (pointers)
          cnt+=l->size(); // count particles in each group
          if (l->front()>end) {
            // erase index from atom tracker
            for (auto i : *l)
              atomTrack.erase(p[i].id, i);
            // move index down
            l->setfront( l->front()-n );
            l->setback( l->back()-n );
            // add updated index to atom tracker
            for (auto i : *l)
              atomTrack.insert(p[i].id, i);
          }
        }
        assert(cnt==p.size() && "Particle mismatch while erasing a group!");
        return true;
      }
      return false;
    }

  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::save(string file) {
      using std::numeric_limits;
      if (checkSanity()) {
        cout << "  Writing space state file '" << file << "'. ";
        std::ofstream fout( file.c_str() );
        if (fout) {
          fout.precision( numeric_limits<double>::digits10 + 1 );
          if (std::is_base_of<Geometry::Cuboid,Tgeometry>::value) 
            fout << geo.len.transpose() << "\n";
          else fout << geo.getVolume() << "\n";
          fout << p.size() << "\n";
          for (auto p_i : p)
            fout << p_i << "\n";
          fout << g.size() << "\n";
          for (auto g_i : g)
            fout << *g_i << "\n";
          fout << "randomstate\n";
          fout << slump.eng;
          cout << "OK!\n";
          return true;
        }
      }
      cout << "FAILED!\n";
      return false;
    }

  /**
   * @param file Filename
   * @param key If set to `RESIZE`, `p` and `trial` will be
   *        expanded if they do not match the file
   *        (for Grand Canonical MC)
   */
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::load(string file, keys key) {
      using namespace textio;
      cout << "Reading space state file '" << file << "'. ";
      if (checkSanity()) {
        fin.close();
        fin.open( file.c_str() );
        if (fin) {
          int n;
          cout << "OK!\n";
          if (std::is_base_of<Geometry::Cuboid,Tgeometry>::value) {
            double x,y,z;
            fin >> x >> y >> z >> n;
            geo.setlen(Point(x,y,z));
          }
          else {
            double vol;
            fin >> vol >> n;
            geo.setVolume(vol);
          }
          if (key==RESIZE && n!=(int)p.size()) {
            cout << indent(SUB) << "Resizing particle vector from "
              << p.size() << " --> " << n << ".\n";
            p.resize(n);
          }
          assert(n==(int)p.size());
          if (n == (int)p.size() ) {
            for (int i=0; i<n; i++)
              p[i] << fin;
            trial=p;
            cout << indent(SUB) << "Read " << n << " particle(s)." << endl;
            // Read groups
            fin >> n;
            if (g.empty()) {  // group list is empty!
              g.resize(n);
              for (auto &i : g)
                i = new Group();
            }
            if (n==(int)g.size()) {
              for ( auto i : g ) {
                *i << fin;
                i->setMassCenter(*this);
                if ( i->name.empty() )
                  i->name = molecule[i->molId].name;
              }
              cout << indent(SUB) << "Read " << n << " group(s)." << endl;
            } else {
              std::cerr << "FAILED! (space groups do not match)." << endl;
              return false;
            }
          }
          string id;
          fin >> id;
          if ( id=="randomstate" ) {
            cout << indent(SUB) << "Restoring random number generator state." << endl;
            fin >> slump.eng;
          }

          initTracker(); // update trackers

          return true;
        }
      }

      std::cerr << "State file not found." << endl;
      return false;
    }

  /**
   * Call to this function is *optional* and may provide better
   * handling of memory in cases where the maximum number of
   * particles in the system is known.
   * This is done by reserving the required amount of memory so
   * that dynamic reallocation is avoided upon vector expansion.
   * Calls to this function is best made immediately after Space
   * construction, but can be safely called at any time.
   */
  template<class Tgeometry, class Tparticle>
    void Space<Tgeometry,Tparticle>::reserve(int MaxNumberOfParticles) {
      p.reserve(MaxNumberOfParticles);
      trial.reserve(MaxNumberOfParticles);
    }

  /**
   * This will register a new group in the Space. If already found,
   * nothing is added. In addition, the following is performed each
   * time a new group is enrolled:
   *
   * - The group mass center is re-calculated and set.
   * - Trial vector is synched with the main (`p`) particle vector.
   *
   * @param newgroup Group to enroll (pointer is saved)
   * @returns Position of the newgroup in the group list.
   * @todo Remove and let space handle all groups internally
  template<class Tgeometry, class Tparticle>
    int Space<Tgeometry,Tparticle>::enroll(Group &newgroup) {
      assert(1==2 && "Obsolete function");
      for (size_t i=0; i<g.size(); ++i)
        if (g[i]==&newgroup)
          return i;
      newgroup.setMassCenter(*this);
      g.push_back(&newgroup);
      trial=p;
      return g.size()-1; //return position of added group
    }*/

  template<class Tgeometry, class Tparticle>
    string Space<Tgeometry,Tparticle>::info() {
      using namespace textio;
      static char w=25;
      double z=netCharge( p.begin(), p.end() );
      std::ostringstream o;
      o << header("Simulation Space and Geometry")
        << geo.info(w)
        << pad(SUB,w,"Number of particles") << p.size() << endl
        << pad(SUB,w,"Electroneutrality")
        << ((fabs(z)>1e-7) ? "NO!" : "Yes") << " "  << z << endl
        << pad(SUB,w,"System sanity check")
        << (checkSanity() ? "Passed" : "Failed") << endl
        << pad(SUB,w,"Number of molecule types") << molecule.size() << endl
        << indent(SUB) << "Groups:" << endl;
      for ( size_t i=0; i<g.size(); i++ ) {
        std::ostringstream range;
        range << "[" << g[i]->front() << "-" << g[i]->back() << "]";
        o << indent(SUBSUB) << std::left
          << setw(6) << i+1
          << setw(17) << range.str()
          << setw(12) <<  g[i]->name;
        if ( g[i]->isAtomic() )
          o << "N/V = " << setw(6) << g[i]->numMolecules()/geo.getVolume()
            << _angstrom+superminus+cubed
            << " = "  << setw(6) << g[i]->size()/geo.getVolume()*1e30/pc::Nav
            << " mM  ";
        o << endl;
      }
      return o.str();
    }

  /**
   * This will insert a particle vector into the end of current space.
   * If `MoleculeData::molId` refers to an atomic molecule (`isAtomic()==true`)
   * the inserter will, if present, insert at the end of the last group with matching molid.
   * For molecular groups, a new group is generated and added to `groupList()` and
   * the particles are added to the end of the `p` and `trial`.
   *
   * During addition, the atom- and molecule trackers are being updated.
   *
   * @param molId Molecule ID as defined by the `MoleculeMap`
   * @param pin Particle vector to insert
   */
  template<class Tgeometry, class Tparticle>
    Group* Space<Tgeometry,Tparticle>::insert(PropertyBase::Tid molId, const p_vec &pin) {
      assert( checkSanity() && "insert(id,pvec) fails at sanity check.");
      if ( !pin.empty() ) {
        // insert atomic groups into existing group, if present
        if (molecule[molId].isAtomic() && !g.empty()) {
          int max=-1, imax=-1; 
          for ( size_t i=0; i<g.size(); i++ ) // search for existing group
            if ( g[i]->molId == molId )          // ...and find last
              if ( g[i]->back() > max ) {
                max=g[i]->back();
                imax=i;                     // index of latest group
              }
          // group exists -- now add particles
          if (imax>=0) {
            // add to particle vectors
            p.insert( p.begin() + g[imax]->back()+1, pin.begin(), pin.end() );
            trial.insert( trial.begin() + g[imax]->back()+1, pin.begin(), pin.end() );
            g[imax]->setback( g[imax]->back() + pin.size() );

            // push forward groups above
            for (auto i : g) {
              if (*i > *g[imax]) {
                // erase index from atom tracker
                for (auto ndx : *i)
                  atomTrack.erase(p[ndx].id, ndx);
                // shift upwards
                i->shift(pin.size());
                // add updated index to atom tracker
                for (auto ndx : *i)
                  atomTrack.insert(p[ndx].id, ndx);
              }
            }
            // add to atom tracker
            for (auto i : *g[imax])
              atomTrack.insert( p[i].id, i );

            return g[imax];
          }
        }
        // create and add molecule to end of particle vector
        Group* x = new Group( molecule[molId].name );
        x->molId = molId;
        x->setrange( p.size(), -1 );
        x->resize( pin.size() );

        assert( (size_t)x->front() == p.size() );
        assert( (size_t)x->back()  == p.size() + pin.size() - 1 );

        // add group and particles
        groupList().push_back(x);
        p.insert( p.end(), pin.begin(), pin.end() );
        trial.insert( trial.end(), pin.begin(), pin.end() );

        // add to atom and molecule trackers
        molTrack.insert( x->molId, x );
        for (auto i : *x)
          atomTrack.insert( p[i].id, i);

        if ( molecule[molId].isAtomic() )
          x->setMolSize(1);
        else
          x->setMolSize( pin.size() );

        x->setMassCenter( *this );

        return x;
      }
      return nullptr;
    }

  /*
   * This will remove the specified group
   * from the space. Later groups will be shuffled down.
   *
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::eraseGroup(Group& group) {
      int n = group.size(); // number of particles in group
      int beg = group.front(); // first particle
      int end = group.back();  // last particle
      bool del = false;
      int i=0;
      auto it = g.begin();

      // remove from trackers
      molTrack.erase(group.molId, &group);
      for (auto i : group)
        atomTrack.erase( p[i].id, i);

      // group isnt in group list
      p.erase( p.begin()+beg, p.begin()+end+1); // remove particles -> doest remove partice pointed by last (+1)
      trial.erase( trial.begin()+beg, trial.begin()+end+1 );

      for (auto* l : g) {
        // group overlaps with l -> all 4 possibilities condensed (== ==, < >=, <= >, < >)
        if( (l->front() == beg && l->back() == end) ) {
          it += i;
          del = true;
        } else {
          if (l->front() <= beg && l->back() >= end) {
            // remove index from atom tracker
            for (auto _i : *l)
              atomTrack.erase( p[_i].id, _i);
            // subtract from end
            l->setback( l->back()-n );
            // add again w. updated index
            for (auto _i : *l)
              atomTrack.insert( p[_i].id, _i);
          }
        }

        // later groups
        if (l->front() > end) {
          // remove index from atom tracker
          for (auto _i : *l)
            atomTrack.erase( p[_i].id, _i);

          l->setfront( l->front()-n );
          l->setback( l->back()-n );

          // add again w. updated index
          for (auto _i : *l)
            atomTrack.insert( p[_i].id, _i);
        }

        i++;
      }

      if(del) {
        assert((*it)->front() == group.front() && (*it)->back() == group.back());
        delete *it;
        g.erase( it );// remove group pointer
      }

#ifndef NDEBUG
      size_t cnt=0;
      for (auto* l : g)
        cnt+=l->size(); // count particles in each group
      assert(cnt==p.size() && "Particle mismatch while erasing a group!");
#endif
      return true;
    }*/

} //namespace
#endif
