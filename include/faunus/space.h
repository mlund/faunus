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
        std::map<Tid, std::vector<T> > map;
        std::map<Tid, Average<double> > Navg;

      public:
        /** @brief Number of elements of type id */
        size_t size(Tid id) const {
          auto i=map.find(id);
          if (i!=map.end())
            return i->second.size();
          return 0;
        }
      
        /**
         * @brief Update average number of particles
         */
        void updateAvg() {
          for (auto &m : map)
            Navg[m.first] += m.second.size();
        }
      
        /** @brief Get average number of particles */
        Average<double> getAvg(Tid id) { return Navg[id]; }

        /**
         * @brief Find random data based on id
         * @param id Id for data
         * @param N Number of unique elements (data) to return
         * @param dst Destination vector -- new elements added to end
         * @return True if found; false otherwise
         */
        bool find(Tid id, size_t N, std::vector<T> &dst) {
          size_t Ninit=dst.size();
          auto i=map.find(id);
          if (i!=map.end()) // id found
            if (i->second.size()>=N) {  // enough elements?
              do {
                auto j=randomElement(i->second.begin(), i->second.end());
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
          auto i = map.find(id);
          if (i!=map.end()) {
            auto pos=std::find(i->second.begin(), i->second.end(), data);
            if (pos==i->second.end())
              i->second.push_back(data);
          }
          else
            map[id].push_back(data);
        }

        /**
         * @brief Erase data
         * @param id Id of data
         * @param data Data to be deleted
         */
        void erase(Tid id, T data) {
          auto i = map.find(id);
          if (i!=map.end()) {
            auto pos=std::find(i->second.begin(), i->second.end(), data);
            if (pos!=i->second.end()) {
              std::swap( *pos, i->second.back() );
              i->second.pop_back();
            }
          }
        } 

    };

  // currently not used for anything!
  template<typename T, typename Tparticle=particle>
    class ParticleVector {
      private:
        int _N;
      public:
        typedef Eigen::Vector3d Tvec;
        typedef Eigen::MatrixXd Tmatrix;
        Tmatrix id;
        Tmatrix r,mu,mu_ind;
        Tmatrix charge;

        int size() const { return _N; }

        void resize(int N) {
          assert(N<=id.size());
          _N=N;
        }

        void allocate(int N) {
          int dim=3;
          r.resize(dim,N);
          mu.resize(dim,N);
          mu_ind.resize(dim,N);
          charge.resize(N,1);
          id.resize(N,1);
        }

        Tparticle operator[](int i) const {
          assert(i<=id.size());
          Tparticle p;
          p = atom[ id[i] ];
          p = r.col(i);
          p.charge=charge[i];
          return p;
        }

    };

  /**
   * @brief Placeholder for particles and groups
   *
   * Every simulation must have a `Space` instance as this contains
   * the particles and information about groups (particle ranges).
   * Space also contains the geometry as given by the template
   * parameter.
   */
  template<typename Tgeometry, typename Tparticle=PointParticle>
    class Space {
      protected:
        std::ifstream fin;
      private:
        slump slp;
        bool overlap_container() const;
        bool overlap() const;
        bool overlap(const Tparticle&) const;//!< Check hardspheres overlap with particle
        bool checkSanity();                 //!< Check group length and vector sync
        std::vector<Group*> g;              //!< Pointers to ALL groups in the system

      public:
        typedef std::vector<Tparticle, Eigen::aligned_allocator<particle> > p_vec;
        typedef p_vec ParticleVector;
        typedef Tparticle ParticleType;
        typedef Tgeometry GeometryType;
        enum keys {OVERLAP_CHECK,NOOVERLAP_CHECK,RESIZE,NORESIZE};
        Tgeometry geo;                      //!< System geometry
        p_vec p;                            //!< Main particle vector
        p_vec trial;                        //!< Trial particle vector. 
        std::vector<Group*>& groupList();   //!< Vector with pointers to all groups

        Space(InputMap&);

        void freeGroups() {
          for(auto it = g.begin(); it != g.end(); ++it) delete *it;
          g.clear();
        }

        bool save(string);                  //!< Save container state to disk
        bool load(string, keys=NORESIZE);   //!< Load container state from disk

        /// \brief Sets molId of groups based on MoleculeMap
        void linkGroupsToTopo();

        /// \brief insert p_vec of MolID to end of p and trial
        Group* insert(PropertyBase::Tid, const p_vec&); // inserts to trial and p

        ///
        /// \brief insert p_vec of MolID
        /// \param enlarge - sets whether to enlarge group of MolID, or to add new isMolecular()==true group
        ///
        Group* insert(const p_vec&, PropertyBase::Tid, bool enlarge=true);

        Group insert(const p_vec&, int=-1);
        bool insert(const Tparticle&, int=-1); //!< Insert particle at pos n (old n will be pushed forward).
        bool insert(string, int, keys=OVERLAP_CHECK); 
        bool erase(int);             //!< Remove n'th particle and downshift/remove groups
        bool eraseGroup(Group& group); ///< find and remove group in g and its particles from p and trial
        bool eraseGroup(int);        //!< Remove n'th group as well as its particles
        int enroll(Group&);          //!< Store group pointer
        void reserve(int);           //!< Reserve space for particles for better memory efficiency

        double charge() const;       //!< Sum all charges
        string info();               //!< Information string
        void displace(const Point&); //!< Displace system by a vector

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
    };

  template<class Tgeometry, class Tparticle>
    Space<Tgeometry,Tparticle>::Space(InputMap &in) : geo(in) {}

  template<class Tgeometry, class Tparticle>
    vector<Group*>& Space<Tgeometry,Tparticle>::groupList() { return g; }

  template<class Tgeometry, class Tparticle>
    double Space<Tgeometry,Tparticle>::charge() const {
      double z=0;
      for (auto &p_i : p)
        z+=p_i.charge;
      return z;
    }

  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::checkSanity() {
      assert(p.size()==trial.size() && "Trial/P vector-size mismatch!");
      assert(p==trial && "Trial/P vector mismatch!");
      bool rc=true;
      if (p.size()!=trial.size())
        rc=false;
      else {
        for (size_t i=0; i<p.size(); i++) {
          if ( geo.collision(p[i]) ) {
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

      if (rc==false) {
        std::cerr << "Space sanity check failed. This is serious!";
        std::exit(1);
      }
      return rc;
    }

  /**
   * This will insert a particle vector into the end current space (p and trial).
   *
   * @param pin Particle vector to insert
   * @param i Insert position (PRESENTLY IGNORED). Default = -1 which means end of current vector
   */
  template<class Tgeometry, class Tparticle>
    Group Space<Tgeometry,Tparticle>::insert(const p_vec &pin, int i) {
      assert(i==-1 && "Vector insertion at random position unimplemented.");
      Group g;
      if ( !pin.empty() ) {
        g.setrange( p.size(), -1);
        assert(g.size()==0 && "Group range broken!");

        p.insert( p.end(), pin.begin(), pin.end() );
        trial.insert( trial.end(), pin.begin(), pin.end() );

        g.resize( pin.size() );
        g.setMassCenter( *this );
        g.setMolSize( pin.size() );
      }
      return g;
    }

  /**
   * @param a Particle to insert
   * @param i Insert position in particle vector. Old i will be pushed forward.
   *          Default is -1 = end of vector.
   *
   * This will insert a particle in both `p` and `trial` vectors and
   * expand or push forward any enrolled groups.
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
      for (auto gj : g) {
        if ( gj->front() > i ) gj->setfront( gj->front()+1  ); // gj->beg++;
        if ( gj->back() >= i ) gj->setback( gj->back()+1 );    //gj->last++; // +1 is a special case for adding to the end of p-vector
      }
      return true;
    }

  /**
   * @param atomname Name if atom to intert
   * @param n Number of atoms to insert
   * @param key Specify `NOOVERLAPCHECK` if overlap is allowed [default: `OVERLAPCHECK`]
   */
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::insert(string atomname, int n, keys key) {
      particle a;
      a=atom[atomname];
      while (n>0) {
        geo.randompos(a);
        if (!overlap(a) || key==NOOVERLAP_CHECK) {
          insert(a,-1);
          n--;
        }
      }
      return true;
    }

  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::erase(int i) {
      if (i<(int)p.size()) {
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
          g.erase( g.begin()+findIndex(is_empty) );
          delete(is_empty);
        }
        return true;
      }
      return false;
    }

  /**
   * This will remove the specified group (given as index in `groupList()`)
   * from the space. Later groups will be shufled down.
   * The group pointer in `groupList` will be destructed.
   *
   * @todo Optimize by swapping groups of similar size to delete always at
   * the end of `g`.
   */
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::eraseGroup(int i) {

      if ( !groupList().empty() ) {

        assert( checkSanity() );
        assert( !g.empty() );
        assert( i < (int)g.size() );

        int n   = g.at(i)->size(); // number of particles in group
        int beg = g[i]->front();   // first particle
        int end = g[i]->back();    // last particle
        delete( g[i] );            // destruct group

        g.erase( g.begin()+i );    // remove group pointer
        p.erase( p.begin()+beg, p.begin()+end+1); // remove particles
        trial.erase( trial.begin()+beg, trial.begin()+end+1 );

        // move later groups down to reflect new particle index
        size_t cnt=0;
        for (auto l : groupList()) { // loop over groups (pointers)
          cnt+=l->size(); // count particles in each group
          if (l->front()>end) {
            l->setfront( l->front()-n );
            l->setback( l->back()-n );
          }
        }
        assert(cnt==p.size() && "Particle mismatch while erasing a group!");
        return true;
      }
      return false;
    }

  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::overlap_container() const {
      for (auto &i : p)
        if (geo.collision(i))
          return true;
      return false;
    }

  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::overlap() const {
      for (auto &i : p)
        if (geo.collision(i))
          return true;
      for (auto i=p.begin(); i!=p.end()-1; ++i)
        for (auto j=i+1; j!=p.end(); ++j) {
          double s=i->radius + j->radius;
          if (geo.sqdist(*i,*j)<s*s)
            return true;
        }
      return false;
    }

  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::overlap(const Tparticle &a) const {
      for (auto &p_i : p) {
        double contact = a.radius + p_i.radius;
        if (geo.sqdist(a,p_i) < contact*contact)
          return true;
      }
      return false;
    }

  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::save(string file) {
      using std::numeric_limits;
      if (checkSanity()) {
        cout << "Writing space state file '" << file << "'. ";
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
          double x, y, z, vol;
          cout << "OK!\n";
          if (std::is_base_of<Geometry::Cuboid,Tgeometry>::value) {
            fin >> x >> y >> z >> n;
            geo.setlen(Point(x,y,z));
          }
          else {
            fin >> vol >> n;
            geo.setVolume(vol);
          }
          if (key==RESIZE && n!=(int)p.size()) {
            cout << indent(SUB) << "Resizing particle vector from " << p.size() << " --> " << n << ".\n";
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
              for (auto g_i : g) {
                *g_i << fin;
                g_i->setMassCenter(*this);
                if (g_i->name.empty())
                  g_i->name = molecule[g_i->molId].name;
              }
              cout << indent(SUB) << "Read " << n << " group(s)." << endl;
              return true;
            } else {
              cout << "FAILED! (space groups do not match).\n";
              return false;
            }
          }
          return true;
        }
      }
      cout << "FAILED!\n";
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
   */
  template<class Tgeometry, class Tparticle>
    int Space<Tgeometry,Tparticle>::enroll(Group &newgroup) {
      for (size_t i=0; i<g.size(); ++i)
        if (g[i]==&newgroup)
          return i;
      newgroup.setMassCenter(*this);
      g.push_back(&newgroup);
      trial=p;
      return g.size()-1; //return position of added group
    }

  template<class Tgeometry, class Tparticle>
    string Space<Tgeometry,Tparticle>::info() {
      using namespace textio;
      static char w=25;
      double z=charge();
      std::ostringstream o;
      o << header("Simulation Space and Geometry")
        << geo.info(w)
        << pad(SUB,w,"Number of particles") << p.size() << endl
        << pad(SUB,w,"Electroneutrality")
        << ((abs(z)>1e-7) ? "NO!" : "Yes") << " "  << z << endl
        << pad(SUB,w,"System sanity check")
        << (checkSanity() ? "Passed" : "Failed") << endl
        << indent(SUB) << "Groups:" << endl;
      for (size_t i=0; i<g.size(); i++) {
        std::ostringstream range;
        range << "[" << g[i]->front() << "-" << g[i]->back() << "]";
        o << indent(SUBSUB) << std::left
          << setw(6) << i+1
          << setw(17) << range.str()
          << "N/V = " << setw(6) << g[i]->numMolecules()/geo.getVolume()
          << _angstrom+superminus+cubed
          << " = "  << setw(6) << g[i]->size()/geo.getVolume()*1e30/pc::Nav
          << " mM  " << g[i]->name << endl;
      }
      return o.str();
    }

  /**
   * This will displace the whole system by a vector while
   * respecting boundary conditions of the simulation container.
   * Group mass centers are also updated and the routine can be used to
   * keep certain groups at a particular, absolute position.
   * Applying this function usually does not cause an energy change
   * unless the Hamiltonian depends on absolute positions (for example
   * electric potentials from planar surfaces).
   */
  template<class Tgeometry, class Tparticle>
    void Space<Tgeometry,Tparticle>::displace(const Point &v) {
      for (auto i : g) {
        i->translate(*this, v);
        i->accept(*this);
        assert( geo.sqdist(i->cm, i->cm_trial) < 1e-6 );
      }
    }

  /**
   * This will insert a particle vector into the end of current space.
   * If `isAtomic()==true` for the molecule, the inserter will, if present,
   * insert at the end of the last group with matching molid.
   */
  template<class Tgeometry, class Tparticle>
    Group* Space<Tgeometry,Tparticle>::insert(PropertyBase::Tid molId, const p_vec &pin) {
      if ( !pin.empty() ) {
        int nold=p.size();

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
            cout << "Adding to existing group..." << endl;
            p.insert( p.begin() + g[imax]->back()+1, pin.begin(), pin.end() );
            trial.insert( trial.begin() + g[imax]->back()+1, pin.begin(), pin.end() );
            g[imax]->setback( g[imax]->back() + pin.size() );
            for (auto i : g)
              if (*i > *g[imax])
                i->shift( pin.size() );

            assert(nold+pin.size()==p.size());
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

        groupList().push_back(x);
        p.insert( p.end(), pin.begin(), pin.end() );
        trial.insert( trial.end(), pin.begin(), pin.end() );

        if ( molecule[molId].isAtomic() )
          x->setMolSize(1);
        else x->setMolSize( pin.size() );

        x->setMassCenter( *this );

        return x;
      }
      return nullptr;
    }

  /**
   * This will insert a particle vector into the current space at the end of its group.
   *
   * @param pin Particle vector to insert
   * @param molId Id of molecule
   * @param enlarge bool
   * @return new Group
   */
  template<class Tgeometry, class Tparticle>
    Group* Space<Tgeometry,Tparticle>::insert(const p_vec &pin, PropertyBase::Tid molId, bool enlarge) {

      assert(!pin.empty());

      bool found = false;
      Group* retGroup = NULL;
      unsigned int last=std::numeric_limits<unsigned int>::max();
      for(unsigned int i=0; i<g.size(); i++) if(g[i]->molId == molId) last=i;

      for(unsigned int i=last; i<g.size(); i++) {
        if(found && ((i!=last+1) || enlarge)) { // when not enlarging group -> dont overwrite back of newly inserted
          g[i]->setfront(g[i]->front() + pin.size());
          g[i]->setback(g[i]->back() + pin.size());
        }

        if(!found && i==last) {

          // add particles to particle vectors
          p.insert(p.begin() + g[i]->back() +1, pin.begin(), pin.end());
          trial.insert(trial.begin() + g[i]->back() +1, pin.begin(), pin.end());
          found = true;

          if(enlarge) {
            // change group range
            g[i]->setback(g[i]->back() + pin.size());
            retGroup = g[i];
            g[i]->setMassCenter(*this);

          } else {
            retGroup = new Group(g[i]->name);
            retGroup->setback(g[i]->back() + pin.size());
            retGroup->setfront(g[i]->back()+1);
            retGroup->setMolSize(g[i]->getMolSize());
            retGroup->name = g[i]->name;
            retGroup->molId = molId;

            retGroup->setMassCenter(*this);

            g.insert(g.begin() + i+1, retGroup);
          }
        }
      }

      if(!found) { // no groups of molId -> adding at end of g
        Group* newGroup = insert(molId, pin);
        newGroup->name = newGroup->name = molecule[molId].name;
        g.push_back(newGroup);
        return newGroup;
      }
      assert(retGroup != NULL);
      return retGroup;
    }


  /**
   * This will remove the specified group
   * from the space. Later groups will be shufled down.
   *
   */
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::eraseGroup(Group& group) {
      int n = group.size(); // number of particles in group
      int beg = group.front(); // first particle
      int end = group.back();  // last particle
      bool del = false;
      int i=0;
      auto it = g.begin();

      // group isnt in group list
      p.erase( p.begin()+beg, p.begin()+end+1); // remove particles -> doest remove partice pointed by last (+1)
      trial.erase( trial.begin()+beg, trial.begin()+end+1 );

      for (auto* l : g) {
        // group overlaps with l -> all 4 possibilities condensed (== ==, < >=, <= >, < >)
        if( (l->front() == beg && l->back() == end) ) {
          it += i;
          del = true;
        } else { if(l->front() <= beg && l->back() >= end) {
          l->setback( l->back()-n );
        }
        }

        // later groups
        if (l->front() > end) {
          l->setfront( l->front()-n );
          l->setback( l->back()-n );
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
    }

} //namespace
#endif
