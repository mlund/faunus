#ifndef FAUNUS_SPACE_H
#define FAUNUS_SPACE_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/slump.h>
#include <faunus/point.h>

//#include <faunus/geometry.h>
#include <faunus/inputfile.h>
#include <faunus/species.h>
#include <faunus/physconst.h>
#include <faunus/group.h>
#include <faunus/point.h>
#include <faunus/space.h>
#include <faunus/textio.h>
#endif

namespace Faunus {
  
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

        bool save(string);                  //!< Save container state to disk
        bool load(string, keys=NORESIZE);   //!< Load container state from disk

        Group insert(const p_vec&, int=-1);
        bool insert(const Tparticle&, int=-1); //!< Insert particle at pos n (old n will be pushed forward).
        bool insert(string, int, keys=OVERLAP_CHECK); 
        bool erase(int);             //!< Remove n'th particle
        bool eraseGroup(int);        //!< Remove n'th group as well as its particles
        int enroll(Group&);          //!< Store group pointer
        void reserve(int);           //!< Reserve space for particles for better memory efficiency

        double charge() const;       //!< Sum all charges
        string info();               //!< Information string
        void displace(const Point&); //!< Displace system by a vector

        /**
         * @brief Find which group given particle index belongs to
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
         * @brief Finds index of particle using addresses
         *
         * Given a particle reference, this function will try
         * to locate the particle position in the particle
         * vector. Both `p` and `trial` and scanned and in
         * that order. If the particle is not from `Space`
         * -1 will be returned.
         */
        int findIndex(const Tparticle &a) const {
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

        inline double dist(int i, int j) const {
          return geo.dist(p[i],p[j]);
        }
    };

  template<class Tgeometry, class Tparticle>
    Space<Tgeometry,Tparticle>::Space(InputMap &in) : geo(in) {}

  template<class Tgeometry, class Tparticle>
    vector<Group*>& Space<Tgeometry,Tparticle>::groupList() {
      return g;
    }

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
      bool rc=false;
      if (p.size()==trial.size())
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
            assert(!"Trial/P vector mismatch!");
            break;
          } else rc=true;
        }
      size_t cnt=0;
      for (auto gi : g)
        if (!gi->empty())
          cnt+=gi->size();

      if (cnt!=p.size()) {
        assert(!"Sum enrolled group sizes does not match particle vector");
        std::cerr << "Space sanity check failed. This is serious!";
        rc=false;
      }

      if (rc==false) {
        assert(!"Space sanity check failed. This is serious!");
        std::cerr << "Space sanity check failed. This is serious!";
      }
      if (rc==false)
        std::exit(1);
      return rc;
    }

  /**
   * This will insert a particle vector into the current space. No overlap checks are performed; this should
   * be done prior to insertion by for example the `Geometry::FindSpace` class.
   *
   * @param pin Particle vector to insert
   * @param i Insert position (PRESENTLY IGNORED). Default = -1 which means end of current vector
   *
   * @todo Implement insertion at random position
   */
  template<class Tgeometry, class Tparticle>
    Group Space<Tgeometry,Tparticle>::insert(const p_vec &pin, int i) {
      assert(i==-1 && "Vector insertion at random position unimplemented.");
      Group g;
      if ( !pin.empty() ) {
        g.setrange( p.size(), -1);
        assert(g.size()==0 && "Group range broken!");

        p.insert(p.end(), pin.begin(), pin.end());
        g.resize(pin.size());

        //for (auto &i : pin) {
        //  p.push_back(i);
        //  trial.push_back(i);
        //  g.resize( g.size()+1 );
        //}

        g.setMassCenter(*this);
        g.setMolSize(pin.size());
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
      assert( !p.empty() && i<(int)p.size() &&
          "Tried to delete non-existing particle or particle vector empty!" );
      if (i>=(int)p.size())
        return false;
      p.erase( p.begin()+i );
      trial.erase( trial.begin()+i );
      for (auto gj : g) {
        if ( i<gj->front() ) gj->setfront( gj->front()-1  ); // gj->beg--;
        if ( i<=gj->back() ) gj->setback( gj->back()-1);     //gj->last--;
        assert( gj->back()>=0 && "Particle removal resulted in empty Group");
      }
      return true;
    }

  /**
   * This will remove the specified group (given as index in `groupList()`)
   * from the space. Later groups will be shufled down.
   * The `groupList()` contains pointers to groups and it is important
   * that you clean up the deleted group only *after* having called this function.
   *
   * @warning untested
   */
  template<class Tgeometry, class Tparticle>
    bool Space<Tgeometry,Tparticle>::eraseGroup(int i) {
      int n=g.at(i)->size(); // number of particles in group
      int beg=g[i]->front(); // first particle
      int end=g[i]->back();  // last particle

      g.erase( g.begin()+i );// remove group pointer
      p.erase( p.begin()+beg, p.begin()+end); // remove particles
      trial.erase( trial.begin()+beg, trial.begin()+end );

      // move later groups down to reflect new particle index
      size_t cnt=0;
      for (auto l : g) { // loop over groups (pointers)
        cnt+=l->size(); // count particles in each group
        if (l->front()>end) {
          l->setfront( l->front()-n );
          l->setback( l->back()-n );
        }
      }
      assert(cnt==p.size() && "Particle mismatch while erasing a group!");
      return true;
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
          fout << geo.getVolume() << "\n"
            << p.size() << "\n";
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
          double vol;
          cout << "OK!\n";
          fin >> vol >> n;
          geo.setVolume(vol);
          if (key==RESIZE && n!=(int)p.size()) {
            cout << indent(SUB) << "Resizing particle vector from " << p.size() << " --> " << n << ".\n";
            p.resize(n);
          }
          if (n == (int)p.size() ) {
            for (int i=0; i<n; i++)
              p[i] << fin;
            trial=p;
            cout << indent(SUB) << "Read " << n << " particle(s)." << endl;
            fin >> n;
            if (n==(int)g.size()) {
              for (auto g_i : g) {
                *g_i << fin;
                g_i->setMassCenter(*this);
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

} //namespace
#endif
