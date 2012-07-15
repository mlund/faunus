#include "faunus/geometry.h"
#include "faunus/inputfile.h"
#include "faunus/species.h"
#include "faunus/physconst.h"
#include "faunus/group.h"
#include "faunus/point.h"
#include <faunus/space.h>
#include <faunus/textio.h>

namespace Faunus {

  Space::Space(Geometry::Geometrybase &geoPtr) {
    assert(&geoPtr!=nullptr && "Space must have a well-defined geometry!");
    geo=&geoPtr;
  }

  Space::~Space() {}

  vector<Group*>& Space::groupList() {
    return g;
  }

  double Space::charge() const {
    double z=0;
    for (auto &p_i : p)
      z+=p_i.charge;
    return z;
  }

  bool Space::checkSanity() {
    assert(geo!=nullptr && "Space must have a well-defined geometry!");
    assert(p.size()==trial.size() && "Trial/P vector mismatch!");
    assert(p==trial);
    bool rc=false;
    if (p.size()==trial.size())
      for (size_t i=0; i<p.size(); i++) {
        if ( geo->collision(p[i]) ) {
          rc=false;
          assert(!"Particle/container collision!");
          break;
        }
        if (p[i].x!=trial[i].x ||
            p[i].y!=trial[i].y ||
            p[i].z!=trial[i].z ||
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
    return rc;
  }

  /*!
   * This will insert a particle vector into the current space. No overlap checks are performed; this should
   * be done prior to insertion by for example the Geometry::FindSpace class.
   *
   * \param pin Particle vector to insert
   * \param i Insert position (PRESENTLY IGNORED). Default = -1 which means end of current vector
   *
   * \todo Implement insertion at random position
   */
  GroupMolecular Space::insert(const p_vec &pin, int i) {
    assert(i==-1 && "Vector insertion at random position unimplemented.");
    GroupMolecular g;
    if ( !pin.empty() ) {
      g.setrange( p.size(), -1);
      assert(g.size()==0 && "Group range broken!");
      for (auto &i : pin) {
        p.push_back(i);
        trial.push_back(i);
        g.resize( g.size()+1 );
      }
      g.setMassCenter(*this);
    }
    return g;
  }

  /*!
   * \param a Particle to insert
   * \param i Insert position in particle vector. Old i will be pushed forward. Default is -1 = end of vector.
   *
   * This will insert a particle in both \c p and \c trial vectors and expand or push forward any enrolled groups.
   */
  bool Space::insert(const particle &a, int i) {
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
    //cout << "!=" << i << endl;
    //for (auto gj : g)
    //  cout << *gj << "     " << gj->size() << endl;
    return true;
  }

  /*!
   * \param atomname Name if atom to intert
   * \param n Number of atoms to insert
   * \param key Ignored for now -- overlap check is always performed
   */
  bool Space::insert(string atomname, int n, keys key) {
    particle a;
    a=atom[atomname];
    while (n>0) {
      geo->randompos(a);
      if (!overlap(a)) {
        insert(a,-1);
        n--;
      }
    }
    return true;
  }

  bool Space::erase(int i) {
    assert( !p.empty() && i<(int)p.size() && "Tried to delete non-existing particle or particle vector empty!" );
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

  bool Space::overlap_container() const {
    for (auto &i : p)
      if (geo->collision(i))
        return true;
    return false;
  }

  bool Space::overlap() const {
    for (auto &i : p)
      if (geo->collision(i))
        return true;
    for (auto i=p.begin(); i!=p.end()-1; ++i)
      for (auto j=i+1; j!=p.end(); ++j) {
        double s=i->radius + j->radius;
        if (geo->sqdist(*i,*j)<s*s)
          return true;
      }
    return false;
  }

  bool Space::overlap(const particle &a) const {
    for (auto &p_i : p) {
      double contact = a.radius + p_i.radius;
      if (geo->sqdist(a,p_i) < contact*contact)
        return true;
    }
    return false;
  }

  bool Space::save(string file) {
    using std::numeric_limits;
    if (checkSanity()) {
      cout << "Writing space state file '" << file << "'. ";
      std::ofstream fout( file.c_str() );
      if (fout) {
        fout.precision( numeric_limits<double>::digits10 + 1 );
        fout << geo->getVolume() << endl
          << p.size() << endl;
        for (auto p_i : p)
          fout << p_i << endl;
        fout << g.size() << endl;
        for (auto g_i : g)
          fout << *g_i << endl;
        fout.close();
        cout << "OK!\n";
        return true;
      }
    }
    cout << "FAILED!\n";
    return false;
  }

  /*!
   * \param file Filename
   * \param key If set to RESIZE the \c p and \c trial will be expanded if they do not match the file (for Grand Canonical MC)
   */
  bool Space::load(string file, keys key) {
    using namespace textio;
    double vol;
    int n;
    cout << "Reading space state file '" << file << "'. ";
    if (checkSanity()) {
      fin.close();
      fin.open( file.c_str() );
      if (fin) {
        cout << "OK!\n";
        fin >> vol >> n;
        geo->setVolume(vol);
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

  /*!
   * This will register a new group in the Space. If already found, nothing is added. In addition, the
   * following is performed each time a new group is enroller:
   * \li The group mass center is re-calculated and set.
   * \li Trial vector is synched with the main (p) particle vector.
   *
   * \param newgroup Group to enroll (pointer is saved)
   * \returns Position of the newgroup in the group list.
   */
  int Space::enroll(Group &newgroup) {
    for (size_t i=0; i<g.size(); ++i)
      if (g[i]==&newgroup)
        return i;
    newgroup.setMassCenter(*this);
    g.push_back(&newgroup);
    trial=p;
    return g.size()-1; //return position of added group
  }

  string Space::info() {
    using namespace textio;
    static char w=25;
    double z=charge();
    std::ostringstream o;
    o << header("Simulation Space and Geometry") 
      << geo->info(w)
      << pad(SUB,w,"Number of particles") << p.size() << endl
      << pad(SUB,w,"Electroneutrality") << ((abs(z)>1e-7) ? "NO!" : "Yes") << " "  << z << endl
      << pad(SUB,w,"System sanity check") << (checkSanity() ? "Passed" : "Failed") << endl
      << indent(SUB) << "Groups:" << endl;
    for (size_t i=0; i<g.size(); i++) {
      std::ostringstream range;
      range << "[" << g[i]->front() << "-" << g[i]->back() << "]";
      o << indent(SUBSUB) << std::left
        << setw(6) << i+1
        << setw(17) << range.str()
        << setw(20) << g[i]->name
        << endl;
    }

    return o.str();
  }

}//namespace
