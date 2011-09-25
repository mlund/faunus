#include "faunus/geometry.h"
#include "faunus/inputfile.h"
#include "faunus/species.h"
#include "faunus/physconst.h"
#include "faunus/group.h"
#include "faunus/point.h"
#include <faunus/space.h>
#include <faunus/faunus.h>

namespace Faunus {

  space::space(Geometry::geometrybase &geoPtr) {
    geo=&geoPtr;
  }

  double space::charge() const {
    double z=0;
    for (auto &p_i : p)
      z+=p_i.charge;
    return z;
  }

  bool space::check_vector() {
    bool rc=false;
    if (p.size()==trial.size())
      for (unsigned short i=0; i<p.size(); i++) {
        if (p[i].x!=trial[i].x ||
            p[i].y!=trial[i].y ||
            p[i].z!=trial[i].z ||
            p[i].charge!=trial[i].charge)
        {
          rc=false;
          break;
        } else rc=true;
      }
    if (rc==false)
      std::cerr << "# Fatal error: Particle vectors corrupted!!\n";
    return rc;
  }

  /*!
   * \param a Particle to insers
   * \param i Position in particle vector
   */
  bool space::insert(particle a, int i) {
    if ( i==-1 || i>p.size() )
      i=p.size();
    p.insert(p.begin()+i, a);
    trial.insert(trial.begin()+i, a);

    for (auto g_j : g) {
      if ( i < g_j->beg )
        g_j->beg++;
      if ( i <= g_j->end )
        g_j->end++;
    }
    return true;
  }

  bool space::insert(string name, int n, spacekeys key) {
    particle a;
    a=atom[name];
    while (n>0) {
      geo->randompos(a);
      if (!overlap(a)) {
        insert(a,-1);
        n--;
      }
    }
    return true;
  }

  bool space::overlap(const particle &a) const {
    for (auto &p_i : p) {
      double contact = a.radius + p_i.radius;
      if (geo->_sqdist(a,p_i) < contact*contact)
        return true;
    }
    return false;
  }

  bool space::remove(unsigned int i) {
    if (i>=p.size())
      return false;
    p.erase( p.begin()+i );
    trial.erase( trial.begin()+i );
    for (int j=0; j<g.size(); j++) {  // move and reduce groups if appropriate
      if ( i<g[j]->beg )
        g[j]->beg--;
      if (i<=g[j]->end)
        g[j]->end--;
    }
    return true;
  }

  bool space::save(string file) {
    std::ofstream fout( file.c_str() );
    if (fout) {
      fout.precision(10);
      fout << p.size() << endl;
      for (auto p_i : p)
        fout << p_i << endl;
      fout << g.size() << endl;
      for (auto g_i : g)
        fout << *g_i << endl;
      fout.close();
      return true;
    }
    return false;
  }

  /*!
   * \param file Filename
   * \param resize True if the current geometry should be resized to match file content (default: false)
   */
  bool space::load(string file, bool resize) {
    unsigned int n;
    fin.close();
    fin.open( file.c_str() );
    if (fin) {
      fin >> n;
      if (resize==true)
        p.resize(n);
      if (n==p.size()) {
        for (int i=0; i<n; i++)
          p[i] << fin;
        trial=p;
        std::cout << "# Read " << n << " space points from " << file << endl;

        fin >> n;
        g.resize(n);
        for (auto g_i : g)
          *g_i << fin;
        std::cout << "# Read " << n << " groups from " << file << endl;
        return true;
      }
    }
    std::cerr << "# Space data NOT read from file " << file << endl;
    return false;
  }

  int space::enroll(group &newgroup) {
    for (int i=0; i<g.size(); ++i)
      if (g[i]==&newgroup)
        return i;
    g.push_back(&newgroup);
    return g.size()-1; //return position if added group
  }

  string space::info() {
    using namespace textio;
    static char w=25;
    double z=charge();
    std::ostringstream o;
    o << header("Simulation Space and Geometry") 
      << geo->info(w);
    o << pad("Number of particles",w,SUB) << p.size() << endl;
    o << pad("Number of groups",w,SUB) << g.size() << endl;
    for (int i=0; i<g.size(); i++) {
      o << std::left << std::setw(w);
      o << "" << std::left << setw(4) << i << setw(20) << g[i]->name
        << "[" << g[i]->beg << "-" << g[i]->end << "]"
        << endl;
    }
    o << pad("Volume (AA^3)",w,SUB) << geo->getvolume() << endl;
    o << pad("Electroneutrality",w,SUB) << ((abs(z)>1e-7) ? "NO!" : "Yes") << " "  << z << endl;
    return o.str();
  }
}//namespace
