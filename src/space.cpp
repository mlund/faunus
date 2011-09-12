#include "faunus/container.h"
#include "faunus/inputfile.h"
#include "faunus/species.h"
#include "faunus/physconst.h"
#include "faunus/group.h"
#include "faunus/point.h"
#include <faunus/space.h>

namespace Faunus {

  space::space(Geometry::geometrybase &geoPtr) {
    geo=&geoPtr;
  }

  double space::charge() const {
    double z=0;
    for (unsigned short i=0; i<p.size(); i++)
      z+=p[i].charge;
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
  bool space::insert(particle a, unsigned int i) {
    if (i>p.size())
      return false;
    p.insert(p.begin()+i, a);
    trial.insert(trial.begin()+i, a);

    for (int j=0; j<g.size(); j++) {  // move and expand groups if appropriate
      if ( i<g[j]->beg )
        g[j]->beg++;
      if ( i<=g[j]->end )
        g[j]->end++;
    }
    return true;
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
      for (int i=0; i<p.size(); i++)
        fout << p[i] << endl;
      fout << g.size() << endl;
      for (int i=0; i<g.size(); i++)
        fout << *g[i] << endl;
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
        for (int i=0; i<g.size(); i++)
          *g[i] << fin;
        std::cout << "# Read " << n << " groups from " << file << endl;
      }
    }
    std::cerr << "# Container data NOT read from file " << file << endl;
    return false;
  }

  string space::info() {
    static char w=25;
    double z=charge();
    std::ostringstream o;
    o << endl
      << "# SIMULATION SPACE:" << endl
      << geo->info(w);
    geo->pad(o,w); o << "Number of particles" << p.size() << endl;
    geo->pad(o,w); o << "Volume (AA^3)" << geo->getvolume() << endl;
    geo->pad(o,w); o << "Electroneutrality" << ((abs(z)>1e-7) ? "NO!" : "Yes") << " "  << z << endl;
    return o.str();
  }
}//namespace
