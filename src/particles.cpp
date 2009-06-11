#include "faunus/particles.h"
namespace Faunus {
  //! \return Size of particle vector after addition
  int particles::push_back(const particle &par) {
    p.push_back(par);
    trial.push_back(par);
    return p.size();
  }

  double particles::charge() {
    double z=0;
    for (unsigned short i=0; i<p.size(); i++)
      z+=p[i].charge;
    return z;
  }

  /*!\param origo Center of the spherical region
   * \param r Radius of the spherical region
   */
  double particles::charge(const point &origo, double r) {
    double q=0,r2=r*r;
    for (int i=0; i<p.size(); i++)
      if (p[i].sqdist(origo) <= r2)
        q += p[i].charge;
    return q;
  }
  /*!\param test particle
   */
  bool particles::overlap(const particle &a) {
    for (unsigned short i=0; i<p.size(); i++)
      if ( clash(p[i],a)==true) return true;
    return false;
  }
  bool particles::overlap(const vector<particle> &v) {
    unsigned short i = v.size();
    for (unsigned int j=0; j<p.size(); j++) {
      for (unsigned int k=0; k<i; k++) {
        if( clash(p[j],v[k])==true ) return true;
      }
    }
    return false;
  }
  bool particles::check_vector() {
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
      std::cout << "# Fatal error: Particle vectors corrupted!!\n";
    return rc;
  }
  int particles::count(unsigned char id, const point &origo, double r) {
    int i,cnt=0,n=p.size();
    double r2=r*r;
    for (i=0; i<n; i++)
      if (p[i].sqdist(origo) < r2) cnt++;
    return cnt;
  }

  vector<unsigned char> particles::list_of_species() const {
    vector<unsigned char> id;
    for (int i=0; i<p.size(); i++) {
      vector<unsigned char>::iterator iter = std::find(id.begin(), id.end(), p[i].id);
      if (iter==id.end())
        id.push_back(p[i].id);
    }
    return id;
  }

  bool particles::insert(particle a, unsigned int i) {
    if (i>p.size())
      return false;
    p.insert(p.begin()+i, a);
    trial.insert(trial.begin()+i, a);
    return true;
  };

  bool particles::remove(unsigned int i) {
    p.erase( p.begin()+i );
    trial.erase( trial.begin()+i );
    return true;
  }

  bool particles::clash(const particle &p1, const particle &p2) {
    return p1.overlap(p2);
  }
}//namespace
