#include "faunus/particles.h"
namespace Faunus {
  //! \return Size of particle vector after addition
  int particles::push_back(particle &par) {
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
  double particles::charge(point &origo, double r) {
    double q=0,r2=r*r;
    for (int i=0; i<p.size(); i++)
      if (p[i].sqdist(origo) <= r2)
        q += p[i].charge;
    return q;
  }
  bool particles::overlap(particle &a) {
    for (unsigned short i=0; i<p.size(); i++)
      if (p[i].overlap(a)==true) return true;
    return false;
  }
  bool particles::overlap(vector<particle> &v) {
    unsigned short i = v.size();
    for (unsigned int j=0; j<p.size(); j++) {
      for (unsigned int k=0; k<i; k++) {
        if(p[j].overlap(v[k])==true) return true;
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
  int particles::count(particle::type id, point &origo, double r) {
    int i,cnt=0,n=p.size();
    double r2=r*r;
    for (i=0; i<n; i++)
      if (p[i].sqdist(origo) < r2) cnt++;
    return cnt;
  }
}//namespace
