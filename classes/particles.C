#include "particles.h"

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
};
bool particles::overlap(particle &a) {
  for (unsigned short i=0; i<p.size(); i++)
    if (p[i].overlap(a)==true) return true;
  return false;
}
