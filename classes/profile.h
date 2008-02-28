#ifndef _PROFILE_H
#define _PROFILE_H

#include "average.h"

class profile : protected xytable<float,average<unsigned int> >{
  protected:
    virtual float volume(double)=0;
  public:
    profile(float min, float max, float res=0.5) :
      xytable<float,average<unsigned int> >(res,min,max) { }
    virtual void add(particle &)=0; //!< Add a particle
    void update(vector<particle> &);//!< Search for and add particles
    float conc(float);              //!< Get concentration at coordinate
    void show();                    //!< Print distribution
};
float profile::conc(float x) {
  return (*this)(x).sum / ( (*this)(x).cnt * volume(x) );
}
void profile::update(vector<particle> &p) {
  for (int i=0; i<p.size(); i++) add(p[i]);
}
void profile::show() {
  cout << "r conc." << endl;
  cout << xmin << " " << xmax() << endl;
  for (float d=xmin; d<xmax(); d+=xres)
    cout << d << " " << conc(d) << endl;
}

/*!\brief Cylindrical particle distribution
 * \author Mikael Lund
 * \date Canberra, 2008
 * Calculates the particle density in a cylinder around
 * the z-axis.
 */
class cylindric_profile : public profile {
  protected:
    float volume(double z) { return xres*acos(-1.)*r*r; }
  public:
    double r;                //!< Radius of the cylinder
    particle::type id;       //!< Particle type to analyse
    cylindric_profile(
        float radius, particle::type type, float min,float max,float res=.5) :
      profile(min,max,res) {
        r=radius;
        id=type;
      }
    void add(particle &);
};
void cylindric_profile::add(particle &p) {
  if (p.id==id)
    if (p.z>=xmin && p.z<=xmax())
      if (p.x*p.x+p.y*p.y<=r*r)
        (*this)(p.z)+=1;
}
#endif
