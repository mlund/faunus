#ifndef _PROFILE_H
#define _PROFILE_H

#include "xytable.h"
#include "io.h"

/*!\brief Particle profile base class
 * \author Mikael Lund
 * \date Canberra, 2008
 *
 * This is a base class for analysing particle distributions
 * along some one-dimensional coordinate. Derived classes should
 * implement an add() function that, given, a particle adds it
 * to a histogram if certain criteria are fulfulled.
 */
class profile : protected xytable<float,unsigned long int>{
  protected:
    unsigned int cnt;
    virtual float volume(double)=0; //!< Get volume at coordinate
  public:
    profile(float min, float max, float res=0.5) :
      xytable<float,unsigned long int>(res,min,max) { cnt=0; }
    virtual void add(particle &)=0; //!< Add a particle
    void update(vector<particle> &);//!< Search for and add particles
    float conc(float);              //!< Get concentration at coordinate
    bool write(string);             //!< Print distribution
};
float profile::conc(float x) { return ((*this)(x)>0) ? (*this)(x)/(cnt*volume(x)) : 0; }
void profile::update(vector<particle> &p) { for (int i=0; i<p.size(); i++) add(p[i]); }
bool profile::write(string name) {
  io fio;
  ostringstream o;
  for (float d=xmin; d<xmax(); d+=xres)
    o << d << " " << conc(d) << endl;
  return fio.writefile(name,o.str());
}

/*!\brief Cylindrical particle distribution
 * \author Mikael Lund
 * \date Canberra, 2008
 *
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
      if (p.x*p.x+p.y*p.y<=r*r) {
        (*this)(p.z)++;
        cnt++;
      }
}
#endif
