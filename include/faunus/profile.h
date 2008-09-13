#ifndef FAU_PROFILE_H
#define FAU_PROFILE_H

#include "faunus/xytable.h"
#include "faunus/io.h"

namespace Faunus {

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
    std::ostringstream o;
    for (float d=xmin; d<xmax(); d+=xres)
      o << d << " " << conc(d) << endl;
    return fio.writefile(name,o.str());
  }

  /*!\brief Cummulative sum around a central particle
   * \author Mikael Lund
   * \date Canberra, 2008
   * \warning The central particle is passed through the constructor and
   *          it is assumed that the particle (usually in a vector) is not
   *          reallocated. Hence, do NOT modify the particle vector after
   *          having called the constructor.
   */
  class cummsum : public profile {
    protected:
      particle *origo;         //!< Central particle (coordinate origo)
      float volume(double z) { return 1.; }
    public:
      particle::type id;       //!< Particle type to analyse
      cummsum(
          particle::type type, particle &center, float max,float res=.5) :
        profile(0,max,res) {
          id=type;
          origo=&center;
        }
      void add(particle &p) { if (p.id==id) (*this)(origo->dist(p))++; }
  };

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
}
#endif
