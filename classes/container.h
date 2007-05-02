#ifndef _CONTAINER_H
#define _CONTAINER_H

#include "particles.h"
#include "slump.h"
#include "point.h"
#include "group.h"
#include "species.h"
#include "hardsphere.h"

/*! \brief Polymorphic class for simulation containers
 *  \author Mikael Lund
 */
class container : public particles,  public species {
  protected:
    slump slp;
  public:
    float volume;                               //!< Volume of the container [AA^3]
    virtual void randompos(point &)=0;          //!< Random point within container
    virtual bool collision(point &)=0;          //!< Check for collision with walls
    virtual group insert(particle::type, short);//!< Insert particles
    virtual string info() {
      float z=charge();
      ostringstream o;
      o << "# Container:" << endl
        << "#   Number of particles  = " << p.size() << endl
        << "#   Volume (AA^3)        = " << volume << endl
        << "#   Electroneutrality    = " 
        << ((z!=0) ? "NO!" : "Yes") << " "  << z << endl;;
      return o.str();
    }
};

/*! \brief Spherical simulation container
 *  \author Mikael Lund
 */
class cell : public container {
  private:
    float r2,diameter;
  public:
    float r;              //!< Radius
    cell(float);
    string info();
    void randompos(point &);
    bool collision(point &p) {
      return 
        (p.x*p.x+p.y*p.y+p.z*p.z > r2) ? true:false;
    }
};

/*! \brief Cubic simulation container
 *  \author Mikael Lund
 *  \todo Not finished!
 */
class box : public container {
  public:
    float len; //!< Side length
    box(float);
    void randompos(point &);
    inline bool collision(point &p) {};
};

/*! \brief Cylindrical simulation container
 *  \author Mikael Lund
 *  \todo Not finished!
 */
class cylinder : public container {
  public:
    float len; //!< Cylinder length
    float r;   //!< Cylinder radius
    cylinder(float);
    void randompos(point &);
    bool collision(point &p) {};
};
#endif
