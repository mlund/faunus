#ifndef _CONTAINER_H
#define _CONTAINER_H

#include "slump.h"
#include "point.h"

/*! \brief Polymorphic class for simulation containers
 *  \author Mikael Lund
 */
class container : public particles {
  protected:
    slump slp;
  public:
    float volume;                               //!< Volume of the container [AA^3]
    virtual void randompos(point &)=0;          //!< Random point within container
    inline virtual bool collision(point &)=0;   //!< Check for collision with walls
    virtual bool insert(particle, short);       //!< Insert particle
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
    void randompos(point &);
    inline bool collision(point &p) {
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
    bool collision(point &p) {};
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
