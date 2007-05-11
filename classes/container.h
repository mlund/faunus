#ifndef _CONTAINER_H
#define _CONTAINER_H

#include "particles.h"
#include "slump.h"
#include "point.h"
#include "species.h"

/*! \brief Polymorphic class for simulation containers
 *  \author Mikael Lund
 */
class container : public particles,  public species {
  protected:
    slump slp;
  public:
    float volume;                               //!< Volume of the container [AA^3]
    inline virtual bool collision(point &)=0;   //!< Check for collision with walls
    virtual void randompos(point &)=0;          //!< Random point within container
    virtual string info();                      //!< Return info string
    virtual string povray();                    //!< POVRAY object representing the cell
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
    string povray();
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
    inline bool collision(point &p) {};
};

/*! \brief "Clutch" like container.
 *  \author Mikael Lund
 *
 *  A spherical cell with a particle inaccessible area shaped
 *  as a disc in the middle of the sphere. The disc is parallel
 *  to the XY-plane and spans two Z-values as specified in the
 *  constructor.
 *
 *  \image html clutch.png
 */
class clutch : public container {
  private:
    double r2;
    float diameter;
  public:
    float r,zmin,zmax;
    clutch(float, float, float);
    void randompos(point &);
    bool collision(point &p) {
      if (p.z<zmax && p.z>zmin)
        return true;
      if (p.x*p.x+p.y*p.y+p.z*p.z > r2)
        return true;
      return false;
    }
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
