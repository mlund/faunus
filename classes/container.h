#ifndef _container_h
#define _contained_h

#include "slump.h"
#include "point.h"

/*! \brief Polymorphic class for simulation containers
 *  \author Mikael Lund
 */
class container {
  protected:
    slump slp;
  public:
    float volume;                               //!< Volume of the container [AA^3]
    virtual void randompos(point &)=0;          //!< Random point within container
    inline virtual bool collision(point &)=0;   //!< Check for collision with walls
};

/*! \brief Spherical simulation container
 *  \author Mikael Lund
 */
class cell : public container {
  public:
    float r,              //!< Radius
          r2,             //!< Squared radius
          diameter;       //!< Diameter
    cell(float);
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
    bool collision(point &p) {};
};

#endif
