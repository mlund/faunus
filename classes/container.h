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
    float volume;                                  //!< Volume of the container [AA^3]
    inline virtual bool collision(point &)=0;      //!< Check for collision with walls
    virtual void randompos(point &)=0;             //!< Random point within container
    virtual string info();                         //!< Return info string
    virtual string povray();                       //!< POVRAY object representing the cell
    inline virtual void boundary(point &)=0;       //!< Apply boundary conditions to a point
    inline virtual double dist(point &a,point &b) {//!< Calculate distance between points
      return a.dist(b);
    }
    //void displace(point &, float);      //!< Displace particle
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
    inline void boundary(point &) {};
    void randompos(point &);
    string povray();
    inline bool collision(point &p) {
      return 
        (p.x*p.x+p.y*p.y+p.z*p.z > r2) ? true:false;
    }
};

//---------------------------------------------------------
/*! \brief Cubic simulation container w. periodic boundaries
 *  \author Mikael Lund
 *  \todo Not finished!
 */
class box : public container {
  private:
    point d;
    double len_half, len_inv;
  public:
    double len; //!< Side length
    box(double);
    string info();
    void randompos(point &);
    void randompos(vector<point> &);
    point randompos();
    inline bool collision(point &p) {return false;};
    string povray();

    //! Calculate distance using the minimum image convention
    inline double dist(point &a, point &b) {
      return a.dist(b, len, len_inv);
    }

    //! Apply periodic boundary conditions
    inline void boundary(point &p) {
      p.x=p.x-len*floor(p.x*len_inv+.5);
      p.y=p.y-len*floor(p.y*len_inv+.5);
      p.z=p.z-len*floor(p.z*len_inv+.5);
    }

    //! Randomly displace particle w. bpc.
    //void displace(point &p, float dp) {
    //  p.x+=dp*slp.random_half();
    //  p.y+=dp*slp.random_half();
    //  p.z+=dp*slp.random_half();
    //  bpc(p);
    //}
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
    void boundary(point &) {};
    bool collision(point &p) {
      if (p.z<zmax && p.z>zmin)
        return true;
      if (p.x*p.x+p.y*p.y+p.z*p.z > r2)
        return true;
      return false;
    }
};

/*! \brief Cylindrical simulation container
 *  \author Mikael Lund/Bjoern Persson
 *  \todo Needs some testing
 */
class cylinder : public container {
  public:
    float len;   //!< Cylinder length
    float r;     //!< Cylinder radius
    float r2;    //!< Cylinder radius squared
    float diameter;
    void boundary(point &) {};
    cylinder(float,float);
    void randompos(point &);
    inline bool collision(point &p) {
      return 
        (p.x*p.x+p.y*p.y>r2 || (p.z<0||p.z>len)) ? true:false;
    };
    string info(); //!< Cylinder info
    string povray();
};
#endif
