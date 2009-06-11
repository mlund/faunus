#ifndef FAU_CONTAINER_H
#define FAU_CONTAINER_H

#include "faunus/particles.h"
#include "faunus/slump.h"
#include "faunus/species.h"
#include "faunus/inputfile.h"

namespace Faunus {
  /*! \brief Polymorphic class for simulation containers
   *  \author Mikael Lund
   */
  class container : public particles {
    protected:
      slump slp;
      double volume;                                      //!< Volume of the container [AA^3]
    public:
      //atoms atom;
      double getvolume() {return volume;}
//      container();
//      container(inputfile &);
      virtual void setvolume(double){}                    //!< Specify new volume
      virtual bool collision(const particle &)=0;            //!< Check for collision with walls
      virtual void randompos(point &)=0;                  //!< Random point within container
      virtual string info();                              //!< Return info string
      virtual string povray();                            //!< POVRAY object representing the cell
      virtual void boundary(point &){};                   //!< Apply boundary conditions to a point
      virtual void scale(point&, const double&) const {}  //!< Scale point to a new box length
      inline virtual double sqdist(const point &a, const point &b) { 
        return a.sqdist(b); }
      inline virtual double dist(const point &a,const point &b) {//!< Calculate distance between points
        return a.dist(b);
      }
      inline virtual point vdist(const point &a, const point &b) { return a-b; }
  };

  /*! \brief Spherical simulation container
   *  \author Mikael Lund
   */
  class cell : public container {
    private:
      double r2,diameter;
      void setradius(double);
    public:
      double r;              //!< Radius
      cell(double);
      cell(inputfile &);
      string info();
      void randompos(point &);
      string povray();
      bool collision(const particle &a) {
        double x,y,z;
        x=std::abs(a.x)+a.radius;
        y=std::abs(a.y)+a.radius;
        z=std::abs(a.z)+a.radius;
        return ( x*x+y*y+z*z > r2 ) ? true:false;
      }
  };

  //---------------------------------------------------------
  /*! \brief Cubic simulation container w. periodic boundaries
   *  
   *  \author Mikael Lund
   */
  class box : public container {
    private:
      bool setlen(double);                 //!< Reset box length
    public:
      double len_half, len_inv, tlen_inv;  //!< tlen is the trial box
      void setvolume(double);
      double len;                          //!< Side length
      box(double);
      box(inputfile &);
      string info();
      void randompos(point &);
      //void randompos(vector<point> &);   // not implemented
      point randompos();
      bool collision(const particle &a) {
        if (std::abs(a.x)>len_half ||
            std::abs(a.y)>len_half ||
            std::abs(a.z)>len_half ) {
		return true;
	}
        return false;
      }
      bool clash(const particle &a, const particle &b) {
        point p;
        p.x=std::abs(a.x-b.x);
        p.y=std::abs(a.y-b.y);
        p.z=std::abs(a.z-b.z);
        if (p.x>len_half) p.x-=len;
        if (p.y>len_half) p.y-=len;
        if (p.z>len_half) p.z-=len;
        return (pow(p.len(),2)<pow(a.radius+b.radius, 2))
          ? true : false;
      }
      string povray();
      //! Calculate distance using the minimum image convention
      inline double dist(const point &a, const point &b) { return a.dist(b, len, len_half); }
      inline double sqdist(const point &a, const point &b) { return a.sqdist(b, len, len_half); }
      inline point vdist(const point &a, const point &b) {
        point r;
        r.x=std::abs(a.x-b.x);
        r.y=std::abs(a.y-b.y);
        r.z=std::abs(a.z-b.z);
        if (r.x>len_half) r.x-=len;
        if (r.y>len_half) r.y-=len;
        if (r.z>len_half) r.z-=len;
        return r;
      }
      inline int anint(double x) { return int(x>0. ? x+.5 : x-.5); }
      //! Apply periodic boundary conditions
      inline void boundary(point &a) {
        if (std::abs(a.x)>len_half) a.x-=len*anint(a.x/len);
        if (std::abs(a.y)>len_half) a.y-=len*anint(a.y/len);
        if (std::abs(a.z)>len_half) a.z-=len*anint(a.z/len);
      }
      inline void scale(point &a, const double &newlen) const { a = a*(newlen/len); }
  };

  //---------------------------------------------------------
  /*! \brief Box with periodic boundaries in the x and y direction.
   *  
   *  \author Mikael Lund
   *  \date Asljunga, 2008
   */
  class slit : public box {
    public:
      slit(inputfile &);
      string info();
      inline void boundary(point &a) {
        a.x=a.x-len*anint(a.x*len_inv);
        a.y=a.y-len*anint(a.y*len_inv);
      }
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
      double diameter;
    public:
      double r,zmin,zmax;
      clutch(double, double, double);
      void randompos(point &);
      bool collision(const particle &a) {
        if (a.z<zmax && a.z>zmin)
          return true;
        if (a.x*a.x+a.y*a.y+a.z*a.z > r2)
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
      double len;   //!< Cylinder length
      double r;     //!< Cylinder radius
      double r2;    //!< Cylinder radius squared
      double diameter;
      cylinder(double,double);
      void randompos(point &);
      bool collision(const particle &a) {
        return 
          (a.x*a.x+a.y*a.y>r2 || (a.z<0||a.z>len)) ? true:false;
      };
      string info(); //!< Cylinder info
      string povray();
  };
}
#endif
