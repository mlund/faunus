#ifndef FAU_CONTAINER_H
#define FAU_CONTAINER_H

#include "faunus/particles.h"
#include "faunus/point.h"
#include "faunus/slump.h"

namespace Faunus {
  class inputfile;
 
  /*! \brief Polymorphic class for simulation containers
   *  \author Mikael Lund
   */
  class container : public particles {
    protected:
      slump slp;
      double volume;                                      //!< Volume of the container [AA^3]
    public:
      double getvolume() {return volume;}
      virtual void setvolume(double);                     //!< Specify new volume
      virtual bool collision(const particle &)=0;         //!< Check for collision with walls
      virtual bool slicecollision(const particle &) { return false; };      //!< Check for collision with the slice borders 
      virtual void randompos(point &)=0;                  //!< Random point within container
      virtual string info();                              //!< Return info string
      virtual string povray();                            //!< POVRAY object representing the cell
      virtual void boundary(point &) const {};            //!< Apply boundary conditions to a point
      virtual void scale(point&, const double&) const {}  //!< Scale point to a new volume length
      virtual bool saveToDisk(string);                    //!< Save container date to disk
      virtual bool loadFromDisk(string,bool=false);       //!< Load container data from disk
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
    public:
      void setradius(double);
      double r;              //!< Radius
      cell(double);
      cell(inputfile &);
      string info();
      void setvolume(double);
      void randompos(point &);
      string povray();
      bool collision(const particle &a) {
        double x,y,z;
        x=std::abs(a.x);//+a.radius;
        y=std::abs(a.y);//+a.radius;
        z=std::abs(a.z);//+a.radius;
        return ( x*x+y*y+z*z > r2 ) ? true:false;
      }
  };

  //---------------------------------------------------------
  /*! \brief Cuboid simulation container with periodic boundaries
   *
   *  \author Chris Evers
   *  \date Lund, nov 2010
   *
   *  The cuboid simulation container has right angles, rectangular faces 
   *  and periodic boundaries. A slice can be introduced to constrain the position
   *  of some of the particles to a part of the cuboid. The function slicecollision
   *  can be used to make sure particles are positioned within in the slice.
   */
  class cuboid : public container {
    private:
      bool setlen(point);                      //!< Reset cuboid sidelengths
      bool setslice(point, point);             //!< Reset slice position

    public:
      cuboid(inputfile &);                    //!< Read input parameters

      point len;                               //!< Sidelengths
      point len_half, len_inv;                 //!< Half and inverse sidelengths
      point slice_min, slice_max;              //!< Position of slice corners

      string info();                           //!< Return info string

      point randompos();                       //!< Get point with random position
      void randompos(point &);                 //!< Move point to random position
//    //void randompos(vector<point> &);       //   not implemented

      bool slicecollision(const particle &a);  //!< Check collision with slice edges
      bool collision(const particle &a) {      //!< Check collision with cuboid edges
        if (std::abs(a.x) > len_half.x ||
            std::abs(a.y) > len_half.y ||
            std::abs(a.z) > len_half.z  )
          return true;
        return false;
      }

      //! Calculate distance using the minimum image convention
      inline double sqdist(const point &p1, const point &p2) {   //!< Squared distance 
        double dx=std::abs(p1.x-p2.x),
               dy=std::abs(p1.y-p2.y),
               dz=std::abs(p1.z-p2.z);
        if (dx>len_half.x) dx-=len.x;
        if (dy>len_half.y) dy-=len.y;
        if (dz>len_half.z) dz-=len.z;
        return dx*dx + dy*dy + dz*dz;
      }

      inline double dist(const point &p1, const point &p2) {     //!< Root squared distance
        return sqrt(sqdist(p1,p2));
      }

      inline point vdist(const point &a, const point &b) {       //!< Distance vector
        point r;
        r.x=std::abs(a.x-b.x);
        r.y=std::abs(a.y-b.y);
        r.z=std::abs(a.z-b.z);
        if (r.x>len_half.x) r.x-=len.x;
        if (r.y>len_half.y) r.y-=len.y;
        if (r.z>len_half.z) r.z-=len.z;
        return r;
      }

      inline int anint(double x) const {
        return int(x>0. ? x+.5 : x-.5);
      }

      //! Apply periodic boundary conditions
      inline void boundary(point &a) const {
        if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
        if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
        if (std::abs(a.z)>len_half.z) a.z-=len.z*anint(a.z*len_inv.z);
      }
  };

  //---------------------------------------------------------
  /*! \brief Cuboidslit: cubuid without periodic boundary in the z direction
   *  
   *  \author Chris Evers
   *  \date Lund, nov 2010
   */
  class cuboidslit : public cuboid {
    public:
      cuboidslit(inputfile &);
      string info();

      //! Calculate distance using the minimum image convention
      inline double sqdist(const point &p1, const point &p2) {   //!< Squared distance 
        double dx=std::abs(p1.x-p2.x),
               dy=std::abs(p1.y-p2.y),
               dz=p1.z-p2.z;
        if (dx>len_half.x) dx-=len.x;
        if (dy>len_half.y) dy-=len.y;
        return dx*dx + dy*dy + dz*dz;
      }

      inline point vdist(const point &a, const point &b) {       //!< Distance vector
        point r;
        r.x=std::abs(a.x-b.x);
        r.y=std::abs(a.y-b.y);
        r.z=a.z-b.z;
        if (r.x>len_half.x) r.x-=len.x;
        if (r.y>len_half.y) r.y-=len.y;
        return r;
      }

      //! Apply periodic boundary conditions
      inline void boundary(point &a) const {
        if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
        if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
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
      void setvolume(double);              //!< Set volume (and sidelength) of box
      double len;                          //!< Side length
      box(double);
      box(inputfile &);
      string info();
      string povray();
      void randompos(point &);
      //void randompos(vector<point> &);   // not implemented
      point randompos();

      bool collision(const particle &a) {
        if (std::abs(a.x)>len_half)
          return true;
        if (std::abs(a.y)>len_half)
          return true;
        if (std::abs(a.z)>len_half)
          return true;
        return false;
      }

      bool clash(const particle &a, const particle &b) {
        point c;
        c.x=std::abs(a.x-b.x);
        c.y=std::abs(a.y-b.y);
        c.z=std::abs(a.z-b.z);
        if (c.x>len_half) c.x-=len;
        if (c.y>len_half) c.y-=len;
        if (c.z>len_half) c.z-=len;
        return (pow(c.len(),2)<pow(a.radius+b.radius, 2))
          ? true : false;
      }

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

      inline int anint(double x) const {
        return int(x>0. ? x+.5 : x-.5);
      }

      /*!
       * \brief Apply periodic boundary conditions
       */
      inline void boundary(point &a) const {
        if (std::abs(a.x)>len_half) a.x-=len*anint(a.x/len);
        if (std::abs(a.y)>len_half) a.y-=len*anint(a.y/len);
        if (std::abs(a.z)>len_half) a.z-=len*anint(a.z/len);
      }

      /*!
       * \brief Linear scaling of particle position after a volume fluctuation
       */
      inline void scale(point &a, const double &newlen) const {
        a = a*(newlen/len);
      }
  };

  //---------------------------------------------------------
  /*! \brief Box with periodic boundaries in the x and y direction.
   *  
   *  \author Mikael Lund
   *  \date Asljunga, 2008
   */
  class slit : public box {
    //private 
      //bool setlen(double);                 //!< Reset box length
    public:
      double xyarea, zlen, zlen_half; //crossecional area 
      slit(inputfile &);
      string info();
      void randompos(point &);

      inline void boundary(point &a) const {
        a.x=a.x-len*anint(a.x*len_inv);
        a.y=a.y-len*anint(a.y*len_inv);
      }

      inline double sqdist(const point &p1, const point &p2) {
        double dz=p1.z-p2.z,
               dx=std::abs(p1.x-p2.x),
               dy=std::abs(p1.y-p2.y);
        if (dx>len_half) dx-=len;
        if (dy>len_half) dy-=len;
        return dx*dx + dy*dy + dz*dz;
      }

      inline double dist(const point &p1, const point &p2) {
        return sqrt(sqdist(p1,p2));
      }
      
      inline point vdist(const point &a, const point &b) {
        point r;
        r.x=std::abs(a.x-b.x);
        r.y=std::abs(a.y-b.y);
        r.z=a.z-b.z;
        if (r.x>len_half) r.x-=len;
        if (r.y>len_half) r.y-=len;
        return r;
      }

      bool collision(const particle &a) {
        if (std::abs(a.x)>len_half ||
            std::abs(a.y)>len_half ||
            std::abs(a.z)>zlen_half ) {
		return true;
	}
        return false;
      }
  };

  class xyplane : public box {
    public:
      xyplane(inputfile &in);
      void randompos(point &);
      void randompos(vector<point> &);
  };

  
  /*! \brief Gridded slit box.
   *  
   *  \author Mikael Lund
   *  \date Asljunga, 2010
   */
  class gridslit : public slit {
  private:
    double l;     //!< Grid spacing
  public:
    double ngrid; //!< Number of grids in each dimension
    gridslit(inputfile &);
    string info();
    bool collision(const particle &);
    void randompos(point &);
    inline void boundary(point &a) const {
      slit::boundary(a);
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

#ifdef HYPERSPHERE
  /*! \brief Hypersphere simulation container
   *  \author Martin Trulsson
   *  \date Lund, 2009
   */
  class hypersphere : public cell {
    private:
      static const double pi;
    public:
      hypersphere(inputfile &);
      string info();
      void randompos(point &);
      bool collision(const particle &);

      inline double dist(const point &a, const point &b) {
        return r*a.geodesic(b); // CHECK!!! (virtual=slow!)
      }

      inline double sqdist(const point &a, const point &b) {
        return pow(dist(a,b),2); // !! SHOULD BE REAL DISTANCE CHECK!! (virtual=slow!)
      }

      bool overlap(const particle &a) {
        for (int i=0; i<p.size(); i++)
          if (hyperoverlap(a,p[i])==true)
            return true;
        return false;
      }

      inline bool hyperoverlap(const particle &a, const particle &b) {
        return (r*a.geodesic(b)<a.radius+b.radius);
      }
  };
#endif
}//namespace
#endif
