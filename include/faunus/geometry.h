#ifndef FAU_GEOMETRY_H
#define FAU_GEOMETRY_H

#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/slump.h"

namespace Faunus {

  class InputMap;

  /*!
   * \brief Namespace for geometric operations.
   */
  namespace Geometry {

    /*!
     * \brief Polymorphic class for simulation geometries
     * \author Mikael Lund
     *
     * This is a polymorph class that handles simulation geometries.
     * It contains distance calculation functions, boundary conditions, volume
     * etc. A number of functions are defined as pure virtual (=0) meaning
     * that these MUST be defined in derived classes.
     *
     * \note This class uses dynamic polymorphism, i.e. virtual functions, which
     *       may have negative impact on performance as function inlining may not be
     *       possible. This is usually a problem only for inner loop distance calculations.
     *       To get optimum performance in inner loops use a derived class directly and do
     *       static, compile-time polymorphism (templates).
     */
    class Geometrybase {
      private:
        virtual string _info(char)=0;
        virtual void _setVolume(double)=0;
        virtual double _getVolume() const=0;
      protected:
        slump slp;
        string name;                                        //!< Name of the geometry
      public:
        enum collisiontype {BOUNDARY,ZONE};                 //!< Types for collision() function
        double getVolume() const;                           //!< Get volume of container (A^3)
        void setVolume(double);                             //!< Specify new volume (A^3)
        double dist(const Point&,const Point&);             //!< Distance between two points (A)
        string info(char=20);                               //!< Return info string
        bool save(string);                                  //!< Save geometry state to disk
        bool load(string,bool=false);                       //!< Load geometry state from disk

        virtual bool collision(const particle&, collisiontype=BOUNDARY)=0;//!< Check for collision with boundaries, forbidden zones, matter,..
        virtual void randompos(Point &)=0;                  //!< Random point within container
        virtual void boundary(Point &) const=0;             //!< Apply boundary conditions to a point
        virtual void scale(Point&, const double&) const;    //!< Scale point to a new volume - for NPT ensemble
        virtual double sqdist(const Point &a, const Point &b) const=0;
        virtual Point vdist(const Point&, const Point&)=0;  //!< Distance in vector form
        virtual ~Geometrybase();
    };

    /*!
     * \brief Spherical geometry
     * \author Mikael Lund
     *
     * This is a spherical simulation container, surrounded by a hard wall.
     */
    class Sphere : public Geometrybase {
      private:
        double r,r2,diameter;
        void _setVolume(double);
        double _getVolume() const;
        string _info(char);
      public:
        void setradius(double);                 //!< Set radius (angstrom)
        Sphere(double);                         //!< Construct from radius (angstrom)
        Sphere(InputMap&, string="Sphere");     //!< Construct from InputMap key \c prefix_radius
        void randompos(Point &);
        void boundary(Point &) const;
        bool collision(const particle &, collisiontype=BOUNDARY);
        inline double sqdist(const Point &a, const Point &b) const {
          register double dx,dy,dz;
          dx=a.x-b.x;
          dy=a.y-b.y;
          dz=a.z-b.z;
          return dx*dx + dy*dy + dz*dz;
        }
        Point vdist(const Point&, const Point&);
        void scale(Point&, const double&) const; //!< Linear scaling along radius (NPT ensemble)
    };

    /*! \brief Cuboid geometry with periodic boundaries
     *
     *  \author Chris Evers
     *  \date Lund, nov 2010
     *
     *  The Cuboid simulation container has right angles, rectangular faces 
     *  and periodic boundaries. A slice can be introduced to constrain the position
     *  of some of the space to a part of the Cuboid. The function slicecollision
     *  can be used to make sure space are positioned within in the slice.
     */
    class Cuboid : public Geometrybase {
      private:
        string _info(char);                      //!< Return info string
        void _setVolume(double);
        double _getVolume() const;
      protected:
        bool setslice(Point, Point);             //!< Reset slice position
        Point len_inv;                           //!< Inverse sidelengths
        inline int anint(double x) const {
          return int(x>0. ? x+.5 : x-.5);
        }
        Point slice_min, slice_max;              //!< Position of slice corners
 
      public:
        Cuboid(InputMap&);                       //!< Read input parameters
        bool setlen(Point);                      //!< Reset Cuboid sidelengths
        Point len;                               //!< Sidelengths
        Point len_half;                          //!< Half sidelength
        Point randompos();           
        void randompos(Point&);      
        bool save(string);           
        bool load(string,bool=false);
        bool collision(const particle&, collisiontype=BOUNDARY);

        inline double sqdist(const Point &a, const Point &b) const {
          double dx=std::abs(a.x-b.x);
          double dy=std::abs(a.y-b.y);
          double dz=std::abs(a.z-b.z);
          if (dx>len_half.x) dx-=len.x;
          if (dy>len_half.y) dy-=len.y;
          if (dz>len_half.z) dz-=len.z;
          return dx*dx + dy*dy + dz*dz;
        }

        inline Point vdist(const Point &a, const Point &b) {
          Point r=a-b;
          if (r.x>len_half.x)
            r.x-=len.x;
          else if (r.x<-len_half.x)
            r.x+=len.x;
          if (r.y>len_half.y)
            r.y-=len.y;
          else if (r.y<-len_half.y)
            r.y+=len.y;
          if (r.z>len_half.z)
            r.z-=len.z;
          else if (r.z<-len_half.z)
            r.z+=len.z;
          return r;
        }

        inline void boundary(Point &a) const {
          if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
          if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
          if (std::abs(a.z)>len_half.z) a.z-=len.z*anint(a.z*len_inv.z);
        }

        void scale(Point&, const double&) const;
    };

    /*!
     * \brief Cubuid with no periodic boundaries in z direction
     * \author Chris Evers
     * \date Lund, nov 2010
     */
    class Cuboidslit : public Cuboid {
      public:
        Cuboidslit(InputMap &);

        inline double sqdist(const Point &a, const Point &b) const {
          double dx=std::abs(a.x-b.x);
          double dy=std::abs(a.y-b.y);
          double dz=a.z-b.z;
          if (dx>len_half.x) dx-=len.x;
          if (dy>len_half.y) dy-=len.y;                                      
          return dx*dx + dy*dy + dz*dz;
        }   

        inline Point vdist(const Point &a, const Point &b) {
          Point r=a-b;
          if (r.x>len_half.x)
            r.x-=len.x;
          else if (r.x<-len_half.x)
            r.x+=len.x;
          if (r.y>len_half.y)
            r.y-=len.y;
          else if (r.y<-len_half.y)
            r.y+=len.y;
          return r;
        }

        inline void boundary(Point &a) const {
          if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
          if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
        }
    };

    /*!
     * \brief Cylindrical simulation container
     * \author Mikael Lund & Bjorn Persson
     * \todo Periodic ends seems more useful - any particular reason for having hard ends?
     *
     * This is a cylindrical simulation container where all walls
     * are HARD. The origin is in the middle of the cylinder.
     */
    class Cylinder : public Geometrybase {
      private:
        string _info(char); //!< Cylinder info
        void _setVolume(double);
        double _getVolume() const;
        double halflen;
        double r2;    //!< Cylinder radius squared
        void init(double,double);
        double len;   //!< Cylinder length
        double r;     //!< Cylinder radius
        double diameter;
      public:
        Cylinder(double, double);
        Cylinder(InputMap &);
        void randompos(Point &);
        bool collision(const particle&, collisiontype=BOUNDARY);
        inline double sqdist(const Point &a, const Point &b) const {
          register double dx,dy,dz;
          dx=a.x-b.x;
          dy=a.y-b.y;
          dz=a.z-b.z;
          return dx*dx + dy*dy + dz*dz;
        }
    };

#ifdef HYPERSPHERE
    /*! \brief HyperSphere simulation container
     *  \author Martin Trulsson
     *  \date Lund, 2009
     */
    class hyperSphere : public Sphere {
      private:
        const double pi;
        string _info(char);
      public:
        hyperSphere(InputMap&);
        void randompos(point &);
        bool collision(const particle &, collisiontype=BOUNDARY);

        double dist(const point &a, const point &b) {
          return r*a.geodesic(b); // CHECK!!! (virtual=slow!)
        }

        inline double sqdist(const point &a, const point &b) const {
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

    /*!
     * \brief Vector rotation routines
     * \note Boundary condition are respected.
     * \author Mikael Lund
     * \date Canberra, 2009
     */
    class VectorRotate {
      private:
        Point origin, u;
        double cosang, sinang;
        double e1mcox, e1mcoy, e1mcoz;
        Geometrybase* geoPtr;
      public:
        void setAxis(Geometrybase&, const Point&, const Point&, double);  //!< Set axis of rotation and degrees to rotate
        Point rotate(const Geometrybase&, Point) const;                   //!< Rotate point around axis
    };

  }//namespace Geometry
}//namespace Faunus
#endif
