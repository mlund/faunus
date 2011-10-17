#ifndef FAU_GEOMETRY_H
#define FAU_GEOMETRY_H

#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/slump.h"

namespace Faunus {

  class InputMap;
  class group;

  /*!
   * \brief Contain geometry related classes for describing simulation cells
   */
  namespace Geometry {

    /*!
     * \brief Polymorphic class for simulation geometries
     * \author Mikael Lund
     *
     * This is an important polymorph class that handles simulation geometries.
     * It contains distance calculation functions, boundary conditions, volume
     * etc.
     *
     * \note This class uses dynamic polymorphism, i.e. virtual functions, which
     *       may have negative impact on performance as function inlining may not be
     *       possible. This is usually a problem only for inner loop distance calculations.
     *       To get optimum performance in inner loops use a derived class directly and do
     *       static compile-time polymorphism.
     * \todo Migrate to NVI (Nonvirtual Interface) pattern -- see "C++ Coding Standards", p.68
     */
    class Geometrybase {
      private:
        virtual string _info(char)=0;
      protected:
        slump slp;
        double volume;                                      //!< Volume of the container [AA^3]
        string name;                                        //!< Name of the geometry
      public:
        enum collisiontype {BOUNDARY,ZONE};                 //!< Types for collision() function
        double getvolume() const;                           //!< Get volume of container
        virtual void setvolume(double);                     //!< Specify new volume
        virtual bool collision(const particle&, collisiontype=BOUNDARY)=0;//!< Check for collision with boundaries, forbidden zones, matter,..
        virtual void randompos(point &)=0;                  //!< Random point within container
        virtual void boundary(point &) const=0;             //!< Apply boundary conditions to a point
        virtual void scale(point&, const double&) const;    //!< Scale point to a new volume - for NPT ensemble
        double dist(const point&,const point&);             //!< Distance between two points
        virtual double sqdist(const point &a, const point &b)=0;
        virtual point vdist(const point&, const point&)=0;
        string info(char=20);                               //!< Return info string
        bool save(string);                                  //!< Save container state to disk
        bool load(string,bool=false);                       //!< Load container state from disk
        virtual ~Geometrybase();
    };

    /*! \brief Spherical geometry
     *  \author Mikael Lund
     *  \todo Implement space scaling for isobaric ensemble
     *
     *  This is a spherical simulation container, surrounded by a hard wall.
     */
    class Sphere : public Geometrybase {
      private:
        double r2,diameter;
      public:
        void setradius(double);
        double r;              //!< Radius
        Sphere(double);
        Sphere(InputMap&);
        string _info(char);
        void setvolume(double);
        void randompos(point &);
        void boundary(point &) const;
        bool collision(const particle &, collisiontype=BOUNDARY);
        inline double sqdist(const point &a, const point &b) {
          register double dx,dy,dz;
          dx=a.x-b.x;
          dy=a.y-b.y;
          dz=a.z-b.z;
          return dx*dx + dy*dy + dz*dz;
        }
        point vdist(const point&, const point&);
    };

    //---------------------------------------------------------
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
      protected:
        bool setslice(point, point);             //!< Reset slice position
        point len_inv;                           //!< Inverse sidelengths

      public:
        Cuboid(InputMap&);                       //!< Read input parameters
        bool setlen(point);                      //!< Reset Cuboid sidelengths
        point len;                               //!< Sidelengths
        point len_half;                          //!< Half sidelength
        point slice_min, slice_max;              //!< Position of slice corners
        point randompos();                       //!< Get point with random position
        void randompos(point &);                 //!< Move point to random position
        bool save(string);                       //!< Save container state to disk
        bool load(string,bool=false);            //!< Load container state from disk
        bool collision(const particle&, collisiontype=BOUNDARY);

        inline double sqdist(const point &a, const point &b) {
          double dx=std::abs(a.x-b.x);
          double dy=std::abs(a.y-b.y);
          double dz=std::abs(a.z-b.z);
          if (dx>len_half.x) dx-=len.x;
          if (dy>len_half.y) dy-=len.y;
          if (dz>len_half.z) dz-=len.z;
          return dx*dx + dy*dy + dz*dz;
        }

        inline point vdist(const point &a, const point &b) {       //!< Distance vector
          point r=a-b;
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

        inline int anint(double x) const {
          return int(x>0. ? x+.5 : x-.5);
        }

        //! Apply periodic boundary conditions
        inline void boundary(point &a) const {
          if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
          if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
          if (std::abs(a.z)>len_half.z) a.z-=len.z*anint(a.z*len_inv.z);
        }
        virtual void scale(point&, const double&) const; //!< Scale point to a new volume - for NPT ensemble
    };

    /*!
     * \brief Cuboidslit: cubuid without periodic boundary in the z direction
     * \author Chris Evers
     * \date Lund, nov 2010
     */
    class Cuboidslit : public Cuboid {
      public:
        Cuboidslit(InputMap &);

        //! Calculate distance using the minimum image convention
        inline double sqdist(const point &a, const point &b) {   //!< Squared distance 
          double dx=std::abs(a.x-b.x);
          double dy=std::abs(a.y-b.y);
          double dz=a.z-b.z;
          if (dx>len_half.x) dx-=len.x;
          if (dy>len_half.y) dy-=len.y;                                      
          return dx*dx + dy*dy + dz*dz;
        }   

        inline point vdist(const point &a, const point &b) {       //!< Distance vector
          point r=a-b;
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

        //! Apply periodic boundary conditions
        inline void boundary(point &a) const {
          if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
          if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
        }
    };

    /*! \brief Cylindrical simulation container
     *  \author Mikael Lund & Bjoern Persson
     *  \todo Needs some testing
     *
     *  This is a Cylinder container where all walls
     *  are HARD. The origin is in the middle of the
     *  Cylinder.
     */
    class Cylinder : public Geometrybase {
      private:
        string _info(char); //!< Cylinder info
        double halflen;
        double r2;    //!< Cylinder radius squared
        void init(double,double);
      public:
        double len;   //!< Cylinder length
        double r;     //!< Cylinder radius
        double diameter;
        Cylinder(double, double);
        Cylinder(InputMap &);
        void randompos(point &);
        bool collision(const particle&, collisiontype=BOUNDARY);
        inline double sqdist(const point &a, const point &b) {
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

    class GeometryList {
      private:
        vector<Geometrybase*> base;
      public:
        void add(Geometrybase &g) {
          if (&g!=NULL)
            if ( std::find(base.begin(), base.end(), &g)==base.end() )
              base.push_back(&g);
        }

        template<typename Tgeo> void deepCopy(const Geometry::Geometrybase& source) {
          assert(&source!=NULL);
          const Tgeo* src=dynamic_cast<const Tgeo*>(&source); //get dereved class pointer to source
          for (auto *g : base) {
            Tgeo* target=dynamic_cast<Tgeo*>(g); // derived class pointer to target
            (*target)=(*src);                    // copy data from source to target
          }
        }

    };
  }//namespace Geometry

}//namespace
#endif
