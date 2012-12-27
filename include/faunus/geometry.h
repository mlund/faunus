#ifndef FAU_GEOMETRY_H
#define FAU_GEOMETRY_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/slump.h>
#include <Eigen/Eigen>
#endif

namespace Faunus {

  /*!
   * \brief Namespace for geometric operations.
   *
   * This namespace contains classes for handling various simulation geometries
   * such as cubes, spheres, cylinder, slits etc. The geometry of a simulation
   * is handled by the base class Geometry::Geometrybase that gives all geometries
   * a common interfaces to handle distance calculations, boundary conditions,
   * volume calculation and so forth. In this way, the geometry of a simulation may be
   * changed without any code changes in for example energy calculations.
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
        inline int anint(double x) const {
          return int(x>0. ? x+.5 : x-.5);
        }
      public:
        enum collisiontype {BOUNDARY,ZONE};                 //!< Types for collision() function
        double getVolume() const;                           //!< Get volume of container (A^3)
        void setVolume(double);                             //!< Specify new volume (A^3)
        double dist(const Point&,const Point&);             //!< Distance between two points (A)
        string info(char=20);                               //!< Return info string
        bool save(string);                                  //!< Save geometry state to disk
        bool load(string,bool=false);                       //!< Load geometry state from disk

        virtual bool collision(const particle&, collisiontype=BOUNDARY) const=0;//!< Check for collision with boundaries, forbidden zones, matter,..
        virtual void randompos(Point &)=0;              //!< Random point within container
        virtual void boundary(Point &) const=0;             //!< Apply boundary conditions to a point
        virtual void scale(Point&, const double&) const;    //!< Scale point to a new volume - for NPT ensemble
        virtual double sqdist(const Point &a, const Point &b) const=0; //!< Squared distance between two points
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
        Sphere(InputMap&, string="sphere");     //!< Construct from InputMap key \c prefix_radius
        void randompos(Point &);
        void boundary(Point &) const;
        bool collision(const particle &, collisiontype=BOUNDARY) const;
        inline double sqdist(const Point &a, const Point &b) const {
          return (a-b).squaredNorm();
        }
        inline Point vdist(const Point &a, const Point &b) { return a-b; }
        void scale(Point&, const double&) const; //!< Linear scaling along radius (NPT ensemble)
    };

    /*! \brief Cuboid geometry with periodic boundaries
     *
     *  The Cuboid simulation container has right angles, rectangular faces 
     *  and periodic boundaries.
     *
     *  \author Chris Evers
     *  \date Lund, nov 2010
     */
    class Cuboid : public Geometrybase {
      private:
        string _info(char);                      //!< Return info string
        void _setVolume(double);
        double _getVolume() const;
        enum scaletype {XYZ,XY};
        scaletype scaledir;                      //!< Scale directions for pressure scaling
      protected:
        Point len_inv;                           //!< Inverse sidelengths
      public:
        Cuboid(InputMap&);                       //!< Construct from input file
        bool setlen(Point);                      //!< Reset Cuboid sidelengths
        Point len;                               //!< Sidelengths
        Point len_half;                          //!< Half sidelength
        Point randompos();           
        void randompos(Point&);      
        bool save(string);           
        bool load(string,bool=false);
        bool collision(const particle&, collisiontype=BOUNDARY) const;

        inline double sqdist(const Point &a, const Point &b) const {
          double dx( std::abs(a.x()-b.x()) );
          double dy( std::abs(a.y()-b.y()) );
          double dz( std::abs(a.z()-b.z()) );
          if (dx>len_half.x()) dx-=len.x();
          if (dy>len_half.y()) dy-=len.y();
          if (dz>len_half.z()) dz-=len.z();
          return dx*dx + dy*dy + dz*dz;
        }

        inline Point vdist(const Point &a, const Point &b) {
          Point r=a-b;
          if (r.x()>len_half.x())
            r.x()-=len.x();
          else if (r.x()<-len_half.x())
            r.x()+=len.x();
          if (r.y()>len_half.y())
            r.y()-=len.y();
          else if (r.y()<-len_half.y())
            r.y()+=len.y();
          if (r.z()>len_half.z())
            r.z()-=len.z();
          else if (r.z()<-len_half.z())
            r.z()+=len.z();
          return r;
        }

        inline void boundary(Point &a) const {
          if (std::abs(a.x())>len_half.x()) a.x()-=len.x()*anint(a.x()*len_inv.x());
          if (std::abs(a.y())>len_half.y()) a.y()-=len.y()*anint(a.y()*len_inv.y());
          if (std::abs(a.z())>len_half.z()) a.z()-=len.z()*anint(a.z()*len_inv.z());
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
          double dx=std::abs(a.x()-b.x());
          double dy=std::abs(a.y()-b.y());
          double dz=a.z()-b.z();
          if (dx>len_half.x()) dx-=len.x();
          if (dy>len_half.y()) dy-=len.y();                                      
          return dx*dx + dy*dy + dz*dz;
        }   

        inline Point vdist(const Point &a, const Point &b) {
          Point r(a-b);
          if (r.x()>len_half.x())
            r.x()-=len.x();
          else if (r.x()<-len_half.x())
            r.x()+=len.x();
          if (r.y()>len_half.y())
            r.y()-=len.y();
          else if (r.y()<-len_half.y())
            r.y()+=len.y();
          return r;
        }

        inline void boundary(Point &a) const {
          if (std::abs(a.x())>len_half.x()) a.x()-=len.x()*anint(a.x()*len_inv.x());
          if (std::abs(a.y())>len_half.y()) a.y()-=len.y()*anint(a.y()*len_inv.y());
        }
    };

    /*!
     * \brief Cylindrical simulation container
     * \author Mikael Lund and Bjorn Persson
     *
     * This is a cylindrical simulation container where all walls
     * are HARD. The origin is in the middle of the cylinder.
     */
    class Cylinder : public Geometrybase {
      private:
        string _info(char); //!< Cylinder info
        void _setVolume(double);
        double _getVolume() const;
        void init(double,double);
        double r2;    //!< Cylinder radius squared
        double r;     //!< Cylinder radius
        double diameter;
      protected:
        double len;   //!< Cylinder length
        double halflen; //!< Cylinder half length
      public:
        Cylinder(double, double);      //!< Construct from length and radius
        Cylinder(InputMap &);          //!< Construct from inputmap
        void randompos(Point &);
        void boundary(Point &) const;
        bool collision(const particle&, collisiontype=BOUNDARY) const;
        inline double sqdist(const Point &a, const Point &b) const {
          return (a-b).squaredNorm();
        }
        inline Point vdist(const Point &a, const Point &b) { return a-b; }
    };

    /*!
     * \brief Cylinder with periodic boundaries in the z direction
     * \author Mikael Lund
     * \warning something seems rotten...
     */
    class PeriodicCylinder : public Cylinder {
      public:
        PeriodicCylinder(double, double);
        PeriodicCylinder(InputMap&);
        void boundary(Point&) const;
        inline double sqdist(const Point &a, const Point &b) const {
          double dx=a.x()-b.x();
          double dy=a.y()-b.y();
          double dz=std::abs(a.z()-b.z());
          if (dz>halflen)
            dz-=len;
          return dx*dx + dy*dy + dz*dz;
        }
        inline Point vdist(const Point &a, const Point &b) {
          Point r=a-b;
          if (r.z()>halflen)
            r.z()-=len;
          else if (r.z()<-halflen)
            r.z()+=len;
          return r;
        }
    };

#ifdef HYPERSPHERE
    /*  \brief HyperSphere simulation container
     *  \author Martin Trulsson
     *  \date Lund, 2009
     */
    class hyperSphere : public Sphere {
      private:
        const double pi;
        string _info(char);
      public:
        hyperSphere(InputMap&);
        void randompos(Point&);
        bool collision(const particle&, collisiontype=BOUNDARY) const;

        // dist() is not virtual...
        double dist(const Point &a, const Point &b) {
          return r*a.geodesic(b);
        }

        inline double sqdist(const Point &a, const Point &b) const FOVERRIDE {
          return pow(dist(a,b),2);
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
     * \brief Calculate mass center of a group
     */
    template<typename Tgeometry, typename Tgroup>
      Point massCenter(const Tgeometry &geo, const p_vec &p, const Tgroup &g) {
        if (g.empty())
          return Point(0,0,0);
        assert(!p.empty());
        assert(g.back() < (int)p.size());
        assert(&geo!=NULL);
        double sum=0;
        Point cm(0,0,0);
        Point o = p[ g.front()+(g.back()-g.front())*0.5 ];  // set origo to middle particle
        for (auto i : g) {
          Point t = p[i]-o;       // translate to origo
          geo.boundary(t);        // periodic boundary (if any)
          cm += t * p[i].mw;
          sum += p[i].mw;
        }
        if (sum<1e-6) sum=1;
        cm=cm/sum + o;
        geo.boundary(cm);
        return cm;
      }

    Point massCenter(const Geometrybase&, const p_vec&); //!< Calculate mass center of a particle vector


    void translate(const Geometrybase&, p_vec&, Point); //!< Translate a particle vector by a vector
    void cm2origo(const Geometrybase&, p_vec&); //!< Translate a particle vector so mass center is in (0,0,0)

    /*!
      \brief Geometric transform of a Point (rotation, translation...)
      */
    template<typename Ttransformer>
      void transform(const Geometrybase &geo, const Ttransformer &t, Point &x) {
        x=t*x;
        geo.boundary(x);
      }

    /*!
     * \brief Find an empty space for a particle vector in a space of other particles
     * \author Mikael Lund
     * \date Asljunga 2011
     * \todo Implement random rotation in addition to current translation scheme.
     *
     */
    class FindSpace {
      private:
        bool containerOverlap(const Geometrybase&, const p_vec&);
        virtual bool matterOverlap(const Geometrybase&, const p_vec&, const p_vec&);
      public:
        FindSpace();
        virtual ~FindSpace();
        Point dir;                  //!< default = [1,1,1]
        bool allowContainerOverlap; //!< default = false;
        bool allowMatterOverlap;    //!< default = false;
        bool find(Geometrybase&, const p_vec&, p_vec&, unsigned int=1000);
    };

    /*!
     * \brief Vector rotation routines
     * \note Boundary condition are respected.
     * \author Mikael Lund
     * \date Canberra, 2009
     * \todo Redundant? See QuaternionRotateEigen.
     */
    class VectorRotate {
      private:
        virtual Point rotate(Point) const;      //!< Rotate point around axis
      protected:
        Point origin, u;
        double cosang, sinang;
        double e1mcox, e1mcoy, e1mcoz;
        Geometrybase *geoPtr;
      public:
        Point operator()(const Point&) const; // Rotate point around axis. 
        virtual ~VectorRotate();
        virtual void setAxis(Geometrybase&, const Point&, const Point&, double);  //!< Set rotation axis and degrees
        double getAngle() const;                                                  //!< Get set rotation angle
    };

    /*!
     * \brief Quaternion rotation routine using the Eigen library
     */
    class QuaternionRotateEigen : public VectorRotate {
      private:
        Eigen::Quaterniond q;
        inline Point rotate(Point a) const {
          a=a-origin;
          geoPtr->boundary(a);
          a=q*a+origin;
          geoPtr->boundary(a);
          return a;
        }
      public:
        inline void setAxis(Geometrybase &g, const Point &beg, const Point &end, double angle) {
          origin=beg;
          u=end-beg;
          g.boundary(u);
          u.normalize();
          geoPtr=&g;
          q=Eigen::AngleAxisd(angle, u);
          cosang=std::cos(angle);
        }
    };

    /*!
     * \brief Quaternion rotation routine
     * \todo Redundant?
     */
    class QuaternionRotate : public VectorRotate {
      private:
        double d1, d2, d3, d4, d5, d6, d7, d8, d9 ;
        Point rotate(Point) const;
        struct quat {             /* Define a quaternion structure */
          double w,x,y,z;
        };
        quat q;
      public:
        void setAxis(Geometrybase&, const Point&, const Point&, double);  //!< Set rotation axis and degrees
    };

    /*! \brief Calculate minimum distance between two line segments
     *
     * Find closest distance between line segments and return its vector
     * gets orientations and lengths of line segments and the vector connecting
     * their center os masses (from vec1 to vec2)
     * Copyright 2001, softSurfer (www.softsurfer.com)
     * This code may be freely used and modified for any purpose
     * providing that this copyright notice is included with it.
     * SoftSurfer makes no warranty for this code, and cannot be held
     * liable for any real or imagined damage resulting from its use.
     * Users of this code must verify correctness for their application.
     */
    Point mindist_segment2segment(const Point&, double, const Point&, double, const Point&);
    Point mindist_segment2point(const Point&, double, const Point&);

    inline Point vec_perpproject(const Point &A, const Point &B) {
      Point x;
      x=A - B* (A.dot(B));
      return x;
    }
    int test_intrpatch(const CigarParticle &, Point &, double , double , double [5]);
    int find_intersect_plane(const CigarParticle &, const CigarParticle &, const Point &, const Point &, double , double , double [5]);
    int find_intersect_planec(const CigarParticle &, const CigarParticle &, const Point &, const Point &, double , double , double [5]);
    int psc_intersect(const CigarParticle &, const CigarParticle &, const Point &, double [5], double );  
    int cpsc_intersect(const CigarParticle &, const CigarParticle &,const Point &, double [5], double );



  }//namespace Geometry
}//namespace Faunus
#endif
