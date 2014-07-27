#ifndef FAU_GEOMETRY_H
#define FAU_GEOMETRY_H

#ifndef SWIG
#include <faunus/point.h>
#include <faunus/slump.h>
#include <faunus/textio.h>
#include <faunus/physconst.h>
#include <Eigen/Geometry>
#endif

namespace Faunus {

  /**
   * @brief Namespace for geometric operations.
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
        double dist(const Point&,const Point&) const;       //!< Distance between two points (A)
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
        void setRadius(double);                 //!< Set radius (angstrom)
        Sphere(double);                         //!< Construct from radius (angstrom)
        Sphere(InputMap&, string="sphere");     //!< Construct from InputMap key \c prefix_radius
        void randompos(Point &);
        void boundary(Point &p) const {};
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
        inline bool collision(const particle &a, collisiontype type=BOUNDARY) const {
          if (std::abs(a.x())>len_half.x()) return true;
          if (std::abs(a.y())>len_half.y()) return true;
          if (std::abs(a.z())>len_half.z()) return true;
          return false;
        }

        /**
         * For reviews of minimum image algorithms,
         * see doi:10/ck2nrd and doi:10/kvs
         */
        double sqdist(const Point &a, const Point &b) const {
          Point d = (a-b).cwiseAbs();
          for (int i=0; i<3; ++i)
            if (d[i]>len_half[i]) d[i]-=2*len_half[i];
          return d.squaredNorm();

          // Alternative algorithm:
          // Eigen::Vector3d d = a-b;
          // Eigen::Vector3i k = (2*d.cwiseProduct(len_inv)).cast<int>();
          // return (d-k.cast<double>().cwiseProduct(len)).squaredNorm();
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
     * \brief Cuboid with no periodic boundaries in z direction
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
        inline Point vdist(const Point &a, const Point &b) FOVERRIDE {
          return a-b;
        }
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

        inline Point vdist(const Point &a, const Point &b) FOVERRIDE {
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

    /**
     * @brief Calculate center of cluster of particles
     * @param geo Geometry
     * @param p Particle vector
     * @param g Range of particle index (`Group`, `vector<int>`, ...)
     * @param weight Weight function
     *
     * The weight function is typically a lambda function that takes a particle
     * as an argument and returns the weight, for example Mw, charge, or unity
     * for mass center, charge center, or geometric center. Functions for
     * these three examples are predefined.
     */
    template<class Tgeo, class Tpvec, class TGroup, class Tweightf>
      Point anyCenter(const Tgeo &geo, const Tpvec &p, const TGroup &g, Tweightf weight) {
        assert(!p.empty());
        assert(g.back() < (int)p.size());
        Point cm(0,0,0);
        if (!g.empty()) {
          Point o = p[ g.front()+(g.back()-g.front())*0.5 ];  // set origo to middle particle
          double sum=0;
          for (auto i : g) {
            Point t = p[i]-o;       // translate to origo
            geo.boundary(t);        // periodic boundary (if any)
            cm += t * weight(p[i]);
            sum += weight(p[i]);
          }
          if (fabs(sum)<1e-6) sum=1;
          cm=cm/sum + o;
          geo.boundary(cm);
        }
        return cm;
      }

    /** @brief Calculate mass center of cluster of particles */
    template<class Tgeo, class Tpvec, class TGroup>
      Point massCenter(const Tgeo &geo, const Tpvec &p, const TGroup &g) {
        return anyCenter(geo,p,g,
            [](const typename Tpvec::value_type &x) {return x.mw;} );
      }

    /** @brief Calculate mass center of a particle vector */
    template<typename Tgeometry, typename Tp_vec>
      Point massCenter(const Tgeometry &geo, const Tp_vec &p) {
        if (p.empty())
          return Point(0,0,0);
        return massCenter(geo,p,Group(0,p.size()-1));
      }

    /** @brief Calculate charge center of cluster of particles */
    template<class Tgeo, class Tpvec, class TGroup>
      Point chargeCenter(const Tgeo &geo, const Tpvec &p, const TGroup &g) {
        return anyCenter(geo,p,g,
            [](const typename Tpvec::value_type &x) { return std::fabs(x.charge); } );
      }

    /** @brief Calculate charge center of a particle vector */
    template<typename Tgeometry, typename Tp_vec>
      Point chargeCenter(const Tgeometry &geo, const Tp_vec &p) {
        if (p.empty())
          return Point(0,0,0);
        return chargeCenter(geo,p,Group(0,p.size()-1));
      }

    /** @brief Calculate geometric center of cluster of particles */
    template<class Tgeo, class Tpvec, class TGroup>
      Point geometricCenter(const Tgeo &geo, const Tpvec &p, const TGroup &g) {
        return anyCenter(geo,p,g,
            [](const typename Tpvec::value_type &x) {return 1.0;} );
      }

    template<class Tspace, class Tgroup>
      Point dipoleMoment(const Tspace &s, const Tgroup &g, double cutoff=1e9,Point mu=Point(0,0,0)) {
        assert(g.size()<=(int)s.p.size());
        for (auto i : g) {
          Point t=s.p[i] - g.cm;
          s.geo.boundary(t);
          if(t.squaredNorm() < cutoff*cutoff)
            mu += t*s.p[i].charge;
        }
        return mu;
      }

    /** @brief Translate a particle vector by a vector */
    template<class Tgeo, class Tpvec>
      void translate(const Tgeo &geo, Tpvec &p, const Point &d) {
        for (auto &pi : p) {
          pi += d;
          geo.boundary(pi);
        }
      }

    /** @brief Translate a particle vector so mass center is in (0,0,0) */
    template<class Tgeo, class Tpvec>
      void cm2origo(const Tgeo &geo, Tpvec &p) {
        translate(geo, p, -massCenter(geo, p) );
      }

    /** @brief Translate a particle vector so charge center is in (0,0,0) */
    template<class Tgeo, class Tpvec>
      void cc2origo(const Tgeo &geo, Tpvec &p) {
        translate(geo, p, -chargeCenter(geo, p) );
      }

    /**
     * @brief Find an empty space for a particle vector in a space of other particles
     * @author Mikael Lund
     * @date Asljunga 2011
     * @todo Implement random rotation in addition to current translation scheme.
     */
    class FindSpace {
      private:
        template<class Tgeometry, class Tpvec>
          bool matterOverlap(const Tgeometry &geo, const Tpvec &p1, const Tpvec &p2) const {
            if (allowMatterOverlap==false)
              for (auto &i : p1)
                for (auto &j : p2) {
                  double max=i.radius+j.radius;
                  if ( geo.sqdist(i,j)<max*max )
                    return true;
                }
            return false;
          }

        template<class Tgeometry, class Tpvec>
          bool containerOverlap(const Tgeometry &geo, const Tpvec &p) const {
            if (allowContainerOverlap==false)
              for (auto &i : p)
                if (geo.collision(i)) return true;
            return false;
          }

      public:
        Point dir;                  //!< default = [1,1,1]
        bool allowContainerOverlap; //!< default = false;
        bool allowMatterOverlap;    //!< default = false;

        inline FindSpace() {
          dir=Point(1,1,1);
          allowContainerOverlap=false;
          allowMatterOverlap=false;   
        }

        /**
         * @param geo Geometry to use
         * @param dst Destination particle vector (will not be touched!)
         * @param p Particle vector to find space for. Coordinates will be changed.
         * @param maxtrials Number of times to try before timeout.
         */
        template<class Tgeometry, class Tpvec>
          bool find(Tgeometry &geo, const Tpvec &dst, Tpvec &p, int maxtrials=1e3) const {
            using namespace textio;
            cout << "Trying to insert " << p.size() << " particle(s)";
            Point v;
            do {
              cout << ".";
              maxtrials--;
              Point cm = massCenter(geo, p);
              geo.randompos(v);
              v = v.cwiseProduct(dir);
              translate(geo, p, -cm+v);
            } while (maxtrials>0 && (containerOverlap(geo,p) || matterOverlap(geo,p,dst)));
            if (maxtrials>0) {
              cout << " OK!\n";
              return true;
            }
            cout << " timeout!\n";
            assert(!"Timeout - found no space for particle(s).");
            return false;
          }
    };

    /**
     * @brief Quaternion rotation routine using the Eigen library
     * @note Boundary condition are respected.
     */
    class QuaternionRotate {
      private:
        double angle_;
        Eigen::Vector3d origin;
        Eigen::Quaterniond q;
        Eigen::Matrix3d rot_mat; // rotation matrix
        Geometrybase *geoPtr;

      public:
        //!< Get rotation origin
        Eigen::Vector3d& getOrigin() { return origin; }

        //!< Get rotation origin
        Eigen::Vector3d getOrigin() const { return origin; }

        //!< Get set rotation angle
        double getAngle() const { return angle_; }

        bool ignoreBoundaries;

        QuaternionRotate() {
          ignoreBoundaries=false;
        }

        /**
         * @brief Set rotation axis and angle
         * @param g Geometry to use for periodic boundaries (if any)
         * @param beg Starting point for vector to rotate around
         * @param end Ending point for vector to rotate around
         * @param angle Radians to rotate
         */
        inline void setAxis(Geometrybase &g, const Point &beg, const Point &end, double angle) {
          assert(&g!=nullptr && "Invalid geometry");
          geoPtr=&g;
          origin=beg;
          angle_=angle;
          Point u(end-beg); //Point u(end-beg);
          assert(u.squaredNorm()>0 && "Rotation vector has zero length");
          g.boundary(u);
          u.normalize(); // make unit vector
          q=Eigen::AngleAxisd(angle, u);

          rot_mat << 0, -u.z(), u.y(),u.z(),0,-u.x(),-u.y(),u.x(),0;
          rot_mat = Eigen::Matrix3d::Identity() + rot_mat*std::sin(angle) + rot_mat*rot_mat*(1-std::cos(angle));
        }

        /** @brief Rotate point - respect boundaries */
        inline Point operator()(Point a) const {
          if(ignoreBoundaries)
            return q*a;
          a=a-origin;
          geoPtr->boundary(a);
          a=q*a+origin;
          geoPtr->boundary(a);
          return a;
        }

        inline Eigen::Matrix3d operator()(Eigen::Matrix3d a) const {
          a = rot_mat*a*rot_mat.transpose();
          return a;
        }
    };

    /**
     * @brief Calculates the volume of a collection of particles
     *
     * This will use a brute force, stochastic hit and miss algorithm to
     * calculate the net volume of a collection of overlapping
     * particles.
     *
     * @param p Particle vector (structure must be whole)
     * @param n Number of iterations (default: 1e7)
     * @param pradius Probe radius (default: 0)
     *
     * @todo More efficient sampling may be achieved by adjusting the
     *       box vs. molecular volume (i.e. 1:1)
     */
    template<typename Tpvec>
      double calcVolume(const Tpvec &p, unsigned int n=1e7, double pradius=0) {
        double L=0;      // size of test box
        Point gc(0,0,0); // geometric center of molecule
        for (auto &i : p)
          gc += i / p.size();
        for (auto &i : p)
          L = std::max(L, 2*((i-gc).norm() + i.radius + pradius));

        // Start shooting!
        Point r;
        unsigned int hit=0, cnt=0;
        while (++cnt<n) {
          r.x() = slp_global.randHalf();
          r.y() = slp_global.randHalf();
          r.z() = slp_global.randHalf();
          r = r*L + gc;
          for (auto &i : p)
            if ((i-r).squaredNorm()<pow(i.radius+pradius,2)) {
              hit++;
              break;
            }
        }
        return hit/double(cnt) * pow(L,3);
      }

  }//namespace Geometry
}//namespace Faunus
#endif
