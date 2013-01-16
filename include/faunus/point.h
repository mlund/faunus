#ifndef FAUNUS_POINT_H
#define FAUNUS_POINT_H

#ifndef SWIG
#include <faunus/common.h>
#include <Eigen/Core>
#endif

namespace Faunus {
  /**
   * @brief Cartesian coordinates
   *
   * This is the base class for all particles and takes care
   * of positions, only. It is derived from Eigen::Vector3D
   * vector and all particles in faunus can hence be freely
   * mised with Eigen objects.
   *
   * @date 2002-2007
   */
  struct Point : public Eigen::Vector3d {
    typedef double Tcoord;                             //!< Floating point type for Point coordinates
    typedef Eigen::Vector3d Tvec;                      //!< 3D vector from Eigen
    typedef std::function<Tvec(const Tvec)> RotFunctor;//!< Rotation functor

    /** @brief Default constructor. Data is *not* zeroed */
    inline Point() {}

    Point(Tcoord x, Tcoord y, Tcoord z) : Tvec(x,y,z) {}

    template<typename OtherDerived>
      Point(const Eigen::MatrixBase<OtherDerived>& other) : Tvec(other) {}

    template<typename OtherDerived>
      Point& operator=(const Eigen::MatrixBase<OtherDerived> &other) {
        Tvec::operator=(other);
        return *this;
      }

    void clear();                           //!< Zero x,y,z

    Tcoord len() const;                     //!< Get scalar (TO BE REMOVED)

    /**
     * @brief Generate a random unit vector
     *
     * Based on the von Neumann method described in
     * *Allen and Tildesley*, page 349, which ensures
     * a uniform distribution on a unit sphere. More information
     * about *Sphere Point Picking* can be found at
     * [MathWorld](http://mathworld.wolfram.com/SpherePointPicking.html).
     *
     * @param ran Randomnumber function object. Must return a uniform
     *            distribution in the range `[0:1[`. This will typically be
     *            a class derived from `RandomBase`.
     */
    template<typename Trandombase>
      void ranunit(Trandombase &ran) {
        Point::Tvec u;
        Tcoord r2;
        do {
          u.x()=2*( ran()-0.5 );
          u.y()=2*( ran()-0.5 );
          u.z()=2*( ran()-0.5 );
          r2=u.squaredNorm();
        } while (r2>1);
        *this = u/std::sqrt(r2);
        assert(std::abs(norm()-1)<1e-7); // is it really a unit vector?
      }

    void rotate(RotFunctor);                  //!< Transform point (rotation etc)

    /**
     * @brief Translate along a vector
     * @param geo Geometry to use for boundary conditions (see Faunus::Geometry) 
     * @param a Vector to translate with
     */
    template<typename Tgeometry>
      void translate(const Tgeometry &geo, const Point &a) {
        assert(&geo!=nullptr);
#ifndef __clang__
        static_assert(
            std::is_base_of<Geometry::Geometrybase, Tgeometry>::value,
            "Tgeo must be derived from Geometrybase" );
#endif
        (*this)+=a;
        geo.boundary(*this);
      }

    /**
     * @brief Coordinate scaling used for NPT ensemble
     *
     * This will perform a volume scaling of the Point by
     * following the algorithm specified in the Geometry.
     */
    template<typename Tgeometry>
      void scale(const Tgeometry &geo, double newvol) {
        geo.scale(*this, newvol);
      }

    Point& operator<<(std::istream&);            //!< Read from stream
  };

  /**
   * @brief Class for isotropic particles
   *
   * Example:
   *
   * ~~~
   * std::vector<PointParticle> p(2);
   * p[0].radius = 2.0;
   * p[1].z() = 10;
   * std::cout << p[0];
   * ~~~
   *
   */
  struct PointParticle : public Point {
    typedef Point::Tcoord Tradius;
    typedef Point::Tcoord Tcharge;
    typedef Point::Tcoord Tmw;
    typedef unsigned char Tid;
    typedef bool Thydrophobic;

    Tcharge charge;                           //!< Charge number
    Tradius radius;                           //!< Radius
    Tmw mw;                                   //!< Molecular weight
    Tid id;                                   //!< Particle identifier
    Thydrophobic hydrophobic;                 //!< Hydrophobic flag

    PointParticle();                          //!< Constructor

    template<typename OtherDerived>
      PointParticle(const Eigen::MatrixBase<OtherDerived>& other) : Point(other) {}

    template<typename OtherDerived>
      PointParticle& operator=(const Eigen::MatrixBase<OtherDerived> &other) {
        Point::operator=(other);
        return *this;
      }

    PointParticle& operator=(const AtomData&);//!< Copy data from AtomData
    PointParticle& operator<<(std::istream&); //!< Copy data from stream
    friend std::ostream &operator<<(std::ostream&, const PointParticle&);//!< Write to stream

    double volume() const;                    //!< Return volume
    void deactivate();                        //!< Deactivate for use w. faster energy loops
    void clear();                             //!< Zero all data
  };

  /**
   * @brief Dipolar particle
   */
  struct DipoleParticle : public PointParticle {
    Point mu;               //!< Dipole moment unit vector
    double muscalar;        //!< Dipole moment scalar
    Eigen::Matrix3d alpha;

    inline DipoleParticle() : mu(0,0,0), muscalar(0) {};

    /** @brief Copy constructor for Eigen derivatives */
    template<typename OtherDerived>
      DipoleParticle(const Eigen::MatrixBase<OtherDerived>& other) : PointParticle(other) {}

    /** @brief Generic copy operator for Eigen derivatives */
    template<typename OtherDerived>
      DipoleParticle& operator=(const Eigen::MatrixBase<OtherDerived> &other) {
        PointParticle::operator=(other);
        return *this;
      }

    /** @brief Copy operator for base class (i.e no casting to Eigen derivatives) */
    inline DipoleParticle& operator=(const PointParticle &p) {
      PointParticle::operator=(p);
      return *this;
    }

    /** @brief Copy properties from AtomData object */
    inline DipoleParticle& operator=(const AtomData &d) {
      PointParticle::operator=(d);
      // copy more atom properties here...
      return *this;
    }

    /* read in same order as written! */
    inline DipoleParticle& operator<<(std::istream &in) {
      PointParticle::operator<<(in);
      mu.operator<<(in);
      in >> muscalar;
      return *this;
    }

    /* write data members to stream */
    friend std::ostream &operator<<(std::ostream &o, const DipoleParticle &p) {
      o << PointParticle(p) << " " << p.mu << " " << p.muscalar;
      return o;
    }

  };

  /**
   * @brief Sphero-cylindrical particle
   * @date Brno, November 2012
   *
   * detailed information here...
   */
  class CigarParticle : public PointParticle {
    public:
      Point dir; //!< Direction of sphero cylinder (unit vector)
      Point patchdir;
      Point patchsides[2];
      Point chdir;
      double patchangle;
      double pcanglsw;
      double pcangl;
      double halfl;

      inline CigarParticle() : halfl(0) {}

      /** @brief Copy constructor for Eigen derivatives */
      template<typename OtherDerived>
        CigarParticle(const Eigen::MatrixBase<OtherDerived>& other) : PointParticle(other) {}

      /** @brief Generic copy operator for Eigen derivatives */
      template<typename OtherDerived>
        CigarParticle& operator=(const Eigen::MatrixBase<OtherDerived> &other) {
          PointParticle::operator=(other);
          return *this;
        }

      /** @brief Copy operator for base class (i.e no casting to Eigen derivatives) */
      inline CigarParticle& operator=(const PointParticle &p) {
        PointParticle::operator=(p);
        return *this;
      }

      /** @brief Copy properties from AtomData object */
      inline CigarParticle& operator=(const AtomData &d) {
        PointParticle::operator=(d);
        // copy more atom properties here...
        return *this;
      }

      inline void rotate(RotFunctor rot) {
        if (halfl>1e-6) {
          dir = rot(dir);
          patchdir = rot(patchdir);
          patchsides[0] = rot(patchsides[0]);
          patchsides[1] = rot(patchsides[1]);
          chdir = rot(chdir);
        } else
          Point::rotate(rot);
      }

      /* read in same order as written! */
      inline CigarParticle& operator<<(std::istream &in) {
        PointParticle::operator<<(in);
        dir.operator<<(in);
        patchdir.operator<<(in);
        in >> halfl;
        return *this;
      }

      /* write data members to stream */
      friend std::ostream &operator<<(std::ostream &o, const CigarParticle &p) {
        o << PointParticle(p)
          << " " << p.dir << " " << p.patchdir
          << " " << p.halfl;
        return o;
      }
  };

#ifdef HYPERSPHERE
  /**
   * @brief Hypersphere particle
   * @author Martin Trulsson
   * @date Lund, 2009
   * @warning Unfinished - need to transfer from jurassic branch
   */
  class Hyperpoint : public PointParticle {
    public:
      double z1,z2,z3,z4;                     //!< Reduced Coordinates on hypersphere
      translate(const Geometry::Geometrybase&, const Point&);
      friend std::ostream &operator<<(std::ostream&, Hyperpoint&);
      Hyperpoint &operator<<(std::istream&);

      void clear() {
        z1=z2=z3=0;
        z4=1;
      }

      Hyperpoint() { clear(); }

      /*!
       * \brief Squared distance between two points.
       * \return \f[ r^2 = z_1z_1' + z_2z_2' + z_3z_3' + z_4z_4' \f]
       */
      inline double sqdist(const Hyperpoint &p) const {
        return z1*p.z1+z2*p.z2+z3*p.z3+z4*p.z4;
      }

      /*!
       * \brief Geodesic distance between two hyperpoints
       * \return \f[ r_{\mbox{\scriptsize{geod}}} = \arccos{ (r^2) } \f]
       */
      inline double geodesic(const Hyperpoint &p) const {
        return std::acos(sqdist(p));
      }
  };
#endif

}//namespace
#endif
