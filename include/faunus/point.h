#ifndef FAUNUS_POINT_H
#define FAUNUS_POINT_H

#ifndef SWIG
#include <faunus/common.h>
#include <Eigen/Core>
#endif

namespace Faunus {

#ifndef FAU_HYPERSPHERE
  typedef Eigen::Vector3d BasePoint;
#else
  typedef Eigen::Vector4d BasePoint;
#endif

  /**
   * @brief Cartesian coordinates
   *
   * This is the base class for all particles and takes care
   * of positions, only. It is derived from Eigen::Vector3d
   * vector and all particles in faunus can hence be freely
   * mixed with Eigen objects.
   *
   * @date 2002-2007
   */
  struct Point : public BasePoint {
    typedef double Tcoord;                             //!< Floating point type for Point coordinates
    typedef BasePoint Tvec;                            //!< 3D vector from Eigen
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

    /** @brief Zero data */
    void clear() {
      setZero();
    }

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

    Tcoord len() const {
      return norm();
    }

    /** @brief Read from stream */
    Point& operator<<(std::istream &in) {
      for (int i=0; i<size(); ++i)
        in >> (*this)[i];
      return *this;
    }

    /**
     * @brief Transform point (rotation etc)
     * @param rotator Functor that rotates a point and returns the rotated Point
     *
     * The functor should take care of simulation boundaries (if any) and typically one
     * would want to pass the Geometry::VectorRotate class as in the following example:
     * @code
     * Point a(1,0,0);
     * VectorRotate rotator;
     * rotator.setAxis(geometry, Point(0,0,0), Point(0,0,1), 3.14 ); // rotate pi around 0,0,1
     * a.rotate(rotator);
     * @endcode
     */
    void rotate(RotFunctor rotator) {
      *this = rotator(*this);
    }
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

  /**
   * @brief Hypersphere particle
   * @date Lund, 2009-2013
   * @warning Unfinished - need to transfer from jurassic branch
   */
  class HyperParticle : public PointParticle {
    private:

      inline const Point::Tcoord& z1() const { return x(); }
      inline const Point::Tcoord& z2() const { return y(); }
      inline const Point::Tcoord& z3() const { return z(); }
      inline const Point::Tcoord& z4() const { return w(); }
      inline Point::Tcoord& z1() { return x(); }
      inline Point::Tcoord& z2() { return y(); }
      inline Point::Tcoord& z3() { return z(); }
      inline Point::Tcoord& z4() { return w(); }

    public:

      inline void clear() {
        setZero();
        z4()=1;
      }

      inline HyperParticle() { clear(); }

      /** @brief Copy constructor for Eigen derivatives */
      template<typename OtherDerived>
        HyperParticle(const Eigen::MatrixBase<OtherDerived>& other) : PointParticle(other) {}

      /** @brief Generic copy operator for Eigen derivatives */
      template<typename OtherDerived>
        HyperParticle& operator=(const Eigen::MatrixBase<OtherDerived> &other) {
          PointParticle::operator=(other);
          return *this;
        }

      /** @brief Copy operator for base class (i.e no casting to Eigen derivatives) */
      inline HyperParticle& operator=(const PointParticle &p) {
        PointParticle::operator=(p);
        return *this;
      }

      /** @brief Copy properties from AtomData object */
      inline HyperParticle& operator=(const AtomData &d) {
        PointParticle::operator=(d);
        return *this;
      }

      /**
       * @brief Squared distance between two points.
       * @return \f[ r^2 = z_1z_1' + z_2z_2' + z_3z_3' + z_4z_4' \f]
       */
      inline double sqdist(const HyperParticle &p) const {
        return z1()*p.z1() + z2()*p.z2() + z3()*p.z3() + z4()*p.z4();
      }

      /**
       * @brief Geodesic distance between two hyperpoints
       * @return \f[ r_{\mbox{\scriptsize{geod}}} = \arccos{ (r^2) } \f]
       */
      inline double geodesic(const HyperParticle &p) const {
        return std::acos(sqdist(p));
      }

      template<typename Tgeometry>
        void translate(const Tgeometry &geo, const Point &a) {
          double du=a.x();
          double dv=a.y();
          double dw=a.z();
          double nz1, nz2, nz3, nz4,
                 tz1, tz2, tz3, tz4,
                 rho=du, omega=dv, fi=dw;
          nz1=std::sqrt(1.-rho*rho);
          nz2=nz1*std::cos(fi);
          nz1=nz1*std::sin(fi);
          nz3=rho*std::sin(omega);
          nz4=rho*std::cos(omega);

          HyperParticle e1,e2,e3,te1,te2,te3;
          double fact1,fact2,fact3,nabla_nb,fi_nb;

          //nabla_nb=slp.random_one()*2.*acos(-1.);
          //fi_nb=std::acos(slp.random_one());

          e1.z1()=std::cos(nabla_nb);
          e1.z2()=std::sin(nabla_nb);
          e1.z3()=0;
          e1.z4()=0;
          e2.z1()=-std::cos(fi_nb)*std::sin(nabla_nb);
          e2.z2()=std::cos(fi_nb)*std::cos(nabla_nb);
          e2.z3()=std::sin(fi_nb);
          e2.z4()=0;
          e3.z1()=std::sin(fi_nb)*std::sin(nabla_nb);
          e3.z2()=-std::sin(fi_nb)*std::cos(nabla_nb);
          e3.z3()=std::cos(fi_nb);
          e3.z4()=0;

          // First create a random orthonormal basis set at North Pole
          fact1=e1.z1()*z1() +e1.z2()*z2() +e1.z3()*z3();
          te1.z1()=e1.z1()-1./(1.+z4())*fact1*z1();
          te1.z2()=e1.z2()-1./(1.+z4())*fact1*z2();
          te1.z3()=e1.z3()-1./(1.+z4())*fact1*z3();
          te1.z4()=e1.z4()-1./(1.+z4())*fact1*(z4()+1.);

          fact2=e2.z1()*z1() +e2.z2()*z2() +e2.z3()*z3();
          te2.z1()=e2.z1()-1./(1.+z4())*fact2*z1();
          te2.z2()=e2.z2()-1./(1.+z4())*fact2*z2();
          te2.z3()=e2.z3()-1./(1.+z4())*fact2*z3();
          te2.z4()=e2.z4()-1./(1.+z4())*fact2*(z4()+1.);

          fact3=e3.z1()*z1()+e3.z2()*z2()+e3.z3()*z3();
          te3.z1()=e3.z1()-1./(1.+z4())*fact3*z1();
          te3.z2()=e3.z2()-1./(1.+z4())*fact3*z2();
          te3.z3()=e3.z3()-1./(1.+z4())*fact3*z3();
          te3.z4()=e3.z4()-1./(1.+z4())*fact3*(z4()+1.);

          // Then move it to point of z1,z2,z3,z4
          tz1=nz1*te1.z1()+nz2*te2.z1()+nz3*te3.z1()+nz4*z1();
          tz2=nz1*te1.z2()+nz2*te2.z2()+nz3*te3.z2()+nz4*z2();
          tz3=nz1*te1.z3()+nz2*te2.z3()+nz3*te3.z3()+nz4*z3();
          tz4=nz1*te1.z4()+nz2*te2.z4()+nz3*te3.z4()+nz4*z4();

          // Update the point
          z1()=tz1;
          z2()=tz2;
          z3()=tz3;
          z4()=tz4;
        }
  };

}//namespace
#endif
