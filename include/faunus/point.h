#ifndef FAUNUS_POINT_H
#define FAUNUS_POINT_H

#ifndef SWIG
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#include <Eigen/Geometry>

#pragma GCC diagnostic pop
#include <vector>

#endif

namespace Faunus
{

  class AtomData;

  /**
   * @brief Cartesian coordinates
   *
   * This is the base class for all particles and takes care
   * of positions, only. It is derived from `Eigen::Vector3d`
   * vector and all particles in faunus can hence be freely
   * mixed with Eigen objects.
   */
  struct PointBase : public Eigen::Vector3d
  {
      typedef double Tcoord;        //!< Floating point type for Point coordinates
      typedef Eigen::Vector3d Tvec; //!< 3D vector from Eigen

      /** @brief Default constructor. Data is *not* zeroed */
      inline PointBase() {}

      PointBase( Tcoord x, Tcoord y, Tcoord z ) : Tvec(x, y, z) {}

      template<typename OtherDerived>
      PointBase( const Eigen::MatrixBase<OtherDerived> &other ) : Tvec(other) {}

      template<typename OtherDerived>
      PointBase &operator=( const Eigen::MatrixBase<OtherDerived> &other )
      {
          Tvec::operator=(other);
          return *this;
      }

      /** @brief Zero data */
      void clear() { setZero(); }

      /**
       * @brief Convert into random unit vector
       *
       * Based on the von Neumann method described in
       * *Allen and Tildesley*, page 349, which ensures
       * a uniform distribution on a unit sphere. More information
       * about *Sphere Point Picking* can be found at
       * [MathWorld](http://mathworld.wolfram.com/SpherePointPicking.html).
       *
       * @param ran Function that takes no arguments and returns a random
       *            float uniformly distributed in the range `[0:1[`.
       */
      template<class Trandombase>
      PointBase &ranunit( Trandombase &ran )
      {
          Tcoord r2;
          do
          {
              for ( size_t i = 0; i < 3; ++i )
                  this->operator[](i) = 2 * ran() - 1;
              r2 = squaredNorm();
          }
          while ( r2 > 1 );
          *this = *this / std::sqrt(r2);
          assert(norm() - 1 < 1e-7); // is it really a unit vector?
          return *this;
      }

      /**
       * @brief Translate along a vector
       * @param geo Geometry to use for boundary conditions (see Faunus::Geometry)
       * @param a Vector to translate with
       */
      template<typename Tgeometry>
      void translate( const Tgeometry &geo, const PointBase &a )
      {
          (*this) += a;
          geo.boundary(*this);
      }

      /**
       * @brief Coordinate scaling used for NPT and NVT ensemble
       *
       * This will perform a volume scaling of the Point by
       * following the algorithm specified in the Geometry.
       */
      template<typename Tgeometry>
      void scale( const Tgeometry &geo, PointBase &s, double xyz = 1, double xy = 1 )
      {
          geo.scale(*this, s, xyz, xy);
      }

      Tcoord len() const
      {
          return norm();
      }

      /** @brief Read from stream */
      PointBase &operator<<( std::istream &in )
      {
          for ( int i = 0; i < size(); ++i )
              in >> (*this)[i];
          return *this;
      }

      /** @brief Read from string */
      PointBase &operator<<( const std::string &in )
      {
          std::istringstream i(in);
          return operator<<(i);
      }

      /*
       * @brief Internal rotation. No effect on isotropic particle.
       */
      template<typename Trotator>
      void rotate( const Trotator &rot )
      {
      }

      /** @brief Fill from std vector */
      template<class T>
      PointBase &operator<<( const std::vector<T> &v )
      {
          // see http://stackoverflow.com/questions/26094379/typecasting-eigenvectorxd-to-stdvector
          assert(v.size() == 3);
          for ( size_t i = 0; i < 3; ++i )
              (*this)[i] = v.at(i);
          return *this;
      }

  };

  template<class T=double>
  struct Tensor : public Eigen::Matrix<T, 3, 3>
  {
      typedef Eigen::Matrix<T, 3, 3> Tmat; //!< Matrix from Eigen

      /** @brief Default constructor. Data is zeroed */
      Tensor()
      {
          (*this).setZero();
      }

      Tensor( T xx, T xy, T xz, T yy, T yz, T zz )
      {
          (*this) << xx, xy, xz, xy, yy, yz, xz, yz, zz;
      }

      template<typename OtherDerived>
      Tensor( const Eigen::MatrixBase<OtherDerived> &other ) : Tmat(other) {}

      template<typename OtherDerived>
      Tensor &operator=( const Eigen::MatrixBase<OtherDerived> &other )
      {
          Tmat::operator=(other);
          return *this;
      }

      /** @brief Read from stream */
      Tensor &operator<<( std::istream &in )
      {
          (*this).setZero();
          for ( int i = 0; i < (*this).rows(); i++ )
          {
              for ( int j = i; j < (*this).cols(); j++ )
              {
                  in >> (*this)(i, j);
                  (*this)(j, i) = (*this)(i, j);
              }
          }
          return *this;
      }

      /** @brief Read from string */
      Tensor &operator<<( const std::string &in )
      {
          std::istringstream i(in);
          return operator<<(i);
      }

      /** @brief Write data members to stream */
      friend std::ostream &operator<<( std::ostream &o, const Tensor<T> &t )
      {
          for ( int i = 0; i != t.rows(); ++i )
              for ( int j = i; j != t.cols(); ++j )
                  o << t(i, j) << " ";
          return o;
      }

      template<typename Trotator>
      void rotate( const Trotator &rot ) { *this = rot(*this); }

      void eye() { *this = Tmat::Identity(3, 3); }
  };

  /**
   * @brief Hypersphere particle
   * @date Lund, 2009-2013
   * @warning Unfinished - need to transfer from jurassic branch
   */
  class HyperPoint : public PointBase
  {
  private:
      double w;

  public:
      inline HyperPoint() {}

      inline HyperPoint( double z1, double z2, double z3, double w ) : PointBase(z1, z2, z3)
      {
          z4() = w;
      }

      inline HyperPoint( double z1, double z2, double z3 ) : PointBase(z1, z2, z3)
      {
          z4() = 1;
      }

      template<typename OtherDerived>
      HyperPoint( const Eigen::MatrixBase<OtherDerived> &other ) : PointBase(other) {}

      template<typename OtherDerived>
      HyperPoint &operator=( const Eigen::MatrixBase<OtherDerived> &other )
      {
          PointBase::operator=(other);
          return *this;
      }

      inline const double &z1() const { return x(); }

      inline const double &z2() const { return y(); }

      inline const double &z3() const { return z(); }

      inline const double &z4() const { return w; }

      inline double &z1() { return x(); }

      inline double &z2() { return y(); }

      inline double &z3() { return z(); }

      inline double &z4() { return w; }

      /** @brief Read from stream */
      HyperPoint &operator<<( std::istream &in )
      {
          PointBase::operator<<(in);
          in >> z4();
          return *this;
      }

      /** @brief Write to stream */
      friend std::ostream &operator<<( std::ostream &o, const HyperPoint &p )
      {
          o << PointBase(p) << " " << p.z4();
          return o;
      }

      /**
       * @brief Translate along a vector
       * @param geo Geometry to use for boundary conditions (see Faunus::Geometry) 
       * @param a Vector to translate with
       */
      template<typename Tgeometry>
      void translate( const Tgeometry &geo, const HyperPoint &a )
      {
          (*this) += a;
          geo.boundary(*this);
      }

      /**
       * @brief Coordinate scaling used for NPT ensemble
       *
       * This will perform a volume scaling of the Point by
       * following the algorithm specified in the Geometry.
       */
      template<typename Tgeometry>
      void scale( const Tgeometry &geo, double newvol )
      {
          geo.scale(*this, newvol);
      }

      inline double sqdist( const HyperPoint &a ) const
      {
          return dot(a) + z4() * a.z4();
      }

      /**
       * @brief Geodesic distance between two hyperpoints
       * @return @f[ r_g = \arccos{ (r^2) } @f]
       */
      inline double geodesic( const HyperPoint &a ) const
      {
          return std::acos(sqdist(a));
      }

      void move( double du, double dv, double dw )
      {
          double nz1, nz2, nz3, nz4,
              tz1, tz2, tz3, tz4,
              rho = du, omega = dv, fi = dw;
          nz1 = std::sqrt(1. - rho * rho);
          nz2 = nz1 * std::cos(fi);
          nz1 = nz1 * std::sin(fi);
          nz3 = rho * std::sin(omega);
          nz4 = rho * std::cos(omega);

          HyperPoint e1, e2, e3, te1, te2, te3;
          double fact1, fact2, fact3, nabla_nb, fi_nb;

          nabla_nb = 0;//slp.random_one()*2.*acos(-1.);
          fi_nb = 0;//std::acos(slp.random_one());

          e1.z1() = std::cos(nabla_nb);
          e1.z2() = std::sin(nabla_nb);
          e1.z3() = 0;
          e1.z4() = 0;
          e2.z1() = -std::cos(fi_nb) * std::sin(nabla_nb);
          e2.z2() = std::cos(fi_nb) * std::cos(nabla_nb);
          e2.z3() = std::sin(fi_nb);
          e2.z4() = 0;
          e3.z1() = std::sin(fi_nb) * std::sin(nabla_nb);
          e3.z2() = -std::sin(fi_nb) * std::cos(nabla_nb);
          e3.z3() = std::cos(fi_nb);
          e3.z4() = 0;

          // First create a random orthonormal basis set at North Pole
          fact1 = e1.z1() * z1() + e1.z2() * z2() + e1.z3() * z3();
          te1.z1() = e1.z1() - 1. / (1. + z4()) * fact1 * z1();
          te1.z2() = e1.z2() - 1. / (1. + z4()) * fact1 * z2();
          te1.z3() = e1.z3() - 1. / (1. + z4()) * fact1 * z3();
          te1.z4() = e1.z4() - 1. / (1. + z4()) * fact1 * (z4() + 1.);

          fact2 = e2.z1() * z1() + e2.z2() * z2() + e2.z3() * z3();
          te2.z1() = e2.z1() - 1. / (1. + z4()) * fact2 * z1();
          te2.z2() = e2.z2() - 1. / (1. + z4()) * fact2 * z2();
          te2.z3() = e2.z3() - 1. / (1. + z4()) * fact2 * z3();
          te2.z4() = e2.z4() - 1. / (1. + z4()) * fact2 * (z4() + 1.);

          fact3 = e3.z1() * z1() + e3.z2() * z2() + e3.z3() * z3();
          te3.z1() = e3.z1() - 1. / (1. + z4()) * fact3 * z1();
          te3.z2() = e3.z2() - 1. / (1. + z4()) * fact3 * z2();
          te3.z3() = e3.z3() - 1. / (1. + z4()) * fact3 * z3();
          te3.z4() = e3.z4() - 1. / (1. + z4()) * fact3 * (z4() + 1.);

          // Then move it to point of z1,z2,z3,z4
          tz1 = nz1 * te1.z1() + nz2 * te2.z1() + nz3 * te3.z1() + nz4 * z1();
          tz2 = nz1 * te1.z2() + nz2 * te2.z2() + nz3 * te3.z2() + nz4 * z2();
          tz3 = nz1 * te1.z3() + nz2 * te2.z3() + nz3 * te3.z3() + nz4 * z3();
          tz4 = nz1 * te1.z4() + nz2 * te2.z4() + nz3 * te3.z4() + nz4 * z4();

          // Update the point
          z1() = tz1;
          z2() = tz2;
          z3() = tz3;
          z4() = tz4;
      }
  };

#ifdef FAU_HYPERSPHERE
  typedef HyperPoint Point;
#else
  typedef PointBase Point;  //!< 3D vector
#endif

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
  struct PointParticle : public Point
  {
      typedef Point::Tcoord Tradius;
      typedef Point::Tcoord Tcharge;
      typedef Point::Tcoord Tmw;
      typedef Point::Tcoord Talphax;
      typedef unsigned char Tid;
      typedef bool Thydrophobic;
      Tid id;                                   //!< Particle identifier
      Tcharge charge;                           //!< Charge number
      Tradius radius;                           //!< Radius
      Talphax alphax;
      Tmw mw;                                   //!< Molecular weight
      Thydrophobic hydrophobic;                 //!< Hydrophobic flag

      PointParticle() { clear(); }              //!< Constructor

      template<typename OtherDerived>
      PointParticle( const Eigen::MatrixBase<OtherDerived> &other ) : Point(other) {}

      template<typename OtherDerived>
      PointParticle &operator=( const Eigen::MatrixBase<OtherDerived> &other )
      {
          Point::operator=(other);
          return *this;
      }

      Tcharge &q() { return charge; }

      Tcharge q() const { return charge; }

      template<class T,
          class = typename std::enable_if<std::is_base_of<AtomData, T>::value>::type>
      PointParticle &operator=( const T &d )
      {
          id = d.id;
          charge = d.charge;
          radius = d.radius;
          alphax = d.alphax;
          mw = d.mw;
          hydrophobic = d.hydrophobic;
          return *this;
      }

      /**
       * @brief Copy from stream
       *
       * This will read all data from stream in the same order as written.
       * Note that a short integer is expected for the particle id
       * since chars (Tid=char) does not print well on screen.
       * Derived classes should expand on this so
       * that *all* data is read.
       */
      PointParticle &operator<<( std::istream &in )
      {
          short tmp; // avoid char output in readable text files
          Point::operator<<(in);
          in >> charge >> radius >> mw >> tmp >> hydrophobic;
          id = (Tid) tmp;
          return *this;
      }

      /**
       * @brief Write to stream
       *
       * This will write all data to given stream. Note that the particle id is converted
       * to a short integer since char output (Tid=char) does not print well on screen.
       * Derived classes should expand on this so that *all* data is written.
       */
      friend std::ostream &operator<<( std::ostream &o, const PointParticle &p )
      {
          o << Point(p).transpose()
            << " " << p.charge << " " << p.radius << " " << p.mw << " "
            << (short) p.id << " " << p.hydrophobic;
          return o;
      }

      /** @brief Volume of particle */
      double volume() const
      {
          return 4 * std::acos(-1) * radius * radius * radius / 3;
      }

      /** @brief Zero data */
      void clear()
      {
          Point::clear();
          charge = mw = radius = alphax = 0;
          hydrophobic = false;
          id = 0;
      }

  };

  /**
   * @brief Dipolar particle
   */
  struct DipoleParticle : public PointParticle
  {
      Point mu;               //!< Dipole moment unit vector (permanent+induced)
      double muscalar;        //!< Dipole moment scalar (permanent+induced)
      Point mup;              //!< Permanent dipole moment vector
      Tensor<double> alpha;   //!< Polarization matrix
      Tensor<double> theta;   //!< Quadrupole matrix

      inline DipoleParticle() : mu(0, 0, 0), muscalar(0), mup(0, 0, 0) {};

      /** @brief Copy constructor for Eigen derivatives */
      template<typename OtherDerived>
      DipoleParticle( const Eigen::MatrixBase<OtherDerived> &other ) : PointParticle(other) {}

      /** @brief Generic copy operator for Eigen derivatives */
      template<typename OtherDerived>
      DipoleParticle &operator=( const Eigen::MatrixBase<OtherDerived> &other )
      {
          PointParticle::operator=(other);
          return *this;
      }

      /** @brief Copy operator for base class (i.e no casting to Eigen derivatives) */
      inline DipoleParticle &operator=( const PointParticle &p )
      {
          PointParticle::operator=(p);
          return *this;
      }

      /** @brief Copy properties from AtomData object */
      template<class T,
          class = typename std::enable_if<std::is_base_of<AtomData, T>::value>::type>
      DipoleParticle &operator=( const T &d )
      {
          PointParticle::operator=(d);
          muscalar = d.muscalar;
          mu = d.mu;
          mup = mu * muscalar;
          alpha = d.alpha;
          theta = d.theta;
          return *this;
      }

      /* read in same order as written! */
      inline DipoleParticle &operator<<( std::istream &in )
      {
          PointParticle::operator<<(in);
          mu.operator<<(in);
          in >> muscalar;
          mup.operator<<(in);
          alpha.operator<<(in);
          theta.operator<<(in);
          return *this;
      }

      /* write data members to stream */
      friend std::ostream &operator<<( std::ostream &o, const DipoleParticle &p )
      {
          o << PointParticle(p) << " " << p.mu.transpose() << " " << p.muscalar
            << " " << p.mup << " " << p.alpha << " " << p.theta;
          return o;
      }

      /**
       * @brief Internal rotation: dipole and polarizability
       */
      template<typename Trotator>
      void rotate( const Trotator &rot )
      {
          assert(rot.getOrigin().squaredNorm() < 1e-6);
          mu = rot(mu);
          mup = rot(mup);
          alpha = rot(alpha);
          theta = rot(theta);
      }
  };

  /**
   * @brief Sphero-cylindrical particle
   * @date Brno, November 2012
   *
   * detailed information here...
   */
  class CigarParticle : public PointParticle
  {
  public:
      Point dir; //!< Direction of sphero cylinder (unit vector)
      Point patchdir;
      Point patchsides[2];
      Point chdir;
      double patchangle;
      double pcanglsw;
      double pcangl;
      double halfl;
      double chiral_angle, panglsw, pangl;  // are these needed? From AtomData

      inline CigarParticle() : halfl(0) {}

      /** @brief Copy constructor for Eigen derivatives */
      template<typename OtherDerived>
      CigarParticle( const Eigen::MatrixBase<OtherDerived> &other ) : PointParticle(other) {}

      /**
       * @brief Initialize patchy spherocylinder - run at start and after patch changes
       *
       * Calculates cosine of angles, patch direction including chirality
       * and vector corresponding to sides of patch that are used in
       * calculations of interactions. 
       * This function must be called at the beginning of calculations and after changes
       * of patch properties.
       * It shall be also after a lot of move to remove accumulated comouptation errors
       */
      inline void init()
      {
          const double zero = 1e-9;
          if ( halfl > zero )
          {
              Point vec;
              Eigen::Quaterniond Q;
              pcangl = std::cos(0.5 * patchangle);
              pcanglsw = std::cos(0.5 * patchangle + panglsw);

              if ( dir.squaredNorm() < zero )
                  dir = {1, 0, 0};

              if ( patchdir.squaredNorm() < zero )
                  patchdir = {0, 1, 0};

              dir.normalize();

              patchdir = patchdir - dir * patchdir.dot(dir); // perp. project
              patchdir.normalize();

              /* calculate patch sides */
              if ( chiral_angle < zero )
                  vec = dir;
              else
              {
                  chdir = dir;
                  Q = Eigen::AngleAxisd(0.5 * chiral_angle, patchdir);
                  chdir = Q * chdir; // rotate
                  vec = chdir;
              }

              /* create side vector by rotating patch vector by half size of patch*/
              /* the first side */
              patchsides[0] = patchdir;
              Q = Eigen::AngleAxisd(0.5 * pangl + panglsw, vec);
              patchsides[0] = Q * patchsides[0]; // rotate
              patchsides[0].normalize();

              /* the second side */
              patchsides[1] = patchdir;
              Q = Eigen::AngleAxisd(-0.5 * pangl - panglsw, vec);
              patchsides[1] = Q * patchsides[1]; // rotate
              patchsides[1].normalize();

              if ( patchsides[0].squaredNorm() < zero )
                  throw std::runtime_error("Patch side vector has zero size.");
          }
      }

      /** @brief Generic copy operator for Eigen derivatives */
      template<typename OtherDerived>
      CigarParticle &operator=( const Eigen::MatrixBase<OtherDerived> &other )
      {
          PointParticle::operator=(other);
          return *this;
      }

      /** @brief Copy operator for base class (i.e no casting to Eigen derivatives) */
      inline CigarParticle &operator=( const PointParticle &p )
      {
          PointParticle::operator=(p);
          return *this;
      }

      /** @brief Copy properties from AtomData object */
      template<class T,
          class = typename std::enable_if<std::is_base_of<AtomData, T>::value>::type>
      CigarParticle &operator=( const T &d )
      {
          PointParticle::operator=(d);
          halfl = d.half_len;
          patchangle = d.pangl;
          pcanglsw = std::cos(0.5 * d.pangl + d.panglsw);
          pcangl = std::cos(0.5 * d.pangl);
          chiral_angle = d.chiral_angle;
          pangl = d.pangl;
          panglsw = d.panglsw;
          init();
          return *this;
      }

      template<typename Trotator>
      void rotate( const Trotator &rot )
      {
          if ( halfl > 1e-6 )
          {
              assert(rot.getOrigin().squaredNorm() < 1e-6);
              dir = rot(dir);
              patchdir = rot(patchdir);
              patchsides[0] = rot(patchsides[0]);
              patchsides[1] = rot(patchsides[1]);
              chdir = rot(chdir);
          }
      }

      /* read in same order as written! */
      inline CigarParticle &operator<<( std::istream &in )
      {
          PointParticle::operator<<(in);
          dir.operator<<(in);
          patchdir.operator<<(in);
          in >> halfl;
          return *this;
      }

      /* write data members to stream */
      friend std::ostream &operator<<( std::ostream &o, const CigarParticle &p )
      {
          o << PointParticle(p)
            << " " << p.dir.transpose() << " " << p.patchdir.transpose()
            << " " << p.halfl;
          return o;
      }
  };

}//namespace
#endif
