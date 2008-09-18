#ifndef FAU_POINT_H
#define FAU_POINT_H

#include "faunus/common.h"
#include "faunus/slump.h"

namespace Faunus {
  /*!
   * \brief Cartesian coordinates
   * \author mikaek lund
   * \date 2002-2007
   */
  class point {
    private:
      inline int anint(double);
    public:
      double x,y,z;                       ///< The coordinates
      point();                            ///< Constructor, zero data.
      point(double,double,double);        ///< Constructor, set vector
      void clear();                       ///< Zero all data.
      double len(); 
      inline double sqdist(const point &) const;      //!< Squared distance to another point
      inline double sqdist(const point &,             //!< -- / / -- 3D minimum image
          const double &, const double &) const;
      inline double dist(const point &) const;        ///< Distance to another point
      inline double dist(const point &, double &, double &) const ; //!< Distance to another point
      void ranunit(slump &);              ///< Generate a random unit vector
      double dot(point &);                ///< Angle with another point
      point operator-();                  ///< Sign reversal
      point operator*(point);             ///< Multiply two vectors
      point operator*(double);            ///< Scale vector
      point operator+(point);             ///< Add two vectors
      point operator-(point);             ///< Substract vector
      point operator+(double);            ///< Displace x,y,z by value
      void operator+=(point);
      friend std::ostream &operator<<(std::ostream &, point &); /// Print x,y,z
      std::string str();
  };

  /*!
   * \brief Class for particles
   * \author mikaek lund
   * \date 2002-2007
   *
   * Example\n
   * \code
   * vector<particle> p(2);
   * p[0].radius = 2.0;
   * p[1].z = 10;
   * p[0].overlap( p[1] ); --> false
   * \endcode
   */
  class particle : public point {
    public:
      particle();
      //! Particle type identifier
      enum type {FIRST=0,GLY,ALA,VAL,LEU,ILE,PHE,TRP,TYR,HIS,SER,THR,MET,CYS,
        ASP,GLN,GLU,ASN,LYS,ARG,PRO,UNK,NTR,CTR,NA,K,F,CL,BR,I,SO4,ION,CATION,ANION,GHOST,
        RNH3,RNH4,RCOOH,RCOO,HYDROPHOBIC,LAST}; 

      double charge;                      //!< Charge number
      double radius;                      //!< Radius
      float mw;                           //!< Molecular weight
      type id;                            //!< Particle identifier
      bool hydrophobic;                   //!< Hydrophobic flag
      inline bool overlap(particle &);    //!< Hardsphere overlap test
      inline bool overlap(particle &, double &);
      inline double potential(const point &);   //!< Electric potential in point
      double volume();                          //!< Return volume of sphere
      double mw2vol(double=1);                  //!< Estimate volume from weight
      double mw2rad(double=1);                  //!< Estimate radius from weight
      void operator=(point);                    //!< Copy coordinates from a point
  };

  /*! \brief Class for spherical coordinates
   *  \author mikaek lund
   *  \date Canberra, 2005-2006
   */
  class spherical : private slump {
    public:
      double r,     //!< Radial
             theta, //!< Zenith angle \f$[0:\pi]\f$
             phi;   //!< Azimuthal angle \f$[0:2\pi]\f$
      spherical(double=0,double=0,double=0);
      point cartesian();                            //!< Convert to cartesian coordinates
      void operator=(point &);                      //!< Convert from cartesian coordinates
      void random_angles();                         //!< Randomize angles
  };

  inline void spherical::operator=(point &p) {
    r=p.len();
    theta=acos(p.z/r);
    phi=asin( p.y/sqrt(p.x*p.x+p.y*p.y) );
    if (p.x<0)
      phi=acos(-1.) - phi;
  }

  inline point spherical::cartesian() {
    point p;
    p.x=r*sin(theta)*cos(phi);
    p.y=r*sin(theta)*sin(phi);
    p.z=r*cos(theta);
    return p;
  }

  //! \todo This function is not completed
  inline void spherical::random_angles() {
    r=1.0 ;
  }

  /*!
   * \return \f$ |r_{12}|^2 = \Delta x^2 + \Delta y^2 + \Delta z^2 \f$
   */
  inline double point::sqdist(const point &p) const {
    register double dx,dy,dz;
    dx=x-p.x;
    dy=y-p.y;
    dz=z-p.z;
    return dx*dx + dy*dy + dz*dz;
  }

  //!\note <cmath> has a round() function -- speed?
  inline int point::anint(double a) { return int(a>0 ? a+.5 : a-.5); }
  inline double point::sqdist(const point &p, const double &len, const double &inv_len) const {
    register double dx,dy,dz;
    /*  dx=x-p.x;
        dy=y-p.y;
        dz=z-p.z;
        dx=dx-len*anint(dx*inv_len);
        dy=dy-len*anint(dy*inv_len);
        dz=dz-len*anint(dz*inv_len);
        return dx*dx + dy*dy + dz*dz;*/
    dx=abs(x-p.x);
    dy=abs(y-p.y);
    dz=abs(z-p.z);
    if (dx>len*0.5) dx-=len;
    if (dy>len*0.5) dy-=len;
    if (dz>len*0.5) dz-=len;
    return dx*dx + dy*dy + dz*dz;
  }
  /*inline double point::sqdist(point &p, double &len, double &inv_len) {
    register double dx,dy,dz;
    dx=x-p.x;
    dy=y-p.y;
    dz=z-p.z;
    return sqdf_(&dx, &dy, &dz, &len, &inv_len);
    }*/
  inline double point::dist(const point &p) const { return sqrt(sqdist(p)); }
  inline double point::dist(const point &p, double &len, double &inv_len) const { 
    return sqrt(sqdist(p, len, inv_len)); }

  /*!
   * \return \f$ \phi = \frac{z}{r_{12}}\f$
   * \note Not multiplied with the Bjerrum length!
   */
  inline double particle::potential(const point &p) { return charge / dist(p); }

  /*!
   * \return True if \f$ r_{12}<(\sigma_1+\sigma_2)/2 \f$ - otherwise false.
   */
  inline bool particle::overlap(particle &p) {
    double r=radius+p.radius;
    return (sqdist(p) < r*r) ? true : false;
  }
  inline bool particle::overlap(particle &p, double &s) {
    double r=radius+p.radius+s;
    return (sqdist(p) < r*r) ? true : false;
  }

}//namespace
#endif
