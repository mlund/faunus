#ifndef FAUNUS_POINT_H
#define FAUNUS_POINT_H

#ifndef SWIG
#include "faunus/common.h"
#endif

namespace Faunus {
  class random;
  class container;
  class specdata;
  /*!
   * \brief Cartesian coordinates
   * \author Mikael Lund
   * \date 2002-2007
   */
  class point {
    private:
      inline int anint(double) const;
    public:
      double x,y,z;                       //!< Cartesian coordinates
      point();                            //!< Constructor, zero data.
      point(double,double,double);        //!< Constructor, set vector
      void clear();                       //!< Zero all data.
      double len() const; 
      inline double sqdist(const point &) const;                                        //!< Squared distance to another point
      inline double sqdist_mi_xyz(const point &, const double &, const double &) const; //!< XYZ cubic minimum image distance  
      inline double sqdist_mi_xyz(const point&, const point&, const point&) const;      //!< XYZ minimum image distance
      inline double sqdist_mi_xy(const point&, const point&, const point&) const;
      inline double sqdist_mi_z(const point&, const double&, const double&) const;
      void ranunit(random &);               //!< Generate a random unit vector
      double dot(const point &) const;      //!< Angle with another point
      point operator-() const;              //!< Sign reversal
      point operator*(double) const;        //!< Scale vector
      point operator+(const point&) const;  //!< Add two vectors
      point operator-(const point&) const;  //!< Subtract vector
      point & operator+=(const point&);     //!< Vector addition
      point & operator*=(const double);     //!< Scaling vector
      bool operator==(const point&) const;
      std::string str();
      friend std::ostream &operator<<(std::ostream &, const point &); //!< Output information
      point &operator<<(std::istream &);                        //!< Get information
  };

  /*!
   * \return \f$ r_{12}^2 = \Delta x^2 + \Delta y^2 + \Delta z^2 \f$
   */
  inline double point::sqdist(const point &p) const {
    register double dx,dy,dz;
    dx=x-p.x;
    dy=y-p.y;
    dz=z-p.z;
    return dx*dx + dy*dy + dz*dz;
  }

  inline int point::anint(double a) const { return int(a>0 ? a+.5 : a-.5); }

  inline double point::sqdist_mi_xyz(const point& p, const double& len, const double& len_half) const {
    double dx=std::abs(x-p.x);
    double dy=std::abs(y-p.y);
    double dz=std::abs(z-p.z);
    if (dx>len_half) dx-=len;
    if (dy>len_half) dy-=len;
    if (dz>len_half) dz-=len;
    return dx*dx + dy*dy + dz*dz;
  }

  inline double point::sqdist_mi_xyz(const point &p, const point& len, const point& len_half) const {   //!< Squared distance 
    double dx=std::abs(x-p.x);
    double dy=std::abs(y-p.y);
    double dz=std::abs(z-p.z);
    if (dx>len_half.x) dx-=len.x;
    if (dy>len_half.y) dy-=len.y;
    if (dz>len_half.z) dz-=len.z;
    return dx*dx + dy*dy + dz*dz;
  }                                   

  inline double point::sqdist_mi_xy(const point &p, const point& len, const point& len_half) const {   //!< Squared distance 
    double dx=std::abs(x-p.x);
    double dy=std::abs(y-p.y);
    double dz=z-p.z;
    if (dx>len_half.x) dx-=len.x;
    if (dy>len_half.y) dy-=len.y;                                      
    return dx*dx + dy*dy + dz*dz;
  }   

  inline double point::sqdist_mi_z(const point &p, const double& len_z, const double& len_half_z) const {   //!< Squared distance 
    double dx=x-p.x;
    double dy=y-p.y;
    double dz=std::abs(z-p.z);
    if (dz>len_half_z) dx-=len_z;
    return dx*dx + dy*dy + dz*dz;
  }   

  /*!
   * \brief Class for particles
   * \author Mikael Lund
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
  class pointparticle : public point {
    public:
      pointparticle();
      double charge;                            //!< Charge number
      double radius;                            //!< Radius
      float mw;                                 //!< Molecular weight
      unsigned char id;                         //!< Particle identifier
      bool hydrophobic;                         //!< Hydrophobic flag
      double volume() const;                    //!< Return volume of sphere
      double mw2vol(double=1) const;            //!< Estimate volume from weight
      double mw2rad(double=1) const;            //!< Estimate radius from weight
      void deactivate();                        //!< Deactivate for use w. faster energy loops
      void clear();                             //!< Clear/reset all data
      bool overlap(const pointparticle&, double) const;                       //!< Check for overlap of two particles
      friend std::ostream &operator<<(std::ostream &, const pointparticle &); //!< Output information
      pointparticle& operator<<(std::istream &);                              //!< Get information
      pointparticle& operator=(const point&);                                 //!< Copy coordinates from a point
      pointparticle operator=(const specdata&) const;
   };

  class cigarparticle : public pointparticle {
    public:
      point omega, patch;
      double patchangle, length;
      void rotate(container&, const point&, double);            //!< Rotate around a vector
      void translate(container&, const point&);                 //!< Translate along a vector
      void scale(container&, double);                           //!< Volume scaling
      bool overlap(const cigarparticle&, double) const;     //!< Check for overlap of two particles
      friend std::ostream &operator<<(std::ostream &, const cigarparticle &); //!< Output information
      cigarparticle &operator<<(std::istream &);                              //!< Get information
      cigarparticle operator+(const point&) const;         //!< Add two vectors
      cigarparticle operator=(const specdata&) const;
      cigarparticle& operator=(const pointparticle&);
  };

#ifdef CIGARPARTICLE
  typedef cigarparticle particle;
#else
  typedef pointparticle particle;
#endif

}//namespace
#endif
