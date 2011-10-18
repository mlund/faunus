#ifndef FAUNUS_POINT_H
#define FAUNUS_POINT_H

#ifndef SWIG
#include "faunus/common.h"
#endif

namespace Faunus {
  class RandomBase;
  class geometrybase;
  class AtomData;
  /*!
   * \brief Cartesian coordinates
   * \author Mikael Lund
   * \date 2002-2007
   */
  class Point {
    private:
      inline int anint(double) const;
    public:
      double x,y,z;                       //!< Cartesian coordinates
      Point();                            //!< Constructor, zero data.
      Point(double,double,double);        //!< Constructor, set vector
      void clear();                       //!< Zero all data.
      double len() const; 
      void ranunit(RandomBase &);               //!< Generate a random unit vector
      double dot(const Point &) const;      //!< Angle with another point
      Point operator-() const;              //!< Sign reversal
      const Point operator*(double) const;  //!< Scale vector
      Point operator+(const Point&) const;  //!< Add two vectors
      Point operator-(const Point&) const;  //!< Subtract vector
      Point & operator+=(const Point&);     //!< Vector addition
      Point & operator*=(const double);     //!< Scaling vector
      bool operator==(const Point&) const;
      std::string str();
      friend std::ostream &operator<<(std::ostream &, const Point &); //!< Output information
      Point &operator<<(std::istream &);                        //!< Get information

  };

  inline int Point::anint(double a) const { return int(a>0 ? a+.5 : a-.5); }

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
  class PointParticle : public Point {
    public:
      PointParticle();
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
      bool overlap(const PointParticle&, double) const;                       //!< Check for overlap of two particles
      friend std::ostream &operator<<(std::ostream &, const PointParticle &); //!< Output information
      PointParticle& operator<<(std::istream &);                              //!< Get information
      PointParticle& operator=(const Point&);                                 //!< Copy coordinates from a point
      PointParticle& operator=(const AtomData&);
  };

  class CigarParticle : public PointParticle {
    public:
      Point omega, patch;
      double patchangle, length;
      void rotate(const geometrybase&, const Point&, double);            //!< Rotate around a vector
      void translate(const geometrybase&, const Point&);                 //!< Translate along a vector
      void scale(const geometrybase&, double);                           //!< Volume scaling
      bool overlap(const CigarParticle&, double) const;     //!< Check for overlap of two particles
      friend std::ostream &operator<<(std::ostream &, const CigarParticle &); //!< Output information
      CigarParticle &operator<<(std::istream &);                              //!< Get information
      CigarParticle operator+(const Point&) const;         //!< Add two vectors
      CigarParticle& operator=(const Point&);
      CigarParticle& operator=(const AtomData&);
      CigarParticle& operator=(const PointParticle&);
  };

}//namespace
#endif
