#ifndef FAUNUS_POINT_H
#define FAUNUS_POINT_H

#ifndef SWIG
#include "faunus/common.h"
#endif

namespace Faunus {
  class random;
  class geometrybase;
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

      template<typename Tgeometry> inline double sqdist(const point &b) {
        return Tgeometry::sqdist(*this,b);
      }

  };

  inline int point::anint(double a) const { return int(a>0 ? a+.5 : a-.5); }

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
      pointparticle& operator=(const specdata&);
  };

  class cigarparticle : public pointparticle {
    public:
      point omega, patch;
      double patchangle, length;
      void rotate(const geometrybase&, const point&, double);            //!< Rotate around a vector
      void translate(const geometrybase&, const point&);                 //!< Translate along a vector
      void scale(const geometrybase&, double);                           //!< Volume scaling
      bool overlap(const cigarparticle&, double) const;     //!< Check for overlap of two particles
      friend std::ostream &operator<<(std::ostream &, const cigarparticle &); //!< Output information
      cigarparticle &operator<<(std::istream &);                              //!< Get information
      cigarparticle operator+(const point&) const;         //!< Add two vectors
      cigarparticle& operator=(const specdata&);
      cigarparticle& operator=(const pointparticle&);
  };

}//namespace
#endif
