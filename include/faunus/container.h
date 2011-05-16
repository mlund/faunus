#ifndef FAU_CONTAINER_H
#define FAU_CONTAINER_H

#include "faunus/particles.h"
#include "faunus/point.h"
#include "faunus/slump.h"

namespace Faunus {
  class inputfile;
 
  /*! \brief Polymorphic class for simulation containers
   *  \author Mikael Lund
   */
  class container : public particles {
    protected:
      slump slp;
      double volume;                                      //!< Volume of the container [AA^3]
    public:
      double getvolume() {return volume;}
      virtual void setvolume(double);                     //!< Specify new volume
      virtual bool collision(const particle &)=0;         //!< Check for collision with walls
      virtual bool slicecollision(const particle&);       //!< Check for collision with the slice borders 
      virtual void randompos(point &)=0;                  //!< Random point within container
      virtual string info();                              //!< Return info string
      virtual string povray();                            //!< POVRAY object representing the cell
      virtual void boundary(point &) const=0;             //!< Apply boundary conditions to a point
      virtual void scale(point&, const double&) const;    //!< Scale point to a new volume length
      virtual bool saveToDisk(string);                    //!< Save container date to disk
      virtual bool loadFromDisk(string,bool=false);       //!< Load container data from disk
      virtual double sqdist(const point&, const point&)=0;//!< Squared distance between two points
      virtual double dist(const point&,const point&);     //!< Distance between two points
      virtual point vdist(const point&, const point&);

      /*!
       * \brief Test if given pair potential is compatible with the container (i.e. same boundary conditions)
       * \returns True if potential is OK - false otherwise.
       *
       * Measuring the distance between two randomly placed particles 1000 times using the 
       * pair potential and the container distance methods it is asserted whether they agree.
       * The test also fails if the randompos() function generates a point that collides with
       * the container boundaries.
       */
      template<typename Tpairpot> bool testpotential(Tpairpot &pair) {
        particle a,b;
        for (int i=0; i<1e3; ++i) {
          randompos(a);
          randompos(b);
          if (collision(a)==true) return false;
          if (collision(b)==true) return false;
          if (std::abs(pair.sqdist(a,b)-sqdist(a,b)) > 1e-6) return false;
        }
        return true;
      }
  };

  /*! \brief Spherical simulation container
   *  \author Mikael Lund
   *  \todo Implement particles scaling for isobaric ensemble
   *
   *  This is a spherical simulation container, surrounded by a hard wall.
   */
  class cell : public container {
    private:
      double r2,diameter;
    public:
      void setradius(double);
      double r;              //!< Radius
      cell(double);
      cell(inputfile &);
      string info();
      void setvolume(double);
      void randompos(point &);
      void boundary(point &) const;
      string povray();
      bool collision(const particle &);
      inline double sqdist(const point &p1, const point &p2) { return p1.sqdist(p2); }
  };

  //---------------------------------------------------------
  /*! \brief Cuboid simulation container with periodic boundaries
   *
   *  \author Chris Evers
   *  \date Lund, nov 2010
   *
   *  The cuboid simulation container has right angles, rectangular faces 
   *  and periodic boundaries. A slice can be introduced to constrain the position
   *  of some of the particles to a part of the cuboid. The function slicecollision
   *  can be used to make sure particles are positioned within in the slice.
   */
  class cuboid : public container {
    protected:
      bool setslice(point, point);             //!< Reset slice position
      point len_inv;                           //!< Inverse sidelengths

    public:
      cuboid(inputfile &);                     //!< Read input parameters
      bool setlen(point);                      //!< Reset cuboid sidelengths

      point len;                               //!< Sidelengths
      point len_half;                          //!< Half sidelength
      point slice_min, slice_max;              //!< Position of slice corners

      string info();                           //!< Return info string

      point randompos();                       //!< Get point with random position
      void randompos(point &);                 //!< Move point to random position

      bool slicecollision(const particle &a);  //!< Check collision with slice edges
      bool collision(const particle &a) {      //!< Check collision with cuboid edges
        if (std::abs(a.x) > len_half.x ||
            std::abs(a.y) > len_half.y ||
            std::abs(a.z) > len_half.z  )
          return true;
        return false;
      }

      //! Calculate distance using the minimum image convention
      inline double sqdist(const point &p1, const point &p2) {   //!< Squared distance 
        return p1.sqdist_mi_xyz(p2, len, len_half);
      }

      inline point vdist(const point &a, const point &b) {       //!< Distance vector
        point r;
        r.x=std::abs(a.x-b.x);
        r.y=std::abs(a.y-b.y);
        r.z=std::abs(a.z-b.z);
        if (r.x>len_half.x) r.x-=len.x;
        if (r.y>len_half.y) r.y-=len.y;
        if (r.z>len_half.z) r.z-=len.z;
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
  };

  //---------------------------------------------------------
  /*! \brief Cuboidslit: cubuid without periodic boundary in the z direction
   *  
   *  \author Chris Evers
   *  \date Lund, nov 2010
   */
  class cuboidslit : public cuboid {
    public:
      cuboidslit(inputfile &);
      string info();

      //! Calculate distance using the minimum image convention
      inline double sqdist(const point &p1, const point &p2) {   //!< Squared distance 
        return p1.sqdist_mi_xy(p2, len, len_half);
      }

      inline point vdist(const point &a, const point &b) {       //!< Distance vector
        point r;
        r.x=std::abs(a.x-b.x);
        r.y=std::abs(a.y-b.y);
        r.z=a.z-b.z;
        if (r.x>len_half.x) r.x-=len.x;
        if (r.y>len_half.y) r.y-=len.y;
        return r;
      }

      //! Apply periodic boundary conditions
      inline void boundary(point &a) const {
        if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
        if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
      }
  };

  /*!
   * \brief Cubic simulation container w. periodic boundaries
   * \author Mikael Lund
   *
   * This is a cubic simulation container with periodic boundaries.
   * It is more or less superseeded by the cuboid container so consider
   * using that instead.
   */
  class box : public container {
    private:
      bool setlen(double);                 //!< Reset box length
    public:
      double len_half, len_inv, tlen_inv;  //!< tlen is the trial box
      void setvolume(double);              //!< Set volume (and sidelength) of box
      double len;                          //!< Side length
      box(double);
      box(inputfile &);
      string info();
      string povray();
      void randompos(point &);
      point randompos();

      bool collision(const particle &a) {
        if (std::abs(a.x)>len_half) return true;
        if (std::abs(a.y)>len_half) return true;
        if (std::abs(a.z)>len_half) return true;
        return false;
      }

      bool clash(const particle &, const particle &);

      //! Calculate distance using the minimum image convention
      inline double sqdist(const point &a, const point &b) {
        return a.sqdist_mi_xyz(b, len, len_half);
      }

      inline point vdist(const point &a, const point &b) {
        point r;
        r.x=std::abs(a.x-b.x);
        r.y=std::abs(a.y-b.y);
        r.z=std::abs(a.z-b.z);
        if (r.x>len_half) r.x-=len;
        if (r.y>len_half) r.y-=len;
        if (r.z>len_half) r.z-=len;
        return r;
      }

      inline int anint(double x) const {
        return int(x>0. ? x+.5 : x-.5);
      }

      /*!
       * \brief Apply periodic boundary conditions
       */
      inline void boundary(point &a) const {
        if (std::abs(a.x)>len_half) a.x-=len*anint(a.x/len);
        if (std::abs(a.y)>len_half) a.y-=len*anint(a.y/len);
        if (std::abs(a.z)>len_half) a.z-=len*anint(a.z/len);
      }

      /*!
       * \brief Linear scaling of particle position after a volume fluctuation
       */
      inline void scale(point &a, const double &newlen) const {
        a = a*(newlen/len);
      }
  };

  //---------------------------------------------------------
  /*! \brief Box with periodic boundaries in the x and y direction.
   *  
   *  \author Mikael Lund
   *  \date Asljunga, 2008
   */
  class slit : public box {
    public:
      double xyarea, zlen, zlen_half; //crossecional area 
      slit(inputfile &);
      string info();
      void randompos(point &);

      inline void boundary(point &a) const {
        a.x=a.x-len*anint(a.x*len_inv);
        a.y=a.y-len*anint(a.y*len_inv);
      }

      inline double sqdist(const point &p1, const point &p2) {
        double dz=p1.z-p2.z,
               dx=std::abs(p1.x-p2.x),
               dy=std::abs(p1.y-p2.y);
        if (dx>len_half) dx-=len;
        if (dy>len_half) dy-=len;
        return dx*dx + dy*dy + dz*dz;
      }

      inline point vdist(const point &a, const point &b) {
        point r;
        r.x=std::abs(a.x-b.x);
        r.y=std::abs(a.y-b.y);
        r.z=a.z-b.z;
        if (r.x>len_half) r.x-=len;
        if (r.y>len_half) r.y-=len;
        return r;
      }

      bool collision(const particle &a) {
        if (std::abs(a.x)>len_half ||
            std::abs(a.y)>len_half ||
            std::abs(a.z)>zlen_half ) { return true; }
        return false;
      }
  };

  class xyplane : public box {
    public:
      xyplane(inputfile &in);
      void randompos(point &);
      void randompos(vector<point> &);
  };

  
  /*! \brief Gridded slit box.
   *  
   *  \author Mikael Lund
   *  \date Asljunga, 2010
   */
  class gridslit : public slit {
  private:
    double l;     //!< Grid spacing
  public:
    double ngrid; //!< Number of grids in each dimension
    gridslit(inputfile &);
    string info();
    bool collision(const particle &);
    void randompos(point &);
  };
  
  
  /*! \brief "Clutch" like container.
   *  \author Mikael Lund
   *
   *  A spherical cell with a particle inaccessible area shaped
   *  as a disc in the middle of the sphere. The disc is parallel
   *  to the XY-plane and spans two Z-values as specified in the
   *  constructor.
   *
   *  \image html clutch.png
   */
  class clutch : public container {
    private:
      double r2;
      double diameter;
    public:
      double r,zmin,zmax;
      clutch(double, double, double);
      void randompos(point &);
      bool collision(const particle &);
  };

  /*! \brief Cylindrical simulation container
   *  \author Mikael Lund/Bjoern Persson
   *  \todo Needs some testing
   *
   *  This is a cylinder container where all walls
   *  are HARD. The origin is in the middle of the
   *  cylinder.
   */
  class cylinder : public container {
    private:
      double halflen;
      double r2;    //!< Cylinder radius squared
    public:
      double len;   //!< Cylinder length
      double r;     //!< Cylinder radius
      double diameter;
      cylinder(double, double);
      cylinder(inputfile &);
      void randompos(point &);
      bool collision(const particle &);
      string info(); //!< Cylinder info
      string povray();
  };

#ifdef HYPERSPHERE
  /*! \brief Hypersphere simulation container
   *  \author Martin Trulsson
   *  \date Lund, 2009
   */
  class hypersphere : public cell {
    private:
      static const double pi;
    public:
      hypersphere(inputfile &);
      string info();
      void randompos(point &);
      bool collision(const particle &);

      double dist(const point &a, const point &b) {
        return r*a.geodesic(b); // CHECK!!! (virtual=slow!)
      }

      inline double sqdist(const point &a, const point &b) {
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
}//namespace
#endif
