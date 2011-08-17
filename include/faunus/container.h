#ifndef FAU_CONTAINER_H
#define FAU_CONTAINER_H

#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/slump.h"

namespace Faunus {
  class inputfile;
  class group;
  
  class space {
  protected:
    std::ifstream fin;
    std::ofstream fout;
  private:
    slump slp;
  public:
    vector<particle> p;                             //!< The main particle vector
    vector<particle> trial;                         //!< Trial particle vector. 
    vector<group*> g;                               //!< Pointers to all groups in the system.
    
    virtual bool saveToDisk(string);                    //!< Save container state to disk
    virtual bool loadFromDisk(string,bool=false);       //!< Load container state from disk
    
    void add(particle);
    bool insert(particle, unsigned int);            //!< Insert particle at pos n.
    bool remove(unsigned int);                      //!< Remove particle n.
    
    double charge() const;                          //!< Sum all charges in particle vector
    bool check_vector();                            //!< Check if p and trial are equal!
  };
  
  /*!
   * \brief Polymorphic class for simulation containers
   * \author Mikael Lund
   */
  class container : public space {
  protected:
    slump slp;
    double volume;                                      //!< Volume of the container [AA^3]
  public:
    double getvolume() {return volume;}
    virtual void setvolume(double);                     //!< Specify new volume
    virtual bool collision(const particle&)=0;          //!< Check for collision with walls
    virtual bool collision_internal(const particle&);   //!< Check for internal collisions - used to define "forbidden" zones
    virtual void randompos(point &)=0;                  //!< Random point within container
    virtual void boundary(point &) const=0;             //!< Apply boundary conditions to a point
    virtual void scale(point&, const double&) const;    //!< Scale point to a new volume - for Npt ensemble

    virtual double sqdist(const point&, const point&)=0;//!< Squared distance between two points
    virtual double dist(const point&,const point&);     //!< Distance between two points
    virtual point vdist(const point&, const point&)=0;
    virtual string info();                              //!< Return info string
 
    bool saveToDisk(string);                    //!< Save container state to disk
    bool loadFromDisk(string,bool=false);       //!< Load container state from disk
    
    
    /*!
     * \brief Test if given pair potential is compatible with the container (i.e. same boundary conditions)
     * \returns True if potential is OK - false otherwise.
     *
     * Measuring the distance between two randomly placed space 1000 times using the 
     * pair potential and the container distance methods it is asserted whether they agree.
     * The test also fails if the randompos() function generates a point that collides with
     * the container boundaries.
     */
    template<typename Tpairpot> bool check_potential(Tpairpot &pair) {
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
   *  \todo Implement space scaling for isobaric ensemble
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
   *  of some of the space to a part of the cuboid. The function slicecollision
   *  can be used to make sure space are positioned within in the slice.
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
    
    bool collision_internal(const particle &a);  //!< Check collision with slice edges
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
      point r=a-b;
      if (r.x>len_half.x)
        r.x-=len.x;
      else if (r.x<-len_half.x)
        r.x+=len.x;
      if (r.y>len_half.y)
        r.y-=len.y;
      else if (r.y<-len_half.y)
        r.y+=len.y;
      if (r.z>len_half.z)
        r.z-=len.z;
      else if (r.z<-len_half.z)
        r.z+=len.z;
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
  
  /*!
   * \brief Cuboidslit: cubuid without periodic boundary in the z direction
   * \author Chris Evers
   * \date Lund, nov 2010
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
      point r=a-b;
      if (r.x>len_half.x)
        r.x-=len.x;
      else if (r.x<-len_half.x)
        r.x+=len.x;
      if (r.y>len_half.y)
        r.y-=len.y;
      else if (r.y<-len_half.y)
        r.y+=len.y;
      return r;
    }
    
    //! Apply periodic boundary conditions
    inline void boundary(point &a) const {
      if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
      if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
    }
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
