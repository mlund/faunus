#ifndef group_h
#define group_h

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "point.h"
#include "average.h"
#include "particles.h"
#include "container.h"
#include "pot_spring.h"
#include "slump.h"
#include "inputfile.h"

/*! \brief Groups set of particles, such as molecules, salt etc.
 *  \author Mikael Lund
 *  \date Lund 2004
 *
 *  A group defines the range of a a protein, the mobile ions,
 *  a chain or whatever. If is used extensively in MC movements,
 *  energy calculations and many, many other places.
 */
class group {
  protected:
    slump slp;
    string title;
  public:
    particle cm, cm_trial; 
    short int beg,end;    ///< Define range in particle vector. [beg;end]
    string name;          ///< Informative (arbitrary) name
    group(int=0);         ///< Constructor, initialize data.

    void set(short int,short int);        ///< Set particle range, "beg" and "end".
    short int size();                     ///< Number of particles in group
    short int random();                   ///< Picks a random particle within this group
    bool find(unsigned int);              ///< Check if particle is part of the group
    virtual double charge(vector<particle> &);//!< Calculate total charge
    point masscenter(vector<particle> &); //!< Calculate center-of-mass
    point masscenter(container &);        //!< Calc. center-of-mass
    virtual string info();                //!< Print information
    void operator+=(group);
    group operator+(group);

    void move(container &, point);                      //!< Translate group
    bool overlap(container &);                          //!< Test overlap w all particles
    void undo(particles &);
    void accept(particles &);                           //!< Accept a move
    void add(container &, vector<particle>, bool=false);//!< Add a particle vector
    void add(container &, particle::type, short);       //!< Add particles w. collision check
    short int count(vector<particle> &, particle::type);//!< Count number of specific particles
    virtual unsigned short displace(container&,double); //!< Displace random particle
    virtual void isobaricmove(container &);             //!< Scale coordinates for a volume fluctuation
};

/*!
 * \brief Group for handling salt particles
 */
class salt : public group {
  private:
    short nanion, ncation;
  public:
    salt(particle::type=particle::NA,
        particle::type=particle::CL);
    particle::type anion, cation; //!< Anion and cation types
    double muex;                  //!< Excess chemical potential
    void add(container &, inputfile &);                 //!< Add salt as specified in config file
    string info(container &);     //!< Show info
};

class macromolecule : public group {
  public:
    macromolecule();
    point mu;            //!< Dipole moment
    average<float> Q;    //!< Total charge. Updated with space::charge()
    average<float> Q2;   //!< Total charge squared.
    average<float> dip;  //!< Dipole moment scalar.

    string info();                      //!< Show info
    void center(container &);           //!< Center group in origo (0,0,0)
    double charge(vector<particle> &);  //!< Calculate total charge
    double radius(vector<particle> &);  //!< Calculate radius
    double dipole(vector<particle> &);  //!< Calculate dipole moment
    void zmove(container &, double);    //!< Move in z-direction, only
    void rotate(container &, double);   //!< Rotate around a point
    void rotate(container &, point, double);
    using group::add;
    void add(container &, inputfile &); //!< Add according to inputfile
    void operator=(group);              //!< Copy from group
    void isobaricmove(container &);     //!< Displace CM with scale difference
};

/*! \brief Freely jointed chain with harmonic spring potentials
 *  \author Mikael Lund
 */
class chain : public group, private pot_spring {
  public:
    chain();
    double k;                   //!< Spring constant
    double req;                 //!< Equilibrium distance between monomers
    short int graftpoint;       //!< Chain is grafted to this point. -1 if free (default)

    double monomerenergy(vector<particle> &, short);    //!< Spring energy of a monomer
    double internalenergy(vector<particle> &);          //!< Internal spring energy
    //!< Spring potential
    inline double quadratic(particle &p1, particle &p2) {
      double r=p1.dist(p2)-req;
      return k*r*r;
    }
};

class planarsurface : public group {
  public:
    point plane;                //!< Plane definition
    point offset;               //!< Plane translation
    planarsurface();
    void project(point &);      //!< Project point onto plane
};

/*! \brief This class is used for zwitter ions present on a surface.
 *  \author Mikael Lund
 *
 *  One end of the zwitter ion (even particles) can move only on the
 *  surface while the other ends (odd particles) are allowed to move
 *  in any direction. The two ends are connected by a spring.
 */
class zwittermembrane : public planarsurface, private pot_spring {
  public:
    short mate(short);                                  //!< Find Zwitter ion partners (two "mates")
    void add(container&, particle, particle, short=1);  //!< Add Zwitter ion(s)
    unsigned short displace(container&, double);        //!< Displace random particle
    double selfenergy(particles &);                     //!< Return internal spring energy
};

#endif

