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
    virtual string info();                //!< Print information
    void operator+=(group);
    group operator+(group);

    void undo(particles &);
    void accept(particles &);                           //!< Accept a move
    void add(particles &, vector<particle>);            //!< Add a particle vector
    void add(container &, particle::type, short);       //!< Add particles w. collision check
    virtual unsigned short displace(container&,double); //!< Displace random particle
};

class macromolecule : public group {
  public:
    macromolecule();
    point mu;            //!< Dipole moment
    average<float> Q;    //!< Total charge. Updated with space::charge()
    average<float> Q2;   //!< Total charge squared.
    average<float> dip;  //!< Dipole moment scalar.

    string info();                      //!< Show info
    double charge(vector<particle> &);  //!< Calculate total charge
    double radius(vector<particle> &);  //!< Calculate radius
    double dipole(vector<particle> &);  //!< Calculate dipole moment
    void move(particles &, point);      //!< Translate group
    void BOXmove(particles &, box &, point);
    void zmove(particles &, double);    //!< Move in z-direction, only
    void rotate(particles &, double);          //!< Rotate around a point
    void rotate(particles &, point, double);
    void operator=(group);              //!< Copy from group
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

