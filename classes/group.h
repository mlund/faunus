#ifndef group_h
#define group_h

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "point.h"
#include "average.h"

using std::ostream;
using namespace std;

/*! \brief Groups set of particles, such as molecules, salt etc.
 *  \author Mikael Lund
 *  \date Lund 2004
 *
 *  A group defines the range of a a protein, the mobile ions,
 *  a chain or whatever. If is used extensively in MC movements,
 *  energy calculations and many, many other places.
 */
class group {
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
    void operator=(group);              //!< Copy from group
};

/*! \brief Freely jointed chain with harmonic spring potentials
 *  \author Mikael Lund
 */
class chain : public group {
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

#endif

