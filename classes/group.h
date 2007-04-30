#ifndef group_h
#define group_h

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
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
  friend ostream &operator<<(ostream &, group &);
  point cm, cm_trial, mu; 

  short int beg,end;    ///< Define range in particle vector. [beg;end]
  string name;          ///< Informative (arbitrary) name
  group(int=0);         ///< Constructor, initialize data.
  
  void set(short int,short int);///< Set particle range, "beg" and "end".
  short int size();             ///< Number of particles in group
  short int random();           ///< Picks a random particle within this group
  bool find(unsigned int);      ///< Check if particle is part of the group
  void operator+=(group);
  group operator+(group);
};

class macromolecule : public group {
  public:
    double radius;        ///< Smallest sphere that can include the group. Set with with space::radius()
    average<double> Q;    ///< Total charge. Updated with space::charge()
    average<double> Q2;   ///< Total charge squared.
    average<double> dip;  ///< Dipole moment scalar.
};

/*! \brief Freely jointed chain with harmonic spring potentials
 *  \author Mikael Lund
 */
class chain : public group {
  public:
    double k;                   //!< Spring constant
    double req;                 //!< Equilibrium distance between monomers
    short int graftpoint;       //!< Chain is grafted to this point. -1 if free (default)
    double springpot(particle &, particle &); //!< Spring energy between two monomers.
};

#endif

