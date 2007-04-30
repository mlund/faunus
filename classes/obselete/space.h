#ifndef space_h
#define space_h

#include <valarray>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "point.h"
#include "group.h"
#include "slump.h"
#include "chain.h"
#include "hardsphere.h"
#include "potentials.h"
#include "container.h"
#include "average.h"

using namespace std;


/*!
 * \brief Class the contains all particles and manipulating methods.
 * \author Mikael Lund
 * \date Lund, 2004
 *
 * Handles all particles in a system and provides manipulating
 * functions like translation, rotation, etc. It contains the
 * main particle vector space::p
 */
class space : private slump {
public:
  enum keys {XYZ, ALL, TRIALMOVE, AUTOACCEPT,MAX,AVERAGE,INTERIOR}; 
  
  vector<particle> p;           //!< The main particle vector
  vector<particle> trial;       //!< Trial particle vector. IMPORTANT: After a MC move this MUST be in sync with space::p
  group append(vector<particle>);//!< Add a particle vector to both space::p and space::trial
  
  bool save(string);            //!< Save particle vector to disk (x,y,z,charge)
  bool load(string);            //!< Load particle vector from disk.

  group insert_chain(chain::chain &);
  group insert_salt(int, double, double,
      container::container &, particle::type=particle::GHOST);    //!< Insert mobile ions
  group insert_salt(int, particle, container::container &);       //!< Insert mobile ions

  int push_back(particle &);                            //!< add particle to both "p" and "trial"
  group sort(species &, group &);                       //!< put charged species first
  void zmove(group &, double, keys=TRIALMOVE);          //!< Translate group in z-direction only
  void zundo(group &);                                  //!< Undo space::zmove
  void move(group &, point, keys=TRIALMOVE);            //!< Translation group in x,y,z direction
  void displace(int, double);                           //!< Random displacement of particle i
  void undo(group &, keys=XYZ);                         //!< Regret MC move. Restores space::trial.
  void accept(group &, keys=XYZ);                       //!< Acceps MC move. Sync space::trial and space::p
  void rotate(group &, double, keys=TRIALMOVE);
  void rotate(group &, point, double=0, keys=TRIALMOVE);//!< Rotate group around point

  bool safetytest_vector();         //!< Test if space::p and space::trial are identical.

  double radius(group &, point &, keys=MAX);  //!<  Calculate radius of group, centered in point
  double charge(group::group &);        //!<  Update netcharge of group
  double charge();                      //!<  Sum all charges in particle vector
  double charge(point &, double);       //!<  Sum all charges within a sphere region
  double rho(point &, group &, double, double, double); ///< Number density
  int count(group &, double);           //!< Count charges
  
  double recalc_dipole(group &);          //!< Calc. unit vector dipole and scalar
  point mass_center(group &);             //!< Calc. center-of-mass
};

#endif
