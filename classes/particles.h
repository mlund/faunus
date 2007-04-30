#ifndef _particles_h
#define _particles_h

#include <vector>
#include <string>
#include "slump.h"
#include "group.h"

/*!
 * \brief Class the contains all particles and manipulating methods.
 * \author Mikael Lund
 * \date Lund, 2004
 *
 * Handles all particles in a system and provides manipulating
 * functions like translation, rotation, etc. It contains the
 * main particle vector space::p
 */
class particles {
private:
  slump slp;
public:
  enum keys {XYZ, ALL, TRIALMOVE, AUTOACCEPT,MAX,AVERAGE,INTERIOR}; 
  
  vector<particle> p;                   //!< The main particle vector
  vector<particle> trial;               //!< Trial particle vector. 
  group append(vector<particle>);       //!< Add particle vector to both "p" and "trial"
  
  int push_back(particle &);                            //!< add particle to both "p" and "trial"
  void zmove(group &, double, keys=TRIALMOVE);          //!< Translate group in z-direction only
  void zundo(group &);                                  //!< Undo space::zmove
  void move(group &, point, keys=TRIALMOVE);            //!< Translation group in x,y,z direction
  void displace(int, double);                           //!< Random displacement of particle i
  void undo(group &, keys=XYZ);                         //!< Regret MC move. Restores space::trial.
  void accept(group &, keys=XYZ);                       //!< Acceps MC move. Sync space::trial and space::p

  void rotate(group &, double, keys=TRIALMOVE);
  void rotate(group &, point, double=0, keys=TRIALMOVE);//!< Rotate group around point

  double radius(group &, point &, keys=MAX);  //!<  Calculate radius of group, centered in point
  double charge(group::group &);        //!<  Update netcharge of group
  double charge();                      //!<  Sum all charges in particle vector
  double charge(point &, double);       //!<  Sum all charges within a sphere region
  
  double recalc_dipole(group &);          //!< Calc. unit vector dipole and scalar
  point mass_center(group &);             //!< Calc. center-of-mass
};

#endif
