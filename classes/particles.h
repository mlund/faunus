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
  enum keys {XYZ, ALL, TRIAL, ACCEPT}; 
  
  vector<particle> p;                    //!< The main particle vector
  vector<particle> trial;                //!< Trial particle vector. 
  macromolecule append(vector<particle>);//!< Add particle vector to both "p" and "trial"
  
  int push_back(particle &);                            //!< add particle to both "p" and "trial"
  void zmove(group &, double, keys=TRIAL);              //!< Translate group in z-direction only
  void zundo(group &);                                  //!< Undo space::zmove
  void move(group &, point, keys=TRIAL);                //!< Translation group in x,y,z direction
  void displace(int, double);                           //!< Random displacement of particle i
  void undo(group &, keys=XYZ);                         //!< Regret MC move. Restores space::trial.
  void accept(group &, keys=XYZ);                       //!< Acceps MC move. Sync space::trial and space::p

  void rotate(group &, double, keys=TRIAL);
  void rotate(group &, point, double=0, keys=TRIAL);    //!< Rotate group around point

  double charge();                      //!<  Sum all charges in particle vector
  double charge(point &, double);       //!<  Sum all charges within a sphere region
};

#endif
