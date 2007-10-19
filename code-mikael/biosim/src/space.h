/*
 * SPACE (M. Lund, 2004)
 *
 * Handles all particles in a system and provides
 * manipulating functions like translation, rotation
 * etc. It contains the particle vector and allows
 * to save/load configurations to disk
 *
 */

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
#include "interact.h"
#include "cell.h"
#include "geometry.h"
#include "average.h"

using namespace std;


//! Class the contains all particles and manipulating methods.
/*!
 * Handles all particles in a system and provides manipulating
 * functions like translation, rotation, etc. It contains the
 * main particle vector space::p
 */
class space : private slump, private geometry {
public:
  enum keys {XYZ, ALL, TRIALMOVE, AUTOACCEPT,MAX,AVERAGE,INTERIOR};         //keywords 
  
  vector<particle> p;           ///< The main particle vector
  vector<particle> trial;       ///< Trial particle vector. IMPORTANT: After a MC move this MUST be in sync with space::p
  group append(vector<particle>);///< Add a particle vector to both space::p and space::trial
  
  bool save(string);            ///< Save particle vector to disk (x,y,z,charge)
  bool load(string);            ///< Load particle vector from disk.

  group insert_chain(chain::chain &);
  group insert_salt(int, double, double, cell::cell &, int=0);
  group insert_species(int, particle, cell::cell &);

  int push_back(particle &);    //add particle to both "p" and "trial"
  group sort(group &);           //put charged species first
  void zmove(group &, double, keys=TRIALMOVE);   ///< Translate group in z-direction only
  void zundo(group &);                           ///< Undo space::zmove
  void move(group &, point, keys=TRIALMOVE);    //general translation
  void displace(int, double);           ///< Random displacement of particle i
  void undo(group &, keys=XYZ);         ///< Regret MC move. Restores space::trial.
  void accept(group &, keys=XYZ);       ///< Acceps MC move. Sync space::trial and space::p
  void rotate(group &, double, keys=TRIALMOVE);
  void rotate(group &, point, double=0, keys=TRIALMOVE);  ///< Rotate group around point

  bool safetytest_vector();         ///< Test if space::p and space::trial are identical.

  point center_of_mass(group &);    // calc. center-of-mass
  void addmasscenter(group &);      // calc. and add mass center to group
  double radius(group &, point &, keys=MAX);  ///<  Calculate radius of group, centered in point
  double charge(group::group &);    ///<  Update netcharge of group
  double charge();                  ///<  Sum all charges in particle vector
  double charge(point &, double);   ///<  Sum charges that lies within a sphere
  double rho(point &, group &, double, double, double); ///< Number density
  int count(group &, double);    ///< Count charges
  int cntHydrophobic(group::group &); ///< Count number of hydrophobic residues
  
  point calcdipole(group &);
  bool adddipole(group &);

};

#endif
