/*
 * Interact class.
 *
 * Calculates intermolecular interactions
 * between particles and groups.
 *  
 * M. Lund, 2004
 *     
 */

#include <vector>
#include <cstdarg>
#include "intrinsics.h"
#include "group.h"
#include "geometry.h"
#include "slump.h"

#ifndef interact_h
#define interact_h


//! Class to test for collision in a particle vector.
/*!
 The main idea is to do the MC move in the space::trial
 particle vector and then parse this to the overlap
 functions. All functions return TRUE if there is a 
 collision.
 If FALSE is returned, no collision is detected and you
 can proceed with energy calcultations.
*/
class collision {
 public:
  bool overlap(vector<particle> &, int);                              ///< all<->particle i.
  bool overlap(vector<particle> &, particle &);                       ///< all<->arbitrary (external) particle.
  bool overlap(vector<particle> &, group::group &);                   ///< all<->group.
  bool overlap(vector<particle> &, group::group &, int);              ///< group<->particle i.
  bool overlap(vector<particle> &, group::group &, particle &);       ///< group<>arbitrary (external) particle
  bool overlap(vector<particle> &, group::group &, group::group &);   ///< group<->group.
  bool overlap(vector<particle> &, vector<short int> &, double);      ///< internal collisions within subset
  bool celloverlap(vector<particle> &, group::group &, double);       ///< group with a spherical cell
  bool chgoverlap(vector<particle> &, group::group &, double);        ///< charge overlap within group
};


//! Class to calculate interactions
/*!
 * This class can calculate interactions within a particle
 * vector. You can do the calculation between groups, particles
 * etc.
 */
class interact {
private:
  slump::slump slmp;  ///< We need random numbers...
public:
  interact(double=7.12591);   ///< Constructor, sets the Bjerrum length (aangstom)
  
  double lB;   ///< Bjerrum length (can be specified in the constructor).
  double vdw;  ///< vdW interaction parameter.
  double k;    ///< quadratic force constant.
  double r0;   ///< quadratic eq. distance.
  double kappa;///< inverse screening length (in aangstroms)
  
  double energy(particle &, particle &);                      ///< particle<->particle (NOT IN KT!).
  double energy(vector<particle> &, int);                     ///< all<->particle i.
  double energy(vector<particle> &, group &);                 ///< all<->group.
  double energy(vector<particle> &, group &, group &);        ///< group<->group.
  double energy_vdw(vector<particle> &, group &, group &);    ///< group<->group (vdW only!).
  double energy(vector<particle> &, group &, int);            ///< group<->particle i.
  double energy(vector<particle> &, group &, particle &);     ///< group<->external particle.
  double energy(vector<particle> &, vector<group> &, int,...);
  double energy(vector<particle> &, int, vector<short int> &); ///< particle<->list of particles.
  double energy(vector<particle> &);                          ///< all<->all (system energy).
  double energy_dh(vector<particle> &, int);                  ///< all<->particle i.
  double energy_dh(vector<particle> &);                       ///< all<->all (system energy).
  double energy_dh(particle &, particle &);
  double energy_dh(vector<particle> &, group &);              ///< all<->group
  double energy_dh(vector<particle> &, group &, group &);     ///< group<->group
  double internal(vector<particle> &, group &);               ///< internal energy in group
  double internal(vector<particle> &, vector<short int> &);   ///< internal energy of particles in vector
  
  inline double energy_vdwchg(particle &, particle &);        ///< group<->group (elec. + vdW).
  
  double potential(vector<particle> &, point &);              ///< Electrostatic potential in a point
  double quadratic(point &, point &);
  double graft(vector<particle> &, group &);
  double chain(vector<particle> &, group &, int);
  double dipdip(point &, point &, double);                    ///< Dipole-dipole energy.
  double iondip(point &, double, double);                     ///< Ion-dipole energy.
  bool metropolis(double);                                    ///< Metropolis MC test criteria
  bool metropolis(double, double);                            ///< Metropolis MC test for GC tit.
  bool metropolis(double, double, int, double);
};

//! Main particle<->particle interaction (this one NOT in kT)
inline double interact::energy(particle &p1, particle &p2) {
  return p1.charge*p2.charge * frsqrte( p1.sqdist(p2) ); 
  //return p1.charge*p2.charge / p1.dist(p2);
};

inline double interact::energy_dh(particle &p1, particle &p2) {
  double u, qq;
  double r = p1.dist(p2);
  qq=p1.charge*p2.charge;
  u=-vdw /(pow(r,6.)*lB);
  if (qq!=0)
    u+=qq/r*exp(-r*kappa);
  //return p1.charge * p2.charge / r * exp(-kappa*r);
  return u;
};

//! Main particle<->particle interaction (vdW + coulomb, NOT in kT)
inline double interact::energy_vdwchg(particle &p1, particle &p2) {
  double u,qq,r2;
  qq = p1.charge * p2.charge;
  r2 = p1.sqdist(p2);
  u = -vdw / (r2*r2*r2*lB); //div. by lB to match qq/r term below...
  //cout << "!";
  if (qq!=0)
    u+=qq * frsqrte(r2);
  return u;
};

inline double interact::quadratic(point &p1, point &p2) {
  double r=p1.dist(p2)-r0;
  return k*r*r;
};

#endif
