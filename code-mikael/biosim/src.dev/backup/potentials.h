#ifndef _potential_h
#define _potential_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "species.h"
#include "xydata.h"
#include "group.h"
#include "slump.h"

/*!
 *  \brief General class for potentials.
 *  \author Mikael Lund
 *  \date Prague, 2007
 *  
 *  This class is inherited by other potential
 *  classes for calculating energies between
 *  particle pairs.
 */
class potential : private slump {
 public:
   double kappa;   //!< Debye screening length 
   double r0,k;    // temporary
   double finv, f; //!< Factor to convert energy to kT.

   //! \param x Factor to convert energy to kT.
   potential(double x) { f=x; finv=1./f; }
   //! Convert energy to kT
   double tokT(double &u) { return f*u; }

   //! Test Metropolis criteria (NVT)
   bool metropolis(double du) {
     if (du > 0)
       if ( random_one()>exp(-du) )
         return false;
     return true;
   }
};

/*!
 *  \brief Lennard-Jones potential
 *  \author Mikael Lund
 *  \year Prague, 2007
 */
class lennardjones {
  private:
    double eps;
  public:
    lennardjones(double epsilon) { eps=epsilon; }
    /*!
     *  L-J pair energy.
     *  \f$ u_{LJ} = \epsilon \left ( \frac{\sigma}{r^{12}} - \frac{\sigma}{r^6} \right ) \f$
     *  (in units of kT/f).
     *  \param r2 Squared distance between particle 1 and 2.
     */
    inline double lj(particle &p1, particle &p2, double &r2) {
      double x=p1.radius+p2.radius, u=x*x/r2;
      x=u*u*u;
      return (x*x-x)*eps;
    }
};

/*!
 *  \brief Coulomb potential
 *  \author Mikael Lund
 *  \date Prague, 2007
 */
class pot_coulomb : public potential, private lennardjones {
  public:
    /*! \param f Bjerrum length
     *  \param eps L-J epsilon parameter */
    pot_coulomb ( double f, double eps=0.2 ) : potential(f=7.13), lennardjones(eps/f) {}

    //! \brief Return Coulomb energy between a pair of particles
    //!
    //! \f$ u = \frac{l_B z_1 z_2}{r} + u_{LJ}(r) \f$ (in units of kT/f)
    inline double pairpot(particle &p1, particle &p2) {
      double r2=p1.sqdist(p2);
      double qq=p1.charge*p2.charge;
      double u=lj(p1,p2,r2);
      return (qq!=0) ? u+qq/sqrt(r2) : u;
    }
};

// Debye-Huckel potential
class pot_debyehuckel : public potential, private lennardjones {
  public:
    pot_debyehuckel( double f, double eps ) : 
      potential(f=7.13), lennardjones(eps/f) {}
    inline double pairpot( particle &p1, particle &p2 ) {
      double r2=p1.sqdist(p2), r=sqrt(r2);
      return lj(p1,p2,r2) + p1.charge*p1.charge/r*exp(-kappa*r);
    }
};

/*!
 * \brief Load pair potentials from disk
 * \author Mikael Lund
 * \date Prague, 2007
 */
class pot_datapmf : public potential, private lennardjones {
  public:
    pot_datapmf(double f, double eps) : potential(f=7.13), lennardjones(eps/f) {}
    xydata<double> pmfd[particle::LAST][particle::LAST];
    bool loadpmf(species &, string);          // load pmf's from disk
    void loadpmf(species &, string,string);   // -//-
    void showpmf(species &);                  // show inter-species pmf's.
    double pairpot (particle &p1, particle &p2) {
      unsigned short i=p1.id,j=p2.id;
      if (i>j) swap(i,j);
      double r2=p1.sqdist(p2);
      if (pmfd[i][j].xmax==0) {               // if no data exists:
        double u,x=p1.charge*p2.charge;       // use Coulomb + lj pot.
        u=lj(p1,p2,r2);
        return (x!=0) ? u+x/sqrt(r2) : u; 
      }; 
      r2=sqrt(r2);
      return (r2>pmfd[i][j].xmax) ? 
        p1.charge*p2.charge/r2 : // use Coulomb pot. outside data 
        pmfd[i][j].x2y(r2);      // ...else use table. 
    }
};

// The pair potential is specified using
// macro definitions at compile time (-D option).
#ifdef POT_COULOMB
#define PAIRPOTENTIAL pot_coulomb     /* plain coulomb potential */
#endif
#ifdef POT_DATAPMF
#define PAIRPOTENTIAL pot_datapmf     /* potential from disk */
#endif
#ifdef POT_DEBYEHUCKEL
#define PAIRPOTENTIAL pot_debyehuckel /* D-H screened potential */
#endif

/*!
 *  \brief Interaction between particles and groups
 *  \author Mikael Lund
 *
 *  Calculates interaction energies between particles and groups. The
 *  pair potential is specified at compile time using the following
 *  macro definitions which is specified by the compiler "-D" option.\n
 *   - POT_COULOMB (Coulomb + LJ)
 *   - POT_DEBYEHUCKEL (Debye-Huckel + LJ)
 *   - POT_DATAPMF (Load pair potential from disk)
 *
 *  Unless other wise specified, all energies will be returned in units of \b kT.
 */
class interaction : public PAIRPOTENTIAL
{
  public:
    interaction(double f, double opt=1) : PAIRPOTENTIAL(f,opt) {};
    double energy(vector<particle> &, int);                     ///< all<->particle i.
    double energy(vector<particle> &, group &);                 ///< all<->group.
    double energy(vector<particle> &, group &, group &);        ///< group<->group.
    double energy(vector<particle> &, group &, int);            ///< group<->particle i.
    double energy(vector<particle> &, group &, particle &);     ///< group<->external particle.
    double energy(vector<particle> &, vector<group> &, int,...);
    double energy(vector<particle> &, int, vector<short int> &); ///< particle<->list of particles.
    double energy(vector<particle> &);                          ///< all<->all (system energy).
    double internal(vector<particle> &, group &);               ///< internal energy in group

    double pot(vector<particle> &, point &);              ///< Electrostatic potential in a point
    double quadratic(point &, point &);
    double graft(vector<particle> &, group &);
    double chain(vector<particle> &, group &, int);
    double dipdip(point &, point &, double);                    ///< Dipole-dipole energy.
    double iondip(point &, double, double);                     ///< Ion-dipole energy.
};

inline double interaction::quadratic(point &p1, point &p2) {
  double r=p1.dist(p2)-r0;
  return k*r*r;
}

class hardsphere {
  public:
    bool overlap(vector<particle> &, int);                              //< all<->particle i.
    bool overlap(vector<particle> &, particle &);                       ///< all<->arbitrary (external) particle.
    bool overlap(vector<particle> &, group::group &);                   ///< all<->group.
    bool overlap(vector<particle> &, group::group &, int);              ///< group<->particle i.
    bool overlap(vector<particle> &, group::group &, particle &);       ///< group<>arbitrary (external) particle
    bool overlap(vector<particle> &, group::group &, group::group &);   ///< group<->group.
    bool overlap(vector<particle> &, vector<short int> &, double);      ///< internal collisions within subset
    bool celloverlap(vector<particle> &, group::group &, double);       ///< group with a spherical cell
    bool chgoverlap(vector<particle> &, group::group &, double);        ///< charge overlap within group
};

#endif

