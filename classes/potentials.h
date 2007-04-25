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
#include "lennardjones.h"

/*!
 *  \brief Setup for potentials.
 *  \author Mikael Lund
 *  \date Prague, 2007
 *
 *  This class is used to pass parameters to classes
 *  that handles particle pair-potentials.
 */
class pot_setup {
 public:
   double kappa,        //!< Inverse Debye screening length
          lB,           //!< Bjerrum length
          eps,          //!< L-J parameter
          r0,           //!< Bond eq. distance
          epsi,         //!< Internal dielectric constant
          epso,         //!< External dielectric constant
          a;            //!< Cavity radius
};

/*!
 *  \brief Coulomb potential
 *  \author Mikael Lund
 *  \date Prague, 2007
 */
class pot_coulomb : private pot_lj {
  public:
    double f; //!< Factor to convert kT/lB to kT (here the Bjerrum length).

    /*! \param pot.lB Bjerrum length
     *  \param pot.eps L-J epsilon parameter */
    pot_coulomb ( pot_setup &pot) : pot_lj(pot.eps/pot.lB) { f=pot.lB; };

    /*! \brief Return Coulomb energy between a pair of particles
     *  \return Energy in units of kT/f (f=lB).
     *
     *  \f$ \beta u/f = \frac{z_1 z_2}{r} + u_{LJ}/f \f$
     */
    inline double pairpot(particle &p1, particle &p2) {
      double r2=p1.sqdist(p2);
      double qq=p1.charge*p2.charge;
      double u=lj(p1,p2,r2);
      return (qq!=0) ? u+qq/sqrt(r2) : u;
    }
};

/*! \brief Debye-Huckel potential
 *  \author Mikael Lund
 */
class pot_debyehuckel : private pot_lj {
  private:
    double k;
  public:
    double f;
    //! \param pot.lB Bjerrum length
    //! \param pot.eps L-J parameter
    //! \param pot.kappa Inverse Debye screening length
    pot_debyehuckel( pot_setup &pot ) : pot_lj(pot.eps/pot.lB) {
      f=pot.lB; 
      k=pot.kappa; 
    };
    //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
    //! \return Energy in kT/f (f=lB)
    inline double pairpot( particle &p1, particle &p2 ) {
      double r2=p1.sqdist(p2), r=sqrt(r2);
      return lj(p1,p2,r2) + p1.charge*p1.charge/r*exp(-k*r);
    }
};

/*!
 * \brief Load pair potentials from disk
 * \author Mikael Lund
 * \date Prague, 2007
 *
 * Class to load potential of mean force (PMF) data from file(s). If
 * the distance between the particles is outside the data range a 
 * simple Coulomb + LJ potential will be applied.
 */
class pot_datapmf : private pot_lj {
  private:
    xydata<double> pmfd[particle::LAST][particle::LAST];
  public:
    double f;
    //! \param pot.lB Bjerrum length
    //! \param pot.eps L-J parameter
    pot_datapmf(pot_setup &pot) : pot_lj(pot.eps/pot.lB) { f=pot.lB; }
    bool loadpmf(species &, string);          // load pmf's from disk
    void loadpmf(species &, string,string);   // -//-
    void showpmf(species &);                  //!< Lists loaded pmf's.
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
 *   - POT_COULOMB (Coulomb + LJ) - see \link pot_coulomb \endlink.
 *   - POT_DEBYEHUCKEL (Debye-Huckel + LJ) - see \link pot_debyehuckel \endlink.
 *   - POT_DATAPMF (Load pair potential from disk) - see \link pot_datapmf \endlink.
 *
 *  Unless otherwise specified, all energies will be returned in units of \b kT.
 */
class interaction : public PAIRPOTENTIAL, private slump
{
  public:
    interaction(pot_setup &pot) : PAIRPOTENTIAL(pot) {};
    double energy(vector<particle> &, int);                     ///< all<->particle i.
    double energy(vector<particle> &, group &);                 ///< all<->group.
    double energy(vector<particle> &, group &, group &);        ///< group<->group.
    double energy(vector<particle> &, group &, int);            ///< group<->particle i.
    double energy(vector<particle> &, group &, particle &);     ///< group<->external particle.
    double energy(vector<particle> &, vector<group> &, int,...);
    double energy(vector<particle> &, int, vector<short int> &);///< particle<->list of particles.
    double energy(vector<particle> &);                          ///< all<->all (System energy).
    double internal(vector<particle> &, group &);               ///< internal energy in group

    double pot(vector<particle> &, point &);              ///< Electrostatic potential in a point
    double quadratic(point &, point &);
    double graft(vector<particle> &, group &);
    double chain(vector<particle> &, group &, int);
    double dipdip(point &, point &, double);                    ///< Dipole-dipole energy.
    double iondip(point &, double, double);                     ///< Ion-dipole energy.

    //! Test Metropolis criteria (NVT ensemble)
    bool metropolis(double du) {
      if (du > 0)
        if ( random_one()>exp(-du) )
          return false;
      return true;
    }
};

inline double interaction::quadratic(point &p1, point &p2) {
  double r0,k;
  double r=p1.dist(p2)-r0;
  return k*r*r;
}

/*!\brief Hardsphere overlap between particles
 * \author Mikael Lund
 */
class hardsphere {
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

#endif

