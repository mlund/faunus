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
#include "lennardjones.h"
#include "inputfile.h"

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
    pot_setup();
    pot_setup(inputfile &);
    double kappa,        //!< Inverse Debye screening length
           lB,           //!< Bjerrum length
           eps,          //!< L-J parameter
           r0,           //!< Bond eq. distance
           epsi,         //!< Internal dielectric constant
           epso,         //!< External dielectric constant
           box,          //!< Cubic box length
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
      register double qq=p1.charge*p2.charge;
      register double u=lj(p1,p2,r2);
      return (qq!=0) ? u+qq/sqrt(r2) : u;
    }
};

/*!
 * \brief Coulomb pot. with minimum image and zero for overlapping spheres.
 * \author Mikael Lund
 * \date 2007
 */
class pot_minimage {
  private:
    double invbox,box;
  public:
    double f;
    pot_minimage(pot_setup &pot) {
      f=pot.lB;
      box=pot.box;
      invbox=1./box;
    }
    inline double pairpot(particle &p1, particle &p2) {
      register double dx,dy,dz,r2;
      dx=p1.x-p2.x;
      dy=p1.y-p2.y;
      dz=p1.z-p2.z;
      dx-=box*floor(dx*invbox+.5);
      dy-=box*floor(dy*invbox+.5);
      dz-=box*floor(dz*invbox+.5);
      r2=dx*dx+dy*dy+dz*dz;
      //r2=p1.sqdist(p2,box,invbox);
      dx=p1.radius+p2.radius;
      return (r2<dx*dx) ? 1e7 : p1.charge*p2.charge/sqrt(r2);
    }
};


class pot_test {
  public:
    double f;
    pot_test( pot_setup &pot ) { f=pot.lB; }
    inline double pairpot(particle &p1, particle &p2) {
      register double r2=p1.sqdist(p2);
      register double a=p1.radius+p2.radius;
      a=a*a;
      if (r2<4*a) {
        a=a/r2;
        a=a*a*a;
        a=a*a/f;
      } else a=0;
      register double qq=p1.charge*p2.charge;
      return (qq!=0) ? qq/sqrt(r2)+a : a;
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
      double r2=p1.sqdist(p2),
             r=sqrt(r2);
      return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
    }
};

/*! \brief Debye-Huckel potential for periodic boundry 
 *         conditions in 3D
 *  \author Mikael Lund/Bjoern Persson
 */
class pot_debyehuckelP3 : private pot_lj {
  private:
    double k;
  public:
    double f;
    double b;
    double inv_b;
    //! \param pot.lB Bjerrum length
    //! \param pot.eps L-J parameter
    //! \param pot.kappa Inverse Debye screening length
    pot_debyehuckelP3( pot_setup &pot ) : pot_lj(pot.eps/pot.lB) {
      f=pot.lB; 
      k=pot.kappa;
      b=pot.box; 
      inv_b=1/pot.box;
    };
    //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
    //! \return Energy in kT/f (f=lB)
    inline double pairpot( particle &p1, particle &p2 ) {
      double r2=p1.sqdist(p2,b,inv_b),
             r=sqrt(r2);
      return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
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
 *  \param T pairpotential
 *
 *  Calculates interaction energies between particles and groups. The
 *  pair potential is specified as a template type which allows inlining.
 *  Unless otherwise specified, all energies will be returned in units of \b kT.
 */

template<class T>
class interaction {
  public:
    T pair;     //!< Pair potential class
    interaction(pot_setup &pot) : pair(pot) {};
    double energy(vector<particle> &, int);                     //!< all<->particle i.
    double energy(vector<particle> &, particle &);              //!< all<->external particle
    double energy(vector<particle> &, group &);                 //!< all<->group.
    double energy(vector<particle> &, group &, group &);        //!< group<->group.
    double energy(vector<particle> &, group &, int);            //!< group<->particle i.
    double energy(vector<particle> &, group &, particle &);     //!< group<->external particle.
    double energy(vector<particle> &, vector<group> &, int,...);
    double energy(vector<particle> &, int, vector<short int> &);//!< particle<->list of particles.
    double energy(vector<particle> &);                          //!< all<->all (System energy).
    double potential(vector<particle> &, unsigned short);       //!< Electric potential at j'th particle
    double internal(vector<particle> &, group &);               //!< internal energy in group
    double pot(vector<particle> &, point &);              //!< Electrostatic potential in a point
    double quadratic(point &, point &);
    double graft(vector<particle> &, group &);
    double chain(vector<particle> &, group &, int);
    double dipdip(point &, point &, double);                    //!< Dipole-dipole energy.
    double iondip(point &, double, double);                     //!< Ion-dipole energy.
};

template<class T>
double interaction<T>::energy(vector<particle> &p, int j) {
  unsigned short ps=p.size();
  double u=0;
  #pragma omp parallel for reduction (+:u)
  for (short i=0; i<j; ++i)
    u+=pair.pairpot( p[i],p[j] );

  #pragma omp parallel for reduction (+:u)
  for (short i=j+1; i<ps; ++i)
    u+=pair.pairpot( p[i],p[j] );
  return pair.f*u;
}

template<class T>
double interaction<T>::energy(vector<particle> &p, group &g) {
  int n=g.end+1, psize=p.size();
  double u=0;
  #pragma omp parallel for reduction (+:u)
  for (int i=g.beg; i<n; ++i) {
    for (int j=0; j<g.beg; j++)
      u += pair.pairpot(p[i],p[j]);
    for (int j=n; j<psize; j++)
      u += pair.pairpot(p[i],p[j]);
  };
  return pair.f*u;
}

template<class T>
double interaction<T>::energy(vector<particle> &p, group &g, int j) {
  double u=0;
  unsigned short len=g.end+1;
  if (g.find(j)==true) {   //avoid self-interaction...
    for (unsigned short i=g.beg; i<j; i++)
      u+=pair.pairpot(p[i],p[j]);
    for (unsigned short i=j+1; i<len; i++)
      u+=pair.pairpot(p[i],p[j]);
  } else                        //simple - j not in g
    for (unsigned short i=g.beg; i<len; i++)
      u+=pair.pairpot(p[i],p[j]);
  return pair.f*u;  
}

template<class T>
double interaction<T>::energy(vector<particle> &p, group &g, particle &a) {
  if (g.beg==-1) return 0;
  double u=0;
  unsigned short i,n=g.end+1;
  for (i=g.beg; i<n; i++)
    u+=pair.pairpot(a, p[i]); 
  return pair.f*u;
}

template<class T>
double interaction<T>::energy(vector<particle> &p) {
  double u=0;
  unsigned short i,j,n = p.size();
  for (i=0; i<n-1; ++i)
    for (j=i+1; j<n; ++j)
      u+=pair.pairpot(p[i], p[j]);
  return pair.f*u; 
}

template<class T>
double interaction<T>::energy(vector<particle> &p, group &g1, group &g2) {
  int ilen=g1.end+1; 
  int jlen=g2.end+1;
  double u=0;
  for (int i=g1.beg; i<ilen; i++)
    for (int j=g2.beg; j<jlen; j++)
      u += pair.pairpot(p[i],p[j]);
  return pair.f*u;
}

/*!
 * ...between the two dipoles a and b, separated by the
 * distance r.
 * \f$ \beta u(r) = l_B \frac{a_x b_x + a_y b_y - 2a_z b_z  }{r^3}\f$
 */
template<class T>
double interaction<T>::dipdip(point &a, point &b, double r) {
  return pair.f*( a.x*b.x + a.y*b.y - 2*a.z*b.z )/(r*r*r);
}
template<class T>
double interaction<T>::iondip(point &a, double q, double r) { return -pair.f*q*a.z/(r*r); }

// Total electrostatic potential in a point
template<class T>
double interaction<T>::pot(vector<particle> &p, point &a) {
  double u=0;
  unsigned short i,n=p.size();  
  for (i=0; i<n; i++) u+=p[i].charge/p[i].dist(a);
  return pair.f*u;
}

// Internal (NON)-electrostatic energy in group
template<class T>
double interaction<T>::internal(vector<particle> &p, group &g) {
  if (g.beg==-1) return 0;
  double u=0;
  unsigned short glen=g.end+1;
  for (unsigned short i=g.beg; i<glen-1; i++)
    for (unsigned short j=i+1; j<glen; j++)
      u+=pair.pairpot(p[i],p[j]);
  return pair.f*u;
}

template<class T>
double interaction<T>::energy(vector<particle> &p, particle &a) {
  double u=0;
  unsigned short n=p.size();
  for (unsigned short i=0; i<n; i++)
    u+=pair.pairpot(p[i], a);
  return pair.f*u;
}

/*! \note If the charge of the j'th particle is 0, ZERO will be returned!
 *  \return \f$ \phi_j = \sum_{i\neq j}^{N} \frac{l_B z_i}{r_{ij}} \f$
 *  \param j The electric potential will be calculated in the point of this particle
 */
template<class T>
double interaction<T>::potential(vector<particle> &p, unsigned short j) {
  if (p[j].charge==0) return 0;
  double u=0;
  unsigned short i,n=p.size();
  for (i=0; i<j; ++i) u+=p[i].charge/p[i].dist(p[j]);
  for (i=j+1; i<n; ++i) u+=p[i].charge/p[i].dist(p[j]);
  return u;
}

#endif

