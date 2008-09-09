#ifndef ewald_h
#define ewald_h

/*!
 * \brief Ewald Summation for long-ranged electrostatics
 * \author Martin Trulsson
 * \date 2007
 * \todo Implement in newer faunus framework. I.e. a pairpot(p1,p2) function
 *
 * Energies are calculated using Ewald
 * the class also contains a optimazation concerning
 * cutoff in the Fourier room
 * Cutoff in realspace is set equal to half the box
 */

#include <vector>
#include <cstdarg>
#include <cmath>
#include <complex>
#include "faunus/group.h"
#include "faunus/slump.h"
#include "faunus/particle.h"
#include "faunus/simbox.h"      // <- Should be obselete
#include "faunus/histogram.h"
#include "faunus/interact.h"
#include "faunus/physconst.h"


class Ewald : public Interact {
private:
  static const double a1=0.254829592, a2=-0.284496736;
  static const double a3=1.421413741, a4=-1.453152027;
  static const double a5=1.061405429, p1=0.3275911;
  
public:
  Ewald(int, double=7.12591, int=5);
  Ewald(int, double, double, int=5);   ///< Constructor, sets the Bjerrum length (aangstom)
  
// Ewald parameters
  double ewaldpres;
  double alpha; 
  double alphasqrt;
  double pisqrt;
  double twopi;
  int totk;
  vector<double> kvec;
  vector<complex<double> > eikr;
  vector<complex<double> > eikrold;
  vector< vector< complex<double> > > eix;
  vector< vector< complex<double> > > eiy;
  vector< vector< complex<double> > > eiz;
  vector< vector< complex<double> > > eixold;
  vector< vector< complex<double> > > eiyold;
  vector< vector< complex<double> > > eizold;
  int kmax;
  int ksqmax;

// Ewald Summation routines 
  inline double erfc(double);                                           //< Error function (comblimentary)
  inline double realSpaceEwald(Particle &, Particle &, Simbox &);       //< Real-Space Ewald particle<->particle (NOT in kT)
  
  double realSpaceEwald(vector<Particle> &, int, Simbox &);             //< Real-Space Ewald all<->particle j
  double realSpaceEwald(vector<Particle> &, Simbox &);                  //< Real-Space Ewald all<->all
  void kSpaceEwald(vector<Particle> &, Simbox &);                       //< k-Space Ewald all
  void kSpaceEwald(vector<Particle> &, Simbox &,int);                   //< k-Space Ewald particle i
  double sumkSpaceEwald(vector<Particle> &);
  double sumkSpaceEwald(vector<Particle> &, int);
  double selfEwald(vector<Particle> &);                                 //< Self-Interaction
  void calcAlphaEwald(int, Simbox &, double=0.001);                     //< Optimize alpha
  void initKSpaceEwald(Simbox &);                                       //< initialize k-Space
};
  
inline double Ewald::realSpaceEwald(Particle &p1, Particle &p2, Simbox &s) {
  double r;
  r=p1.sqdist(p2,s);
  if(r > s.sqlen_half) return 0.;
  r = sqrt(r);
  return p1.charge*p2.charge*erfc(r*alphasqrt)/r;
};

inline double Ewald::erfc(double x) {
//Reference for this approximation is found in Abramowitz and Stegun, Handbook of mathematical functions, eq. 7.1.26
  double t,tp,xsq;

  t = 1.0/(1.0+p1*x);
  xsq=x*x;
  tp = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
  return tp*exp(-xsq);
};

#endif
