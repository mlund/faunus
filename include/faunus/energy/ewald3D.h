#ifndef faunus_ewald3D_h
#define faunus_ewald3D_h

/*!
 * \brief Ewald Summation for long-ranged electrostatics
 * \author Martin Trulsson
 * \date 2007
 *
 * Energies are calculated using Ewald
 * the class also contains a optimization concerning
 * cutoff in the Fourier room
 * Cutoff in realspace is set equal to half the box
 */

#include "faunus/common.h"
#include "faunus/energy/base.h"

namespace Faunus {
  class Ewald {
    private:
      static const double a1=0.254829592, a2=-0.284496736;
      static const double a3=1.421413741, a4=-1.453152027;
      static const double a5=1.061405429, p1=0.3275911;
      double boxlen, halfboxlen;
      double lB; // Bjerrum length (AA)

      // Ewald parameters
      double twopi;
      double ewaldpres;
      double alpha; 
      double alphasqrt;
      double pisqrt;
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
      inline double erfc(double);                                 //< Error function (comblimentary)
      void initKSpaceEwald();                                     //< initialize k-Space

    public:
      Ewald(inputfile &);
      string info();

      // Ewald Summation routines 
      inline double realSpaceEwald(particle &, particle &);       //< Real-Space Ewald particle<->particle (NOT in kT)
      double realSpaceEwald(vector<particle> &, int);             //< Real-Space Ewald all<->particle j
      double realSpaceEwald(vector<particle> &);                  //< Real-Space Ewald all<->all
      void kSpaceEwald(vector<particle> &);                       //< k-Space Ewald all
      void kSpaceEwald(vector<particle> &,int);                   //< k-Space Ewald particle i
      double sumkSpaceEwald(vector<particle> &);
      double sumkSpaceEwald(vector<particle> &, int);
      double selfEwald(vector<particle> &);                       //< Self-Interaction
      void calcAlphaEwald(int, double=0.001);                     //< Optimize alpha
  };

  inline double Ewald::realSpaceEwald(particle &p1, particle &p2) {
    double r=p1.sqdist(p2,boxlen,halfboxlen);
    if (r > halfboxlen*halfboxlen)
      return 0.; // redundant?
    r = sqrt(r);
    return p1.charge*p2.charge*erfc(r*alphasqrt)/r;
  }

  /* \note Reference for this approximation is found in Abramowitz and Stegun,
   *       Handbook of mathematical functions, Eq. 7.1.26
   */
  inline double Ewald::erfc(double x) {
    double t = 1./(1.+p1*x),
           tp = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
    return tp*exp(-x*x);
  }
} // namespace
#endif
