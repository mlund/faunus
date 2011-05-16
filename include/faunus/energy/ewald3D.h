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
      static const double a1, a2, a3, a4, a5, p1;
      double boxlen, halfboxlen;
      double lB; // Bjerrum length (AA)

      // Ewald parameters
      double twopi;
      double ewaldpres;
      double alpha; 
      double alphasqrt;
      double pisqrt;
      double Ewapre;
      int totk;
      vector<double> kvec;
      vector<complex<double> > eikr, eikrold;
      vector< vector< complex<double> > > eix, eiy, eiz, eixold, eiyold, eizold;
      int kmax, ksqmax;
      inline double erfc(double);                                 //< Error function (comblimentary)
      void initKSpaceEwald();                                     //< initialize k-Space

    public:
      double f;                                                   //!< Factor to convert to kT
      Ewald(inputfile &);
      string info();

      // Ewald Summation routines 
      void setvolume(double);
      inline double pairpot(const particle &, const particle &);  //< particle<->particle (real space)
      inline double sqdist(const point &, const point &);

      void kSpaceEwald(vector<particle> &);                       //< k-Space Ewald all
      void kSpaceEwald(vector<particle> &,int);                   //< k-Space Ewald particle i
      double sumkSpaceEwald(vector<particle> &);
      double sumkSpaceEwald(vector<particle> &, int);
      double selfEwald(vector<particle> &);                       //< Self-Interaction
      void calcAlphaEwald(int, double=0.001);                     //< Optimize alpha
  };

  void Ewald::setvolume(double vol) {
    boxlen=pow(vol,1/3.);
    halfboxlen=boxlen/2;
    initKSpaceEwald();
  }

  inline double Ewald::sqdist(const point &p1, const point &p2) {
    return p1.sqdist_mi_xyz(p2,boxlen,halfboxlen);
  }

  inline double Ewald::pairpot(const particle &p1, const particle &p2) {
    double r=p1.sqdist_mi_xyz(p2,boxlen,halfboxlen);
    if (r > halfboxlen*halfboxlen)
      return 0.;
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

  /*
   * This is the actual interaction class used by MC move classes (these cannot be modified!).
   *
   */
  class interaction_ewald : public interaction<Ewald> {
    using interaction<Ewald>::pair;
    public:
      interaction_ewald(inputfile &in, container &con) : interaction<Ewald>(in) {
        interaction<Ewald>::name+=" with ewald long range correction";
        pair.selfEwald(con.p); // update for each volume update (or when alpha changes)
      }

      double energy(vector<particle> &p) {
        pair.kSpaceEwald(p); //?
        return interaction<Ewald>::energy(p) + pair.sumkSpaceEwald(p);
      }

      // should return TOTAL interaction of i'th particle with the rest of the system (real and imaginary)
      double energy(vector<particle> &p, int i) {
        pair.kSpaceEwald(p,i);
        return interaction<Ewald>::energy(p,i) + pair.sumkSpaceEwald(p,i); 
      }

      double energy(vector<particle> &p, const group &g) {
        for (int i=g.beg; i<=g.end; i++)
          pair.kSpaceEwald( p, i );
        double u=interaction<Ewald>::energy(p,g); // real space
        for (int i=g.beg; i<=g.end; i++)
          u+=pair.sumkSpaceEwald( p, i );         // kspace
        return u;
      }

      double energy(vector<particle> &p, const particle &a) {
        return interaction<Ewald>::energy(p,a) + 0; // kspace energy function missing for external particle
      }

  };
} // namespace
#endif
