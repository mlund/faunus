#ifndef FAU_POT_HYPERSPHERE_H
#define FAU_POT_HYPERSPHERE_H
#include "faunus/potentials/base.h"
#ifdef HYPERSPHERE
namespace Faunus {
  /*! \brief Coulomb/LJ potential on a hypersphere
   *  \author Martin Trulsson
   *  \date Lund, 2009
   *
   *  Detailed description goes here...
   *  \li bjerrum - The Bjerrum length, \f$l_B\f$ [AAngstrom]
   *  \li cellradius - Radius of hypersphere [AAngstrom]
   *  \li LJeps - Lennard-Jones interaction parameter [kT] (see Faunus::pot_lj::eps)
   */
  class pot_hypersphere {
    private:
      double pi;
      string name;
    public:
      double f;       //!< Bjerrum length
      double r;       //!< Cell radius
      double sigma;   //!< LJ sigma
      double epsilon; //!< LJ epsilon (kJ/mol)

      pot_hypersphere(inputfile &in) {
        pi=acos(-1.);
        f=in.getflt("bjerrum",7.1);
        r=in.getflt("cellradius",0);
        sigma=in.getflt("LJsigma",2.8863);
        epsilon=4*in.getflt("LJepsilon",1.97023) / f;
        name="LJ/Coulomb potential on a hypersphere";
      }

      inline double LJ(double chi) {
        double chi2,chi6,chi12;
        //double chi=geodesic(p);
        chi=sigma/(chi*r);
        chi2=chi*chi;
        chi6=chi2*chi2*chi2;
        chi12=chi6*chi6;
        return epsilon*(chi12-chi6);
      }

      /*! \brief Return Coulomb energy between a pair of particles
       *  \return Energy in units of \f$kT/l_B\f$ (lB=f). \f[ \beta u/l_B = \frac{z_1 z_2}{r} + \frac{u_{lj}}{l_B} \f]
       */
      inline double pairpot(const particle &p1, const particle &p2) {
        if (overlap(p1,p2)==true)
          return 1e3;
        double chi = p1.geodesic(p2);
        return p1.charge * p2.charge / (pi*r) * ((pi-chi) / tan(chi) -0.5);
      }

      inline double sqdist(const hyperpoint &p1, const hyperpoint &p2) {
        return p1.hypsqdist(p2)*r;
      }

      inline double elpotential(particle &a, hyperpoint &p) {
        double chi=a.geodesic(p);
        return a.charge / (pi*r) * ((pi-chi) / tan(chi) -0.5);
      }

      inline bool overlap(const particle &a, const particle &b) {
        double diameter=a.radius+b.radius;
        double delta=diameter/r;
        if(acos(a.hypsqdist(b)) < delta)
          return true;
        return false;
      }
 
      string info() {
        std::ostringstream o;
        o << "#   Potential type    = " << name << endl
          << "#   Bjerrum length    = " << f << endl
          << "#   Cell radius       = " << r << endl;
        return o.str();
      }
  };
}//namespace
#endif
#endif
