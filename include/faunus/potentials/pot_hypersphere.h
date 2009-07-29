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
        double d=r*p1.geodesic(p2);
        return LJ(d) + p1.charge*p2.charge / d;  // !!!!!!!!!!
      }

      inline double sqdist(const hyperpoint &p1, const hyperpoint &p2) {
        return p1.hypsqdist(p2)*r;
      }

      inline double elpotential(particle &a, hyperpoint &p) {                                                           
        double chi=a.geodesic(p);                                                                                            
        return a.charge / (pi*r) * ((pi-chi) / tan(chi) -0.5);                                                               
      }    

      inline bool overlap(particle &a, particle &b) {
        double diameter=a.radius+b.radius;
        double delta=diameter/r;
        if(acos(a.hypsqdist(b)) < delta)
          return true;
        return false;
      }

      inline hyperpoint eldipfield(particle a, hyperpoint &p) {
        hyperpoint meldipfield, geodir;
        double chi=a.geodesic(p),
               invsinchi,invsinchisq,cotchi,scaldirpoint,fact;
        fact = 1./(pi*r*r*r)*a.mu;
        invsinchi = 1./sin(chi);
        invsinchisq = invsinchi*invsinchi;
        cotchi = 1./tan(chi);
        scaldirpoint = a.dirmu.hypsqdist(p);
        geodir.z1 = -invsinchi*a.z1+cotchi*p.z1;
        geodir.z2 = -invsinchi*a.z2+cotchi*p.z2;
        geodir.z3 = -invsinchi*a.z3+cotchi*p.z3;
        geodir.z4 = -invsinchi*a.z4+cotchi*p.z4;
        meldipfield.z1 = fact*(2.*invsinchi*scaldirpoint*geodir.z1+invsinchi*(cotchi+(pi-chi)*invsinchisq)*(-a.dirmu.z1+scaldirpoint*p.z1
              +3.*cotchi*scaldirpoint*geodir.z1));
        meldipfield.z2 = fact*(2.*invsinchi*scaldirpoint*geodir.z2+invsinchi*(cotchi+(pi-chi)*invsinchisq)*(-a.dirmu.z2+scaldirpoint*p.z2
              +3.*cotchi*scaldirpoint*geodir.z2));
        meldipfield.z3 = fact*(2.*invsinchi*scaldirpoint*geodir.z3+invsinchi*(cotchi+(pi-chi)*invsinchisq)*(-a.dirmu.z3+scaldirpoint*p.z3
              +3.*cotchi*scaldirpoint*geodir.z3));
        meldipfield.z4 = fact*(2.*invsinchi*scaldirpoint*geodir.z4+invsinchi*(cotchi+(pi-chi)*invsinchisq)*(-a.dirmu.z4+scaldirpoint*p.z4
              +3.*cotchi*scaldirpoint*geodir.z4));
        return meldipfield; // field in point p!!
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
