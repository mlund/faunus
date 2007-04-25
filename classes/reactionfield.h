#include "point.h"
#include "legendre.h"
#include "potentials.h"
#include "lennardjones.h"

/*! \brief Charge interactions in- and outside a spherical cavity
 *  \author Mikael Lund
 *  \date 2005-2006 Canberra
 *  \todo A sketch (image) of the system would be nice.
 *
 *  Calculate the pair interaction between charges locate in-
 *  and/or outside a spherical dielectric discontinuity.
 *  See Woodward+Svensson, J.Phys.Chem, 1991, 95, p7471.
 */
class pot_rfield : private pot_lj {
  private:
    int steps;            // # of steps in summation (precision)
    double costheta;      // cos(theta)
    legendre l;           // legendre polynomium 
    double a;             // cavity radius
    double eo, ei;        // dielectric constants (solvent, inside sphere)
  public:
    double f;             //!< Factor to convert to kT

    //! \param pot.epso Dielectric constant outside sphere
    //! \param pot.epsi Dielectric constant inside sphere
    //! \param pot.a Radius of sphere (origo = [0,0,0]).
    //! \param pot.lB Bjerrum length
    //! \param pot.eps L-J parameter
    pot_rfield( pot_setup &pot ) : pot_lj( pot.eps/(pot.epso*pot.lB) )
    {
      f=pot.epso*pot.lB;
      a=pot.a;
      eo=pot.epso;
      ei=pot.epsi;
      steps=50; 
      l.resize(steps);
    };

    inline double pairpot(particle &p1, particle &p2) {
      double r2=p1.sqdist(p2);
      return lj(p1,p2,r2) +
        (p2.charge!=0) ? p2.charge * phi(p1, p2) : 0;
    }

    /*! \brief Self energy (interaction with image charge)
     *  \return \f$ \frac{1}{2} q \phi \f$
     */
    inline double selfenergy(particle &p) {
      return 0.5 * p.charge * phi(p,p,true);
    }

    /*! \brief Calculate potential in point "p" from charge in "p0".
     *  \note Multiply w. lB*eo (=f) to get kT.
     *  \param p0 Charged particle
     *  \param p  Calculate potential in this point
     *  \param self Set to true to calculate self-energy of p0
     */
    double phi(particle &p0, point &p, bool self=false) {
      if (p0.charge==0)
        return 0;
      double sum = 0;
      double r   = p.len();
      double r0  = p0.len();
      double dn;
      costheta=p.dot(p0) / (r * r0);
      l.eval(costheta); 

      //charge outside, potential outside
      if (r0>a && r>a) {
        for (int n=1; n<steps; n+=1) {
          dn   = double(n);
          sum += pow(a*a/(r*r0),n+1) * l.p[n]
            / (  (ei+eo*(1+1/dn))   );
        };
        sum *= (eo-ei) / (eo * a );
        if (self==false)
          sum += 1/(eo*p.dist(p0));
      };

      //charge inside, potential inside
      if (r0<a && r<a) {
        for (int n=0; n<steps; n+=1) {
          dn = double(n);
          sum += pow(r*r0/(a*a),n) * l.p[n]
            / ( eo-ei*(dn/(dn+1)) )  ;
        };
        sum *= (ei-eo) / (ei*a);
        if (self==false)
          sum += 1/(ei*p.dist(p0));
      };

      //charge inside, potential outside
      if (r0<a && r>a) {
        for (int n=0; n<steps; n+=1) {
          dn=double(n);
          sum += (2*dn+1) * pow(r0/r,n) * l.p[n]
            / (  ( dn*ei + eo*(dn+1) )  );
        };
        sum = sum / r;
      };

      //charge outside, potential inside (as above...)
      if (r0>a && r<a) {
        for (int n=1; n<steps; n+=1) {
          dn = double(n);
          sum += pow(r/r0,n) * l.p[n]
            / ( (ei+eo*(1+1/dn))   ) ;
        };
        sum *= (eo-ei) / (eo * r0) ;
        sum += 1/(eo*p.dist(p0));
      };
      return sum * p0.charge;
    }
};
