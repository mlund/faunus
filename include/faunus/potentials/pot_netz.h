#ifndef POT_NETZ_H
#define POT_NETZ_H
#include "faunus/potentials/base.h"
#include "faunus/potentials/pot_coulomb.h"

namespace Faunus {

  /*! \brief Coulomb potential w. an extra empirical PMF
   *  \author Mikael Lund 
   */
  class pot_netz : public pot_coulomb {
    private:
      unsigned char id_I, id_NA, id_CL;

      inline double simple(unsigned char id) {
        if (id==id_I)
          return 8.;
        else if (id==id_CL)
          return .5;
        else
          return 1.;
      }

      inline double sam(double z, unsigned char id) {
        z=0.1*z; // AA->nm
        if (z>1.2) return 0;
        register double A,B,zn,C1,C2,C3,D1,D2,D3;
        if (id==id_I)
          A=50.73,B=18.99,zn=-1.09,C1=-10.22,C2=.54,C3=200,D1=-2.51,D2=.21,D3=60.24;
        else if (id==id_CL)
          A=1.2,B=-1.05,zn=-.99,C1=-7.74,C2=.48,C3=200,D1=-2.26,D2=.12,D3=60.61;
        else if (id==id_NA)
          A=-2.22,B=-7.33,zn=-.89,C1=10.28,C2=.37,C3=100,D1=-.28,D2=.12,D3=2.58;
        else
          return 0.;
        return ( A/pow(z-zn,12)
            - B/pow(z-zn,8)
            + C1*(z-C2)*exp(-C3*pow(z-C2,2) )
            + D1*exp(-D3*pow(z-D2,2) ) ) / f;
      }

      //!< \params z distance from SURFACE.
      //!< \params id ion type
      inline double air(double z, unsigned char id) const { // optimize, please
        z=0.1*z; // AA->nm
        if (z>1.5)
          return 0; // fit not valid beyond 1.5nm
        register double A,B,zn,C1,C2,C3,D1,D2,D3,n;
        if (id==id_I)
          A=.066,B=.977,zn=2.39,C1=-5.1,C2=.7,C3=8.7,D1=-7.32,D2=-.011,D3=2.49,n=1;
        else if (id==id_CL)
          A=15.16,B=4.13,zn=-.37,C1=C2=C3=0,D1=.68,D2=.7,D3=23.99,n=2;
        else if (id==id_NA)
          A=14.62,B=4.35,zn=-.084,C1=4.13,C2=1.1,C3=50.,D1=-.37,D2=.7,D3=10.,n=2;
        else
          return 0.;
        return (  A * ( pow(  exp(-B*(z-zn)) + pow(-1,n)   ,2) - 1. )
            + C1*(z-C2)*exp( -C3*pow(z-C2,2)  )
            + D1*exp( -D3*pow(z-D2,2)  ) ) / f;
      }

    public:
      pot_netz( inputfile &in ) : pot_coulomb(in) {
        name+="/Empirial PMF (AIR)";
        cite="PRL 2007, 99, 226104";
        id_NA=atom["NA"].id;
        id_CL=atom["CL"].id;
        id_I=atom["I"].id;
      }

      inline double hypairpot(const particle &p1, const particle &p2, double r) const {
        if (p1.hydrophobic==true)
          return air(r-p1.radius, p2.id); // c2c -> surface2center
        else if (p2.hydrophobic==true)
          return air(r-p2.radius, p1.id);
        return 0;
      }
  };
}
#endif

