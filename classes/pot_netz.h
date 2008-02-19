#ifndef _POT_NETZ_H
#define _POT_NETZ_H

#include "potentials.h"

/*! \brief Coulomb potential w. an extra empirical PMF
 *  \author Mikael Lund 
 */
class pot_netz : public pot_lj {
  private:
    double A,B,zn,C1,C2,C3,D1,D2,D3,n;
    inline double sam(double z, particle::type) {
      return ( A/pow(z-zn,12)
        - B/pow(z-zn,8)
        + C1*(z-C2)*exp(-C3*pow(z-C2,2) )
        + D1*exp(-D3*pow(z-D2,2) ) ) / f;
    }
    inline double air(double z, particle::type &id) {
      if (id==particle::I)
        A=.066,B=.977,zn=2.39,C1=-5.1,C2=.7,C3=8.7,D1=-7.32,D2=-.011,D3=2.49,n=1;
      if (id==particle::CL)
        A=15.16,B=4.13,zn=-.37,C1=C2=C3=0,D1=.68,D2=.7,D3=23.99,n=2;
      if (id==particle::NA)
        A=14.62,B=4.35,zn=-.084,C1=4.13,C2=1.1,C3=50.,D1=-.37,D2=.7,D3=10.,n=2;
      return (  A * ( pow(  exp(-B*(z-zn)) + pow(-1,n)   ,2) - 1. )
        + C1*(z-C2)*exp( -C3*pow(z-C2,2)  )
        + D1*exp( -D3*pow(z-D2,2)  ) ) / f;
    }
  public:
    double f;             //!< Factor to convert returned energy to kT
    pot_netz( pot_setup &pot ) : pot_lj( pot.eps/pot.lB ) { f=pot.lB; }
    inline double pairpot(particle &p1, particle &p2) {
      register double r2=p1.sqdist(p2), u=lj(p1,p2,r2), r=sqrt(r);
      if (p2.hydrophobic==true)
        if (p1.id==particle::I || p1.id==particle::CL || p1.id==particle::NA)
          u+=air(0.1*r, p1.id); // AA->nm
      if (p1.hydrophobic==true)
        if (p1.id==particle::I || p1.id==particle::CL || p1.id==particle::NA)
          u+=air(0.1*r, p1.id);
      return u + p1.charge*p2.charge/r;
    }
    string info();
};

string pot_netz::info() {
  ostringstream o;
  o << "#   Type               = LJ/Coulomb + empirical PMF" << endl
    << "#   Reference          = PRL (2007),99,226104" << endl
    << "#   Bjerrum length     = " << f << endl
    << "#   LJ epsilon (kT)    = " << eps*f << endl
    << "#   Parameters:"           << endl
    << "#     A,B,zn,n         = " <<  A << "," <<  B << "," << zn << "," << n << endl
    << "#     C1,C2,C3         = " << C1 << "," << C2 << "," << C3 << endl
    << "#     D1,D2,D3         = " << D1 << "," << D2 << "," << D3 << endl;
  return o.str();
}
#endif
