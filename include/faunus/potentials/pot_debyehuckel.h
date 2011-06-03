#ifndef FAU_POT_DEBYEHUCKEL_H
#define FAU_POT_DEBYEHUCKEL_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Debye-Huckel potential
   *  \author Mikael Lund
   *  
   *  Faunus::inputfile is scanned for "bjerrum", "debyelen", "LJeps".
   */
  class pot_debyehuckel : public pot_lj {
    public:
      double I,k;
      pot_debyehuckel( inputfile &in ) : pot_lj(in) {
        name+="/Debye-Huckel";
        f=in.getflt("bjerrum",7.1);
        k=1./in.getflt("debyelen",1.1e4);
        if ( 1/k>=1e4) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 4*std::acos(-1.)*f*6.022e23/1e27*2*I );
        }
        eps=eps/f;
      }
      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in kT/f (f=lB)
      inline double pairpot( const particle &p1, const particle &p2 ) {
        double r2=p1.sqdist(p2),
               r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }
      
      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2);
      }
      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f << endl
          << "#   Kappa             = " << k << endl
          << "#   Debye length      = " << 1./k << endl
          << "#   Ionic strength (M)= " << k*k*1e27/(8*std::acos(-1.)*f*6.022e23) << endl;
        return o.str();
      }
  };
  /*!
   * \brief Screened Columb/LJ potential with minimum image in XY, hard in Z.
   * \author Mikael Lund
   * \date 2009
   */
  class pot_debyehuckelXY : public pot_lj {
    private:
      double halfbox,box;
    public:
      string name;
      double k,I;
      pot_debyehuckelXY( inputfile &in ) : pot_lj(in) {
        name+="Screened Columb/LJ w. minimum image (XY, only)";
        f=in.getflt("bjerrum",7.1);
        box=in.getflt("boxlen");
        halfbox=box/2.;
        k=1./in.getflt("debyelen",1.1e4);
        if ( 1/k>=1e4) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 4*std::acos(-1.)*f*6.022e23/1e27*2*I );
        }
        eps=eps/f;
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);
        halfbox=box/2.;
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        double r2=sqdist(p1, p2),
               r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }
      inline double sqdist(const point &p1, const point &p2) {
        double dz=p1.z-p2.z;
        double dx=std::abs(p1.x-p2.x);
        double dy=std::abs(p1.y-p2.y);
        if (dx>halfbox) dx-=box;
        if (dy>halfbox) dy-=box;
        return dx*dx + dy*dy + dz*dz;
      }
      string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   Bjerrum length    = " << f << std::endl
          << "#   Screening length  = " << 1./k << std::endl
          << "#   4*LJ epsilon      = " <<f*eps<<std::endl
          << "#   Image length      = " << box << std::endl;
        return o.str();
      }
  };
  class pot_debyehuckelXYcc : public pot_lj, public pot_harmonic {
    private:
      double halfbox,box;
      double zbox, zhalfbox;
    public:
      string name;
      double k,I,c,c2,rho2,scd;
      double hphoba,hphobr, hphobmax;
      pot_debyehuckelXYcc( inputfile &in ) : pot_lj(in), pot_harmonic(in) {
        name+="Screened Columb/LJ w. minimum image (XY, only) and cylindrical cutoff\n#            external correction and hyrdophobic interaction ";
        f=in.getflt("bjerrum",7.1);
        box=in.getflt("boxlen");
        halfbox=box/2.;
        zbox=in.getflt("zboxlen");
        zhalfbox=zbox*0.5;
        c=in.getflt("c_cut", halfbox);
        c2=c*c;
        k=1./in.getflt("debyelen",1.1e4);
        if ( 1/k>=1e4) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 4*std::acos(-1.)*f*6.022e23/1e27*2*I );
        }
        scd=0;
        eps=eps/f;
        hphoba=in.getflt("hydrophobic_amp");
        hphobr=in.getflt("hydrophobic_range");
      }
      void setchargedens(double s) {
        scd=s*std::acos(-1.)*2*f/k;
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);
        halfbox=box/2.;
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        double r2=sqdist(p1, p2, rho2);
        if (rho2>c2)
          return 0;
        double r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }
      inline double expot(const particle &p) {  //Returns interaction in kT!
        double z=p.z+zhalfbox;
        return p.charge*scd*exp(-k*sqrt(c2+z*z));
      }
      inline double hphobpot(const particle &p) { //Returns interaction in kT!
        if(p.hydrophobic==true)
          return (p.z+zhalfbox<p.radius*hphobr) ? hphoba : 0.; 
        return 0.;
      }
      inline double sqdist(const point &p1, const point &p2, double rho2) {
        double dz=p1.z-p2.z;
        double dx=std::abs(p1.x-p2.x);
        double dy=std::abs(p1.y-p2.y);
        if (dx>halfbox) dx-=box;
        if (dy>halfbox) dy-=box;
        rho2=dx*dx+dy*dy;
        return rho2 + dz*dz;
      }
      inline double sqdist(const point &p1, const point &p2) {
        double dz=p1.z-p2.z;
        double dx=std::abs(p1.x-p2.x);
        double dy=std::abs(p1.y-p2.y);
        if (dx>halfbox) dx-=box;
        if (dy>halfbox) dy-=box;
        return dx*dx + dy*dy + dz*dz;
      }
      inline double bond(particle &p1, particle &p2) {
        double r=sqrt(sqdist(p1,p2));
        double u=harmonicbond(p1,p2,r)/f;
        p1.radius*=0.5;
        p2.radius*=0.5;
        u+=pairpot(p1,p2);
        p1.radius*=2;
        p2.radius*=2;
        return u;
      }
      string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   Cutoff            = " << c<<std::endl
          << "#   Surf. Charge Dens.= " <<scd*k/f/2/std::acos(-1.)<<std::endl
          << "#   Hydrophobic amp.  = " <<hphoba<<std::endl
          << "#   Hydrophobic range = " <<hphobr<<" (sigma/2)"<<std::endl
          << "#   Bjerrum length    = " << f << std::endl
          << "#   Screening length  = " << 1./k << std::endl
          << "#   4*LJ epsilon      = " <<f*eps<<std::endl
          << "#   Image length      = " << box << std::endl
          << pot_harmonic::info();
        return o.str();
      }
  };
}
#endif
