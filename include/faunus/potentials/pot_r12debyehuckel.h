#ifndef FAU_POT_R12DEBYEHUCKEL_H
#define FAU_POT_R12DEBYEHUCKEL_H
#include "faunus/potentials/base.h"
#include "faunus/legendre.h"
//#include "faunus/fortran.h"

namespace Faunus {
  /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in 3D with r12 repulsion 
   *         usuable in cuboid containers
   *  \author Chris Evers / Mikael Lund
   *  \date Lund February 2011
   */
  class pot_r12debyehuckel : public pot_harmonic {
    private:
    protected:
      point len_inv;                //!< inverse side lengths
      bool setlen(point l) {
        assert(l.x>0);              // debug information
        assert(l.y>0);              // 
        assert(l.z>0);              // 
        if (l.x<=0  ||              // check non-negative value
            l.y<=0  ||              // 
            l.z<=0  )               // 
          return false;             // 
        len = l;                    // cuboid sidelength
        len_half=l*0.5;             // half cuboid sidelength
        len_inv.x=1./len.x;         // inverse cuboid side length
        len_inv.y=1./len.y;         // 
        len_inv.z=1./len.z;         // 
        return true;
      }
    public:
      point len;                    //!< side lengths
      point len_half;               //!< half side lengths
      double I;                     //!< ionicstrenght (M)
      double k;                     //!< inverse debye length
      double eps;                   //!< interaction parameter
      double f;                     //!< bjerrum length
      string name;                  //!< potential name
      pot_r12debyehuckel( inputfile &in ) : pot_harmonic(in) {
        name="r12 + Debye-Huckel w. minimum image"; 
        double cubelen=in.getflt("cuboid_len",-1);  // Read side lengths
        if (cubelen<=0) {
          len.x=in.getflt("cuboid_xlen",0);
          len.y=in.getflt("cuboid_ylen",0);
          len.z=in.getflt("cuboid_zlen",0);
          if (len.x<=0 || len.y<=0 || len.z<=0)
            cubelen=in.getflt("boxlen",-1);
        } else {
          len.x=cubelen;
          len.y=cubelen;
          len.z=cubelen;
        }
        setlen(len);

        f=in.getflt("bjerrum",7.1);                 // Bjerrum length
        k=1/in.getflt("debyelen",1.1e6);            // Inverse debye length k_D
        if ( 1/k>=1e6) {                            //   not read?
          I=in.getflt("ionicstr",0);                //   read ionicstrenght instead
          k=sqrt( 8*pyc.pi*f*pyc.Nav/1e27*I );      //   k_D=\sqrt{8 \pi \lambda_B[A] I[A-3]}
        }
        eps=in.getflt("lj_epsilon", 0.2 );          // Interaction parameter
        eps*=4;                                     //! do we want this?
        eps=eps/f;                                  //! and this?
      }

      //! Define Debye-Huckel energy function
      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + eps \frac{sigma^12}{r^12}/f \f$
      //! \return Energy in units of kT/lB
      inline double pairpot( const particle &p1, const particle &p2 ) {
        double r2=sqdist(p1,p2),              // squared distance between particles
               r=sqrt(r2),                    // distance between particles
               s=p1.radius+p2.radius,         // distance between particles at contact
               u=s*s/r2;                      // (s/r)^2
        u=u*u*u;                              // (s/r)^6
        return p1.charge * p2.charge / r * exp(-k*r) + eps*(u*u);
      }

      //! Calculate distance using the minimum image convention
      inline double sqdist(const point &p1, const point &p2) {   //!< Squared distance 
        double dx=std::abs(p1.x-p2.x),
               dy=std::abs(p1.y-p2.y),
               dz=std::abs(p1.z-p2.z);
        if (dx>len_half.x) dx-=len.x;
        if (dy>len_half.y) dy-=len.y;
        if (dz>len_half.z) dz-=len.z;
        return dx*dx + dy*dy + dz*dz;
      }

      void setvolume(double vol) {
        len.x=pow(vol, 1./3.);
        len.y=pow(vol, 1./3.);
        len.z=pow(vol, 1./3.);
        len_half=len*0.5;
      }

      inline int anint(double x) const {
        return int(x>0. ? x+.5 : x-.5);
      }

      //! Apply periodic boundary conditions
      inline void boundary(point &a) const {
        if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
        if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
        if (std::abs(a.z)>len_half.z) a.z-=len.z*anint(a.z*len_inv.z);
      }

      virtual double energy(const vector<particle> &p, const vector<int> &pairs) {
        int n=pairs.size()/2, l=0;
        double u=0;
        for (int i=0; i<n; i++)
          u+=pairpot(p[ pairs[l++] ], p[ pairs[l++]] );
        return f*u;
      }

      string info() {
        std::ostringstream o;
        o << "#     Name              = " << name  << endl
          << "#     Bjerrum length    = " << f     << endl
          << "#     Debye length      = " << 1./k  << endl
          << "#     Ionic strength    = " << k*k*1e30/(8*pyc.pi*f*pyc.Nav) << " mM" << endl
          << pot_harmonic::info();
        return o.str();
      }

      /*!
       * Calculated the bonded + non-bonded energy between two particles
       * connected with a harmonic bond. To avoid numerical trouble with
       * the non-electrostatic, non-bonded interactions, the radii of
       * the two particles are temporarily reduced.
       */
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
  };

  /*! \brief Tabulated Debye-Huckel potential with r12 repulsion 
   *         usuable in cuboid containers
   *  \author Chris Evers
   *  \date Kastrup, March 2011
   */
  class pot_r12debyehuckel_tab : public pot_r12debyehuckel {
  private:
    xytable< double, double > U;            //!< Debye-Huckel potential at r [kT/z^2]
    double dr2;                             //!< Table resolution
    double r2min;                           //!< Minimum tabulated potential, below exact
    double r2max;                           //!< Maximum tabulated potential, higher zero

  public:
    pot_r12debyehuckel_tab( inputfile &in ) : pot_r12debyehuckel(in) {
      name += ", interpolated table"; 
      dr2 = in.getflt("tabpot_dr2", 0.1);
      double invdr2 = 1/dr2;
      double Umax=in.getflt("tabpot_Umax", .1),
             Umin=in.getflt("tabpot_Umin", 1e-5);
      for (double r2=0; r2<=len.z*len.z; r2+=dr2) {
        double Utest=calcPotential(r2);
        if ( Utest > Umax )
          r2min = r2;
        if ( Utest < Umin ) {
          r2max = r2;
          break;
        }
      }
      U.init(dr2, r2min-dr2, r2max+dr2);
      for (double r2=r2min-dr2; r2<=r2max+dr2; r2+=dr2) {
        U(r2)=calcPotential(r2);
      }
    }

    virtual double calcPotential(double r2) { 
      double r=sqrt(r2);                    // distance between particles
      return 1 / r * exp(-k*r);
    }

    //! Define Debye-Huckel energy function
    //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + eps \frac{sigma^12}{r^12}/f \f$
    //! \return Energy in units of kT/lB
    inline double pairpot( const particle &p1, const particle &p2 ) {
      double r2=sqdist(p1,p2);              // squared distance between particles
      if ( r2 >= r2max )
        return 0;
      double s=p1.radius+p2.radius,         // distance between particles at contact
             u=s*s/r2;                      // (s/r)^2
      u=u*u*u;                              // (s/r)^6
      if ( r2 <= r2min ) {
        double r=sqrt(r2);
        return p1.charge * p2.charge / r * exp(-k*r)  + eps*(u*u);
      }
      else
        return p1.charge * p2.charge * getPotential(r2) + eps*(u*u);
    }

    //! \return Potential in units of kT / (lB z^2)
    double getPotential(const double &r2) { 
      return U.interpolate(r2); 
    }

    string info() {
      std::ostringstream o;
      o << pot_r12debyehuckel::info()
        << "#     Table resolution  = " << sqrt(dr2)   << " AA ( " << U.y.size() << " slits )" << endl
        << "#     Exact pot cut off = " << sqrt(r2min) << " AA ( " << calcPotential(r2min) << " kT/z^2 )" << endl
        << "#     Tab pot cut off   = " << sqrt(r2max) << " AA ( " << calcPotential(r2max) << " kT/z^2 )" << endl;
      return o.str();
    }

    void dump(string filename) {
      U.dumptodisk(filename);
    }
  };

  /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in XY only with r12 repulsion 
   *         usuable in cuboid containers
   *  \author Chris Evers
   *  \date Lund February 2011
   */
  class pot_r12debyehuckelXY : public pot_r12debyehuckel {
    public:
      pot_r12debyehuckelXY( inputfile &in ) : pot_r12debyehuckel(in) {
      }

      //! Calculate distance using the minimum image convention
      inline double sqdist(const point &p1, const point &p2) {   //!< Squared distance 
        double dx=std::abs(p1.x-p2.x),
               dy=std::abs(p1.y-p2.y),
               dz=p1.z-p2.z;
        if (dx>len_half.x) dx-=len.x;
        if (dy>len_half.y) dy-=len.y;
        return dx*dx + dy*dy + dz*dz;
      }

      //! Apply periodic boundary conditions
     inline void boundary(point &a) const {
        if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
        if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
      }

      string info() {
        std::ostringstream o;
        o << pot_r12debyehuckel::info()
          << "#     Periodicity       = slit: xy periodicity only" << endl;
        return o.str();
      }
  };

    /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in XY only with r12 repulsion 
   *         and a tabulated potential, usuable in cuboid 
   *         containers
   *  \author Chris Evers
   *  \date Vierlingsbeek March 2011
   */
  class pot_r12debyehuckelXY_tab : public pot_r12debyehuckel_tab {
    public:
      pot_r12debyehuckelXY_tab( inputfile &in ) : pot_r12debyehuckel_tab(in) {
      }

      //! Calculate distance using the minimum image convention
      inline double sqdist(const point &p1, const point &p2) {   //!< Squared distance 
        double dx=std::abs(p1.x-p2.x),
               dy=std::abs(p1.y-p2.y),
               dz=p1.z-p2.z;
        if (dx>len_half.x) dx-=len.x;
        if (dy>len_half.y) dy-=len.y;
        return dx*dx + dy*dy + dz*dz;
      }

      //! Apply periodic boundary conditions
     inline void boundary(point &a) const {
        if (std::abs(a.x)>len_half.x) a.x-=len.x*anint(a.x*len_inv.x);
        if (std::abs(a.y)>len_half.y) a.y-=len.y*anint(a.y*len_inv.y);
      }

      string info() {
        std::ostringstream o;
        o << pot_r12debyehuckel_tab::info()
          << "#     Periodicity       = slit: xy periodicity only" << endl;
        return o.str();
      }
  };
}
#endif
