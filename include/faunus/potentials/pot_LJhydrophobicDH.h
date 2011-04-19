#ifndef FAU_POT_LJHYDROPHOBICDH_H
#define FAU_POT_LJHYDROPHOBICDH_H
#include "faunus/potentials/base.h"
#include "faunus/legendre.h"
//#include "faunus/fortran.h"

#ifdef FAU_FAST_EXP
#define _EXPFUNC(x) f_exp2(x)
#else
#define _EXPFUNC(x) std::exp(x)
#endif

namespace Faunus {
  /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in 3D, it is extended to preform 
   *         under conditions of constant pressure.
   *         See class isobaric->markovmove.h
   *  \author Mikael Lund/Bjoern Persson
   *  \date Lund/Prag 2008
   *  Lenard-Jones and Debye Huckel poteintial with
   *  square well hydrophobic interaction.
   *  The debth of the well is defined by sticky_u
   *  and the width by sticky_r
   *
   */
  class pot_LJhydrophobicDH : public pot_lj {
    protected:
      double I, k, sticky_r, sticky_u;
    public:
      double box, halfbox;
      pot_LJhydrophobicDH( inputfile &in ) : pot_lj(in) {
        f=in.getflt("bjerrum",7.1);
        k=1/in.getflt("debyelen",1.1e6);
        if ( 1/k>=1e6) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 8*pyc.pi*f*pyc.Nav/1e27*I );
        }
        box=in.getflt("boxlen");
        halfbox=box/2;
        eps=eps/f;
        name+="/Debye-Huckel and hydrophobic w. minimum image";
        sticky_r=in.getflt("sticky_r", 0);   // sticky threshold between *surfaces*
        sticky_u=in.getflt("sticky_u", 0)/f; // sticky energy in units of kT
      }

      /*!
       * \brief   Fast appoximate exponential function
       * \warning On certain linux machines (64 bit but perhaps also 32 bit) this
       *          function returns ZERO when compiled using GCC. The intel compiler
       *          seem OK, though. Solve with -fno-strict-aliasing ???
       * \note    http://firstclassthoughts.co.uk/misc/optimized_power_method_for_java_and_c_and_cpp.html
       */
      double inline f_exp(double x) {
        static double e=2.718281828; //base
        int tmp = (*(1 + (int *) &e));
        int tmp2 = (int)(x * (tmp - 1072632447) + 1072632447);
        double p = 0.0;
        *(1 + (int * ) &p) = tmp2;
        return p;
      }

      /*!
       * \brief Fast appoximate exponential function
       * \url http://www.wilmott.com/messageview.cfm?catid=34&threadid=79761&STARTPAGE=1
       *      http://www.mitpressjournals.org/doi/abs/10.1162/089976699300016467
       * \warning untested!
       *
       * The relative error in the interval [-10;10] is around 0-2%. Seem to perform better
       * than the above f_exp() function.
       */
      inline double f_exp2(double y) {
        double d;
        *((int*)(&d) + 0) = 0;
        *((int*)(&d) + 1) = (int)(1512775 * y + 1072632447);
        return d;
      }

      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in units of kT/lB
      inline double pairpot( const particle &p1, const particle &p2 ) {
        double r2=p1.sqdist(p2,box,halfbox),
               x=p1.radius+p2.radius, u=x*x/r2, r=sqrt(r2);
        x=u*u*u;
        //cout << p1.hydrophobic << endl; 
        if (p1.hydrophobic==true)
          if (p2.hydrophobic==true)
            if (r-p1.radius-p2.radius < sticky_r){
              //cout << "sticky" << endl;                                                 
              return (x*x-x)*eps + p1.charge * p2.charge / r * _EXPFUNC(-k*r) - sticky_u;
            }
          return (x*x-x)*eps + p1.charge * p2.charge / r * _EXPFUNC(-k*r);
        }



      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2,box,halfbox);
      }

      void setvolume(double vol) {
        box=pow(vol, 1./3.);
        halfbox=box/2.;
      }

      int anint(double x) const { return int(x>0. ? x+.5 : x-.5); }
      void boundary(point &a) const {
        if (std::abs(a.x)>halfbox) a.x-=box*anint(a.x/box);
        if (std::abs(a.y)>halfbox) a.y-=box*anint(a.y/box);
        if (std::abs(a.z)>halfbox) a.z-=box*anint(a.z/box);
      }

      virtual double energy(const vector<particle> &p, const vector<int> &pairs) {
        int n=pairs.size()/2, l=0;
        double u=0;
//#pragma omp parallel for reduction (+:u)
        for (int i=0; i<n; i++)
          u+=pairpot(p[ pairs[l++] ], p[ pairs[l++]] );
        return f*u;
      }

      double VectorEnergy( double *r2, double *qq, int *len) {
        int n=*len;
        //std::cout << len << " ";
        double u=0;
        for (int i=1; i<n; i++) {
          double r = sqrt(r2[i]);
          double w = 4.0/r2[i];
          w = w*w*w; //6
          u += qq[i] / r * exp(-k*r) + eps*(w*w-w);
        }
        return f*u;
      }

      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f     << endl
          << "#   Debye length      = " << 1./k  << endl
          << "#   Ionic strength    = " << k*k*1e30/(8*pyc.pi*f*pyc.Nav) << " mM" << endl
          << "#   exp(-lB/D) error  = " << (_EXPFUNC(-f*k)-std::exp(-f*k))/std::exp(-f*k)*100 << " \%" << endl
          << "#   Sticky threshold  = " << sticky_r << " AA\n"                
          << "#   Sticky energy     = " << sticky_u*f << " kT\n";            
        return o.str();
      }
  };

}
#endif
