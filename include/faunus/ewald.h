#ifndef FAUNUS_EWALD_H
#define FAUNUS_EWALD_H

/*!
 * @brief Ewald Summation for long-ranged electrostatics
 * @author Martin Trulsson
 *
 * Energies are calculated using Ewald
 * the class also contains a optimization concerning
 * cutoff in the Fourier room
 * Cutoff in realspace is set equal to half the box
 */

#include "faunus/inputfile.h"
#include <complex>

namespace Faunus {

  template<typename T=double>
    class Ewald {
      private:
        const T a1=0.254829592, a2=-0.284496736, a3=1.421413741,
              a4=-1.453152027, a5=1.061405429, p1=0.3275911;
        T boxlen, halfboxlen;
        T lB; // Bjerrum length (AA)
        T twopi;
        T ewaldpres;
        T alpha; 
        T alphasqrt;
        T pisqrt;
        T Ewapre;
        int totk;
        vector<T> kvec;
        vector<complex<T> > eikr, eikrold;
        vector< vector< complex<T> > > eix, eiy, eiz, eixold, eiyold, eizold;
        int kmax, ksqmax;
        inline T erfc(T);                                 //< Error function (comblimentary)
        void initKSpaceEwald();                                     //< initialize k-Space

      public:
        Ewald(InputMap&, std::string="ewald_");
        string info();

        void setvolume(T);

        template<class Tparticle>
          T pairpot(const Tparticle&, const Tparticle&, T);  //< particle<->particle (real space)

        template<class Tpvec>
          void kSpaceEwald(const Tpvec&);       //< k-Space Ewald all

        template<class Tpvec>
          void kSpaceEwald(const Tpvec&, int);  //< k-Space Ewald particle i

        template<class Tpvec>
          T sumkSpaceEwald(const Tpvec&);

        template<class Tpvec>
          T sumkSpaceEwald(const Tpvec&, int);

        template<class Tpvec>
          T selfEwald(const Tpvec&);            //< Self-Interaction

        void calcAlphaEwald(int, T=0.001);      //< Optimize alpha
    };

  template<class T>
    void Ewald<T>::setvolume(T vol) {
      boxlen=pow(vol,1/3.);
      halfboxlen=boxlen/2;
      initKSpaceEwald();
    }

  template<class T>
    template<class Tparticle>
    T Ewald<T>::pairpot(const Tparticle &p1, const Tparticle &p2, T r2) {
      if (r2 > halfboxlen*halfboxlen)
        return 0.;
      r2 = sqrt(r2);
      return p1.charge*p2.charge*erfc(r2*alphasqrt)/r2;
    }

  /**
   * @note Reference for this approximation is found in Abramowitz and Stegun,
   *       Handbook of mathematical functions, Eq. 7.1.26
   */
  template<class T>
    T Ewald<T>::erfc(T x) {
      T t = 1./(1.+p1*x),
        tp = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
      return tp*exp(-x*x);
    }

  template<class T>
    Ewald<T>::Ewald(InputMap &in, std::string pfx) {
      int size  = in(pfx+"size",0);
      kmax      = in(pfx+"kmax",5);
      alpha     = in(pfx+"alpha",0.001);
      Ewapre    = in(pfx+"Ewapre",0.01);
      lB        = in("bjerrum",7.12591);
      boxlen    = in("cuboid_len",0);
      halfboxlen= boxlen/2;
      alphasqrt=sqrt(alpha);
      T pi=std::acos(-1.);
      twopi = 2*pi;
      pisqrt = sqrt(pi);
      ksqmax=kmax*kmax+1;
      eix.resize(size);
      eiy.resize(size);
      eiz.resize(size);
      eixold.resize(size);
      eiyold.resize(size);
      eizold.resize(size);
      for (int i=0; i<size; i++) {
        eix[i].resize(kmax+1);
        eiy[i].resize(2*kmax+1);
        eiz[i].resize(2*kmax+1);
        eixold[i].resize(kmax+1);
        eiyold[i].resize(2*kmax+1);
        eizold[i].resize(2*kmax+1);
      }
      initKSpaceEwald();
    }

  template<class T>
    string Ewald<T>::info() {
      std::ostringstream o;
      o << "# kmax                    = " << kmax << endl
        << "# Number of wavefunctions = " << totk << endl;
      return o.str();
    }

  template<class T>
    void Ewald<T>::calcAlphaEwald(int size, T inewaldpression) {
      ewaldpres=inewaldpression;
      int i=0;
      T incr=1e-3, A=10, rcut2=halfboxlen*halfboxlen;

      while( A > (1.+1e-10) || A < (1-1e-10)) {
        A=sqrt(alpha) * ewaldpres / exp( -alpha*rcut2 );

        if (A>1) {
          while (A > 1.) {
            alpha-= incr;
            A = sqrt(alpha) * ewaldpres / exp(-alpha*rcut2);
            i++;
            assert(i<=10000);
          }
          assert(i<=10000);

          if (A < (1.-1e-10)) {
            alpha += 2*incr;
            A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
            incr = incr*0.1;
          }  
        } else {
          while (A < 1.) {
            alpha += incr;
            A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
            assert(i<=10000);
          }
          assert(i<=10000);

          if (A < (1.-1e-10)) {
            alpha += 2*incr;
            A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
            incr = incr*0.1;
          }
        }
      }

      assert(alpha>=0 && "Alpha negative!");

      cout << "# New alpha fit = " << alpha << endl;
      alphasqrt=sqrt(alpha);

      kmax = 1. + alpha*halfboxlen*boxlen/twopi*2.;
      //ksqmax = (1.+alpha*halfboxlen*boxlen/twopi*2.)*(1.+alpha*halfboxlen*boxlen/twopi*2.);
      ksqmax = kmax*kmax;
      cout << "# New kmax fit = " << kmax << endl;

      for (int i=0; i<size; i++) {
        eix[i].resize(kmax+1);
        eiy[i].resize(2*kmax+1);
        eiz[i].resize(2*kmax+1);
        eixold[i].resize(kmax+1);
        eiyold[i].resize(2*kmax+1);
        eizold[i].resize(2*kmax+1);
      }
    }

  template<class T>
    template<class Tpvec>
    T Ewald<T>::selfEwald(const Tpvec &p) {
      T u=0;
      for (auto &i : p)
        u += i.charge * i.charge;
      return -u*lB*alphasqrt/pisqrt;
    }

  template<class T>
    void Ewald<T>::initKSpaceEwald() {
      int maxk=2000, ksq;
      T b;
      T rkx,rky,rkz;
      T rksq;
      b=1./4./alpha;
      T twopii=twopi/boxlen;
      T vol=pow(boxlen,3);
      totk=0;

      for (int kx=0; kx<(kmax+1); kx++) {      // uses symmetry
        rkx = twopii*T(kx);
        for (int ky=-kmax; ky<(kmax+1); ky++) {
          rky = twopii*T(ky);
          for (int kz=-kmax; kz<(kmax+1); kz++) {
            rkz = twopii*T(kz);
            ksq = kx*kx+ky*ky+kz*kz;
            if (ksq < ksqmax && ksq!=0) {
              totk+=1;
              assert(totk<=maxk && "K vector too small");
              rksq = rkx*rkx + rky*rky + rkz*rkz;
              kvec.push_back( twopi * exp( -b*rksq ) / rksq / vol );
            }
          }
        }
      }
      kvec.resize(totk);
      eikr.resize(totk);
      eikrold.resize(totk);
      cout << "# Number of wavefunctions :" << totk << endl;
    }

  template<class T>
    template<class Tpvec>
    void Ewald<T>::kSpaceEwald(const Tpvec &p) {
      int size=p.size();
      T twopii=twopi/boxlen;

      for (int i=0; i<size; i++) {
        eix[i][0]=complex<T>(1,0);
        eiy[i][kmax]=complex<T>(1,0);      //kmax = defines the '0'
        eiz[i][kmax]=complex<T>(1,0);

        eix[i][1] = complex<T>( cos( twopii*p[i].x() ), sin( twopii*p[i].x() ) );
        eiy[i][kmax+1] = complex<T>( cos( twopii*p[i].y() ), sin( twopii*p[i].y() ) );
        eiz[i][kmax+1] = complex<T>( cos( twopii*p[i].z() ), sin( twopii*p[i].z() ) );

        eiy[i][kmax-1] = conj(eiy[i][kmax+1]);
        eiz[i][kmax-1] = conj(eiz[i][kmax+1]);
      }

      for (int kk=2; kk<(kmax+1); kk++) {
        for (int i=0; i<size; i++) {
          eix[i][kk]=eix[i][kk-1]*eix[i][1];
          eiy[i][kmax+kk]=eiy[i][kmax+kk-1]*eiy[i][kmax+1];
          eiy[i][kmax-kk]=conj(eiy[i][kmax+kk]);
          eiz[i][kmax+kk]=eiz[i][kmax+kk-1]*eiz[i][kmax+1];
          eiz[i][kmax-kk]=conj(eiz[i][kmax+kk]);
        }
      }
    }

  template<class T>
    template<class Tpvec>
    void Ewald<T>::kSpaceEwald(const Tpvec &p, int j) {
      T twopii=twopi/boxlen;

      eix[j][0]    = complex<T>(1,0);
      eiy[j][kmax] = complex<T>(1,0);
      eiz[j][kmax] = complex<T>(1,0);

      eix[j][1]      = complex<T>( cos(twopii*p[j].x() ), sin(twopii*p[j].x()) );
      eiy[j][kmax+1] = complex<T>( cos(twopii*p[j].y() ), sin(twopii*p[j].y()) );
      eiz[j][kmax+1] = complex<T>( cos(twopii*p[j].z() ), sin(twopii*p[j].z()) );

      eiy[j][kmax-1] = conj( eiy[j][kmax+1] ); //1+kmax = -1   2+kmax = -2 osv.
      eiz[j][kmax-1] = conj( eiz[j][kmax+1] );

      for (int kk=2; kk<(kmax+1); kk++) {
        eix[j][kk]      = eix[j][kk-1] * eix[j][1];
        eiy[j][kmax+kk] = eiy[j][kmax+kk-1] * eiy[j][kmax+1];
        eiy[j][kmax-kk] = conj(eiy[j][kmax+kk]);
        eiz[j][kmax+kk] = eiz[j][kmax+kk-1] * eiz[j][kmax+1];
        eiz[j][kmax-kk] = conj( eiz[j][kmax+kk] );
      }
    }

  template<class T>
    template<class Tpvec>
    T Ewald<T>::sumkSpaceEwald(const Tpvec &p) {
      int size=p.size();
      complex<T> sum;
      int ksq;
      T u=0, fact;
      totk=0;
      for (int kx=0; kx<(kmax+1); kx++) {
        if (kx==0)
          fact=1.;
        else
          fact=2.;
        for (int ky=0; ky<(2*kmax+1); ky++) {
          for (int kz=0; kz<(2*kmax+1); kz++) {
            ksq=kx*kx+(ky-kmax)*(ky-kmax)+(kz-kmax)*(kz-kmax);
            if ( ksq < ksqmax && ksq!=0 ) {
              sum = complex<T>(0,0);
              for (int i=0; i<size; i++)
                sum += p[i].charge * eix[i][kx]*eiy[i][ky]*eiz[i][kz];
              eikr[totk] = sum;
              u += fact*kvec[totk]*real(sum*conj(sum));
              totk+=1;
            }
          }
        }
      }
      return u*lB;
    }

  template<class T>
    template<class Tpvec>
    T Ewald<T>::sumkSpaceEwald(const Tpvec &p, int j) {
      complex<T> sum;
      int ksq;
      T u=0, fact;
      totk=0;
      for (int kx=0; kx<(kmax+1); kx++) {
        if (kx==0)
          fact=1.0;
        else
          fact=2.0;
        for (int ky=0; ky<(2*kmax+1); ky++) {
          for (int kz=0; kz<(2*kmax+1); kz++) {
            ksq=kx*kx+(ky-kmax)*(ky-kmax)+(kz-kmax)*(kz-kmax);
            if (ksq < ksqmax && ksq!=0) {
              eikr[totk] = eikrold[totk]+p[j].charge*(eix[j][kx]*eiy[j][ky]*eiz[j][kz]
                  -eixold[j][kx]*eiyold[j][ky]*eizold[j][kz]);
              sum = eikr[totk];
              u += fact*kvec[totk]*real(sum*conj(sum));
              totk+=1;
            }
          }
        }
      }
      return u*lB;
    }

}//namespace
#endif
