#ifndef FAUNUS_EWALD_H
#define FAUNUS_EWALD_H

/*!
 * @brief Ewald Summation for long-ranged electrostatics
 *
 * Energies are calculated using Ewald
 * the class also contains a optimization concerning
 * cutoff in the Fourier room
 * Cutoff in realspace is set equal to half the box
 *
 * @author Martin Trulsson
 */

#include "faunus/inputfile.h"
#include <complex>

namespace Faunus {

  template<typename T=double>
    class Ewald {
      private:
        T a1,a2,a3,a4,a5,p1;
        T boxlen, halfboxlen;
        T lB; // Bjerrum length (AA)
        T twopi;
        T ewaldpres;
        T alpha; 
        T alphasqrt;
        T pisqrt;
        T Ewapre;
        T rcut2;   // Squared real space cut-off
        int totk;
        vector<T> kvec;
        vector<complex<T> > eikr, eikrold;
        vector< vector< complex<T> > > eix, eiy, eiz, eixold, eiyold, eizold;
        int kmax, ksqmax;
        void kSpaceInit();//< initialize k-Space

      public:
        Ewald(InputMap&, std::string="ewald_");
        string info();

        void setVolume(T);

        void store() {
          eikrold = eikr;
          eixold  = eix;
          eiyold  = eiy;
          eizold  = eiz;
        }

        void restore() {
          swap(eikr,eikrold);
          swap(eix,eixold);
          swap(eiy,eiyold);
          swap(eiz,eizold);
        }

        /** @brief Real-space particle-particle energy */
        T rSpaceEnergy(T, T);

        /** @brief Update k-vectors */
        template<class Tpvec>
          void kSpaceUpdate(const Tpvec&);

        /** @brief Update k-vectors */
        template<class Tpvec>
          void kSpaceUpdate(const Tpvec&, int);

        /** @brief k-space energy */
        template<class Tpvec>
          T kSpaceEnergy(const Tpvec&);

        /** @brief k-space energy */
        template<class Tpvec>
          T kSpaceEnergy(const Tpvec&, int);

        template<class Tpvec>
          T selfEwald(const Tpvec&);            //< Self-Interaction

        /** @brief Optimize alpha */
        void calcAlpha(int, T=0.001);
    };

  template<class T>
    void Ewald<T>::kSpaceInit() {
      int ksq;
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
              assert(totk<=2000 && "K vector too small");
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
    void Ewald<T>::setVolume(T vol) {
      assert(vol>0);
      boxlen=pow(vol,1/3.);
      halfboxlen=boxlen/2;
      kSpaceInit();
      store();
    }

  /**
   * @param qq particle-particle charge product
   * @param r particle-particle distance
   */
  template<class T>
    T Ewald<T>::rSpaceEnergy(T qq, T r) {
      return  (r*r>rcut2 ? 0 :
          qq*std::erfc(r*alphasqrt)/r );
    }

  template<class T>
    Ewald<T>::Ewald(InputMap &in, std::string pfx) {
      a1=0.254829592;
      a2=-0.284496736;
      a3=1.421413741;
      a4=-1.453152027;
      a5=1.061405429;
      p1=0.3275911;

      int size  = in(pfx+"size",0);
      kmax      = in(pfx+"kmax",5);
      alpha     = in(pfx+"alpha",0.001);
      Ewapre    = in(pfx+"Ewapre",0.01);
      lB        = in("bjerrum",7.12591);
      boxlen    = in("cuboid_len",0);
      rcut2     = in(pfx+"cutoff", boxlen/2);
      rcut2 = pow(rcut2,2);
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
      kSpaceInit();
    }

  template<class T>
    string Ewald<T>::info() {
      std::ostringstream o;
      o << "# kmax                    = " << kmax << endl
        << "# Number of wavefunctions = " << totk << endl;
      return o.str();
    }

  template<class T>
    void Ewald<T>::calcAlpha(int size, T inewaldpression) {
      ewaldpres=inewaldpression;
      int i=0;
      T incr=1e-3, A=10;

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
    template<class Tpvec>
    void Ewald<T>::kSpaceUpdate(const Tpvec &p) {
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
    void Ewald<T>::kSpaceUpdate(const Tpvec &p, int j) {
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
    T Ewald<T>::kSpaceEnergy(const Tpvec &p) {
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
              for (int i=0; i<(int)p.size(); i++)
                sum += p[i].charge * eix[i][kx]*eiy[i][ky]*eiz[i][kz];
              eikr[totk] = sum;
              u += fact*kvec[totk]*std::real( sum*conj(sum) );
              totk+=1;
            }
          }
        }
      }
      return u*lB;
    }

  template<class T>
    template<class Tpvec>
    T Ewald<T>::kSpaceEnergy(const Tpvec &p, int j) {
      complex<T> sum;
      T u=0, fact;
      totk=0;
      for (int kx=0; kx<(kmax+1); kx++) {
        (kx==0 ? fact=1 : fact=2);
        for (int ky=0; ky<(2*kmax+1); ky++) {
          for (int kz=0; kz<(2*kmax+1); kz++) {
            int ksq=kx*kx+(ky-kmax)*(ky-kmax)+(kz-kmax)*(kz-kmax);
            if (ksq < ksqmax && ksq!=0) {
              eikr[totk] = eikrold[totk]+p[j].charge *
                ( eix[j][kx]*eiy[j][ky]*eiz[j][kz]
                  - eixold[j][kx]*eiyold[j][ky]*eizold[j][kz] );
              sum = eikr[totk];
              u += fact*kvec[totk] * std::real( sum*conj(sum) );
              totk+=1;
            }
          }
        }
      }
      return u*lB;
    }

}//namespace
#endif
