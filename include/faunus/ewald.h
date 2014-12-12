#ifndef FAUNUS_EWALD_H
#define FAUNUS_EWALD_H

#include "faunus/inputfile.h"
#include <complex>

namespace Faunus {

  namespace Potential {

    /**
     * @brief Pair potential for real-space part of Ewald
     */
    template<bool useIonIon=true, bool useIonDipole=false, bool useDipoleDipole=false> 
      struct EwaldReal : public Potential::Coulomb {

        typedef Potential::Coulomb base;
        double alpha, rc2i;

        EwaldReal(InputMap &in, string pfx="coulomb") : base(in,pfx), alpha(0), rc2i(0) {
          base::name="Ewald Real";
        }

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
            double r2i = 1.0/r.squaredNorm();
            if (r2i < rc2i)
              return 0;
            double E=0;
            double r1i = sqrt(r2i);
            double r1i_d = erfc_x(alpha/r1i)*r1i;

            if(useIonIon)
              E = r1i_d*a.charge*b.charge;
#ifdef DIPOLEPARTICLE
            double T1;
            if(useIonDipole) {
              expK = constant*exp(-alpha2/r2i);
              T1 = (expK + r1i_d)*r2i;
              E += a.charge*r.dot(b.mu)*b.muscalar*T1;
              E -= b.charge*r.dot(a.mu)*a.muscalar*T1;
            }
            if(useDipoleDipole) {
              if(!useIonDipole) {
                expK = constant*exp(-alpha2/r2i);
                T1 = (expK + r1i_d)*r2i;
              }
              double t3 = -b.mu.dot(a.mu)*T1;
              double t5 = b.mu.dot(r)*a.mu.dot(r)*(3.*r1i_d*r2i + (3.*r2i + 2.*alpha2)*expK)*r2i;
              E += -(t5 + t3)*b.muscalar*a.muscalar;
            }
#endif
            return base::bjerrumLength() * E; 
          }
      };
  }//namespace

  namespace Energy {

    /**
     * @brief Ewald summation for electrostatic interactions
     * @date Lund 2014
     *
     * This will....(description, formulas, refs. etc.)
     *
     *  Keyword          |  Description
     * :--------------   | :---------------
     * `ewald_precision` |  Precision ...
     *
     * @todo arbitrary box dimensions; documentation; optimization; automatic alpha calc.
     *      prettify _info() output (see other energy classes)
     */
    template<class Tspace, class Tpairpot, bool useIonIon=true, bool useIonDipole=false, bool useDipoleDipole=false>
      class NonbondedEwald : public NonbondedVector<Tspace,
      Potential::CombinedPairPotential<Potential::EwaldReal<useIonIon,useIonDipole,useDipoleDipole>, Tpairpot> > {
        private:
          typedef NonbondedVector<Tspace,
                  Potential::CombinedPairPotential<
                    Potential::EwaldReal<useIonIon,useIonDipole,useDipoleDipole>, Tpairpot> > base;

          using base::spc;
          typedef typename base::Tparticle Tparticle;
          typedef typename base::Tpvec Tpvec;

          void updateKVectors() {
            kVectors.setZero();
            k2s.setZero();
            Aks.setZero();
            wavefunctions = 0;
            for (int kx = 0; kx <= maxK; kx++) 
              for (int ky = -maxK; ky <= maxK; ky++) 
                for (int kz = -maxKz; kz <= maxKz; kz++) {
                  Point kv = 2*pc::pi*Point(kx,ky,kz)/L;
                  double k2 = kv.dot(kv);
                  if (k2 < check_k2_zero) // Check if k2 != 0
                    continue;
                  if(k2 > k2max)
                    continue;
                  kVectors.col(wavefunctions) = kv; 
                  k2s[wavefunctions] = k2;
                  Aks[wavefunctions] = exp(-k2/(4*alpha2))/k2;
                  wavefunctions++;
                }
          }

          template<class Tpvec, class Tgroup>
            double getSelfEnergy(const Tpvec &p, const Tgroup &g) {
              double Eq = 0.0;
              double Emu = 0.0;
              for (auto n : g) {
                if (useIonIon || useIonDipole)
                  Eq = Eq + p[n].charge*p[n].charge;
#ifdef DIPOLEPARTICLE
                if (useIonDipole || useDipoleDipole)
                  Emu = Emu + p[n].mu.dot(p[n].mu)*p[n].muscalar*p[n].muscalar;
#endif
              }
              return -(alpha/sqrt(pc::pi))*( Eq + (2*alpha2/3)*Emu );
            }

          template<class Tpvec, class Tgroup>
            double getSurfEnergy(const Tpvec &p, const Tgroup &g) {
              Point mus(0,0,0), qrs(0,0,0);
              for (auto n : g) {
                if (useIonIon || useIonDipole)
                  qrs = qrs + p[n].charge*p[n];
#ifdef DIPOLEPARTICLE
                if (useIonDipole || useDipoleDipole)
                  mus = mus + p[n].mu*p[n].muscalar;
#endif
              }
              return const_inf*(2*pc::pi/(( 2*eps_surf + 1)*V))*(qrs.dot(qrs) + 2*qrs.dot(mus) +  mus.dot(mus) );
            }

          /**
           * @brief ... 
           */
          template<class Tpvec, class Tgroup>
            double getQ2(const Tpvec &p, const Tgroup &g, const Point &kv) {
              std::complex<double> Q2(0, 0);
              for (auto i : g) {
                double dotP = kv.dot( p[i] );
                if (useIonIon || useIonDipole) {
                  std::complex<double> temp0(cos(dotP),sin(dotP));
                  Q2 = Q2 + p[i].charge*temp0;
                }
#ifdef DIPOLEPARTICLE
                if (useIonDipole || useDipoleDipole) {
                  std::complex<double> temp1(-sin(dotP),cos(dotP));
                  Q2 = Q2 + kv.dot(p[i].mu)*p[i].muscalar*temp1;
                }
#endif
              }
              double dQ2 = std::abs(Q2);
              return dQ2*dQ2;
            }

          template<class Tpvec, class Tgroup>
            double getReciEnergy(const Tpvec &p, const Tgroup &g) {
              double E = 0.0;
              if(alpha2 == 0)
                return E;
              double Q2;
              for(int i = 0; i < kVectorsLength; i++) {
                Q2 = getQ2(p,g,kVectors.col(i));
                E += Aks[i]*Q2;
              }
              return (2*pc::pi/V)*E;
            }


          void calcMaxK(double alpha2_in, double L_in) {
            maxK  = std::ceil(1.0 + alpha2_in*Rc*L_in/pc::pi);
            maxKz = std::ceil(1.0 + alpha2_in*Rc*lamsep/pc::pi);
            kVectorsLength = (maxK + 1)*(2*maxK + 1)*(2*maxKz + 1) - 1;
            kVectors.resize(3, kVectorsLength); 
            k2s.resize(kVectorsLength);
            Aks.resize(kVectorsLength);
          }

          /**
           * @brief ....
           * @param rcut2 Real space cut-off
           * @param ....
           */
          void calcAlpha(double inewaldpression, double rcut2) {
            double ewaldpres = inewaldpression;
            int i = 0;
            double incr = 1e-5;
            double A = 10.0;

            while( A > (1.+1e-10) || A < (1-1e-10)) {
              A=sqrt(alpha) * ewaldpres / exp(-alpha * rcut2);

              if(A>1){
                while(A > 1.){
                  alpha-= incr;
                  A = sqrt(alpha) * ewaldpres / exp(-alpha*rcut2);
                  i++;
                  if( i > 10000){
                    cout << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
                    exit(1);
                    return;
                  };
                };
                if( i > 10000){
                  cout << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
                  exit(1);
                  return;
                };

                if( A < (1.-1e-10)) {
                  alpha += 2*incr;
                  A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
                  incr = incr*0.1;
                };  
              } else {
                while( A < 1.){
                  alpha += incr;
                  A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
                  if( i > 10000){
                    cout << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
                    exit(1);
                    return;
                  };
                };
                if( i > 10000){
                  cout << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
                  exit(1);
                  return;
                };

                if( A < (1.-1e-10)) {
                  alpha += 2*incr;
                  A = sqrt(alpha)*ewaldpres / exp(-alpha*rcut2);
                  incr = incr*0.1;
                };

              }; 
            };

            assert( alpha >= 0 && "Alpha negative!");

            alpha2 = alpha;
            alpha = sqrt(alpha);
            constant = 2*alpha/sqrt(pc::pi);
          }

          string _info() {
            std::ostringstream o;
            o << base::_info()
              << "maxK                    = " << maxK << endl
              << "wavefunctions           = " << wavefunctions << endl
              << "alpha                   = " << alpha << endl
              << "R_c                     = " << Rc << endl
              << "<Self energy>           = " << selfEA.avg() << endl
              << "<Surf energy>           = " << surfEA.avg() << endl
              << "<Real energy>           = " << realEA.avg() << endl
              << "<Reci energy>           = " << reciEA.avg() << endl;
            return o.str();
          }

          int N, maxK, maxKz, kVectorsLength, k2max, wavefunctions;
          double alpha, alpha2, V, eps_surf, eps_r, Rc, Rc2, rc2i,
                 constant, L, lB, T, const_inf,
                 check_k2_zero, ewapre, lamsep, lamell;

          bool user_alpha, user_maxK;

          Average<double> selfEA, surfEA, realEA, reciEA;
          Eigen::MatrixXd kVectors;
          Eigen::VectorXd k2s, Aks;  // explain...

        public:
          NonbondedEwald(InputMap &in, string pfx="ewald") : base(in) {
            base::name += " (Ewald)";
            ewapre = in.get<double>(pfx+"_precision",0.01);

            lB = base::pairpot.first.bjerrumLength();
            L = typename Tspace::GeometryType(in).len.norm(); //in.get<double>("cuboid_len",pc::infty);
            lamell = in.get<double>(pfx+"_lamell",0.0);
            lamsep = lamell + L;
            Rc = in.get<double>(pfx+"_cutoff", L/2);
            if (Rc>L/2)
              Rc=L/2.;
            Rc2 = Rc*Rc;
            rc2i = 1.0/(Rc*Rc);
            base::pairpot.first.rc2i=rc2i;
            eps_surf = in.get<double>(pfx+"_eps_surf",pc::infty);
            const_inf = 1;
            if(eps_surf > 1e10)  // \epsilon_{Surface} = \infty
              const_inf = 0.0;

            alpha = 0.01;
            user_alpha = false;
            if (in.get<double>(pfx+"_alpha", -1.0  ) != -1.0 ) {
              user_alpha = true;
              alpha = in.get<double>(pfx+"_alpha", -1.0  );
              constant = 2*alpha/sqrt(pc::pi);
            }

            user_maxK = false;
            if (in.get<double>(pfx+"_maxK", -1.0  ) != -1.0 ) {
              user_maxK = true;
              maxK = in.get<double>(pfx+"_maxK", -1.0  );
              maxKz = maxK;
              kVectorsLength = (maxK + 1)*(2*maxK + 1)*(2*maxKz + 1) - 1;
              kVectors.resize(3, kVectorsLength); 
              k2s.resize(kVectorsLength);
              Aks.resize(kVectorsLength);
            }
          }

          double g_external(const Tpvec &p, Group &g) FOVERRIDE {
            return lB * ( getSelfEnergy(p,g) + getReciEnergy(p,g) );
          }

          double external(const Tpvec &p) FOVERRIDE {
            Group g(0, p.size());
            return lB * getSurfEnergy(p,g);
          }

          /**
           * @brief Set space and update k-vectors and alpha
           */
          void setSpace(Tspace &s) {
            base::setSpace(s);
            V = spc->geo.getVolume();
            L = std::cbrt(V);
            N = 0;
            for (auto &i : spc->p)
              if (fabs(i.charge)>1e-6)
                N++;
            if(!user_alpha)
              calcAlpha(ewapre,Rc2);
            alpha2 = alpha*alpha;
            if(!user_maxK)
              calcMaxK(alpha2,L);
            check_k2_zero = 0.1*L*L/(4*pc::pi*pc::pi);

            k2max = std::ceil(4.0*pow(((pc::pi/L) + alpha2*Rc),2));  // From Martin's program
            int k2maxtemp = std::ceil(4*pow(((pc::pi/lamsep) + alpha2*Rc),2));
            if(k2maxtemp > k2max) 
              k2max = std::ceil(k2maxtemp);

            updateKVectors();

            base::pairpot.first.alpha = alpha;
          }
      };

  }//namespace
}//namespace
#endif
