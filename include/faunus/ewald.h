#ifndef FAUNUS_EWALD_H
#define FAUNUS_EWALD_H

#include "faunus/inputfile.h"
#include <complex>

namespace Faunus {

  namespace Potential {

    /**
     * @brief Pair potential for real-space part of Ewald
     * 
     * @f[
     * E_{Real} = \sum_{i=1}^{N-1}\sum_{j=i+1}^N \left[ \left(\bar{\mu}_i \cdot \bar{\mu}_j\right) B(|\bar{r}_{ij}|) - \left(\bar{\mu}_i \cdot \bar{r}_{ij}\right)\left(\bar{\mu}_j \cdot \bar{r}_{ij}\right)  C(|\bar{r}_{ij}|)  \right]
     * @f]
     * 
     * where
     * 
     * @f[
     * B(r) = \frac{erfc(\alpha r) }{r^3} + \frac{2\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2r^2}}{r^2}
     * @f]
     * 
     * and
     * 
     * @f[
     * C(r) = 3\frac{erfc(\alpha r) }{r^5} + \frac{2\alpha}{\sqrt{\pi}}\left( 2\alpha^2 + \frac{3}{r^2} \right)\frac{e^{-\alpha^2r^2}}{r^2}
     * @f]
     * 
     * @note Optimal parameters for ion-dipole interactions are assumed to be the same as in ion-ion interactions.
     *
     */
    template<bool useIonIon=false, bool useIonDipole=false, bool useDipoleDipole=true> 
      struct EwaldReal : public Potential::Coulomb {

        typedef Potential::Coulomb Tbase;
        double alpha, alpha2, constant, rc, rc2i;

        EwaldReal(Tmjson &j, string sec="ewald") : Tbase(j,sec), alpha(0), rc2i(0) {
          Tbase::name="Ewald Real";
        }

        void updateAlpha(double alpha_in) {
          alpha = alpha_in;
          alpha2 = alpha*alpha;
          constant = 2*alpha/sqrt(pc::pi);
        }

        void updateRcut(double rc_in) {
          rc = rc_in;
          rc2i = 1.0/(rc*rc);
        }

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
            if(!useIonIon && !useIonDipole && !useDipoleDipole)
              return 0.0;

            double r2i = 1.0/r.squaredNorm();
            if (r2i < rc2i)
              return 0.0;

            double E = 0.0;
            double r1i = sqrt(r2i);
            double r1i_d = erfc_x(alpha/r1i)*r1i;

            if(useIonIon)
              E += r1i_d*a.charge*b.charge;
            if(!useIonDipole && !useDipoleDipole)
              return E*Tbase::bjerrumLength();

#ifdef DIPOLEPARTICLE
            double expK = constant*exp(-alpha2/r2i);
            double T1 = (expK + r1i_d)*r2i;

            if(useIonDipole) {
              E += a.charge*r.dot(b.mu)*b.muscalar*T1;
              E -= b.charge*r.dot(a.mu)*a.muscalar*T1;
            }

            if(useDipoleDipole) {
              double t3 = -b.mu.dot(a.mu)*T1;
              double t5 = b.mu.dot(r)*a.mu.dot(r)*(3.0*T1 + 2.0*alpha2*expK )*r2i;
              E += -(t5 + t3)*b.muscalar*a.muscalar;
            }
#endif

            return E*Tbase::bjerrumLength();
          }
      };
  }//namespace

  namespace Energy {

    using namespace Faunus::Potential;

    /**
     * @brief Ewald summation for electrostatic interactions
     * @date Lund 2014
     *
     * This will take care of long-ranged electrostatic interactions using the Ewald summation scheme[@article http://onlinelibrary.wiley.com/doi/10.1002/andp.19213690304/abstract]. Currently it is possible to include ions and dipoles. 
     * 
     *  Keyword          |  Description
     * :--------------   | :---------------
     * `delta`           |  Accuracy in system for energy (\f$ e^2/\AA \f$), force(\f$ e^2/\AA^2 \f$) and torque(\f$ e^2/\AA \f$).           (Default: 5e-5)
     * `eps_surf`        |  Dielectric constant of the surroundings                                                                          (Default: \f$ \varepsilon_r = \infty \f$)
     * `spherical_sum`   |  Spherical (or ellipsoid) summation in reciprocal space if set to be true. Cubic summation is applied if false.   (Default: true)
     * `cutoff`          |  Real space cut-off.                                                                                              (Default: Half the minimum box-length)              
     * `alpha`           |  Damping parameter.                                                                                               (Default: According to DOI: 10.1080/08927029208049126 )  
     * `cutoffK`         |  Maximum number of vectors in any axis in k-space.                                                                (Default: According to DOI: 10.1080/08927029208049126 )
     * `cutoffK_x`       |  Maximum number of vectors in x-axis in k-space for ions. Is overridden if 'cutoffK' or 'cutoffK_ion' is set.     (Default: According to DOI: 10.1080/08927029208049126 )
     * `cutoffK_y`       |  Maximum number of vectors in y-axis in k-space for ions. Is overridden if 'cutoffK' or 'cutoffK_ion' is set.     (Default: According to DOI: 10.1080/08927029208049126 )
     * `cutoffK_z`       |  Maximum number of vectors in z-axis in k-space for ions. Is overridden if 'cutoffK' or 'cutoffK_ion' is set.     (Default: According to DOI: 10.1080/08927029208049126 )
     * 
     * @note Tested and implemented through DOI: 10.1063/1.481216
     * 
     * @f[
     * E_{Total} = E_{Real} + E_{Reciprocal} + E_{Self} + E_{Surface} + E_{Pol}
     * @f]
     * 
     * Where the first of these terms comes from 'EwaldReal' while the others are described by
     * 
     * @f[
     * E_{Reciprocal} = \frac{2\pi}{V}\sum_{ {\bf k} \ne {\bf 0}} A_k\left|Q^{q\mu}\right|^2 \;\;,\;\; A_k = \frac{e^{-k^2/4\alpha^2}}{k^2} \;\;,\;\; Q^{q\mu} = \sum_{j}\left(q_j + i({\boldsymbol \mu}_j\cdot {\bf k})  \right)e^{i({\bf k}\cdot {\bf r}_j)}
     * @f]
     * 
     * @f[
     * E_{Self} = -\sum_{j} \left( \frac{\alpha}{\sqrt{\pi}}q_j^2 + \frac{2\alpha^3}{3\sqrt{\pi}}|{\boldsymbol \mu}_j|^2   \right)
     * @f]
     * 
     * @f[
     * E_{Surface} = \frac{2\pi}{(2\varepsilon_{sur} + 1)V}\left[  \left| \sum_{j}q_j{\boldsymbol r}_j   \right|^2 + 2\left(\sum_{j}q_i{\boldsymbol r}_j \right)\cdot\left(\sum_{j}{\boldsymbol \mu}_j \right) + \left|\sum_{j}{\boldsymbol \mu}_j \right|^2   \right]
     * @f]
     * 
     * @f[
     * E_{Pol} = ...  
     * @f]
     * 
     * where 
     * 
     * @f[
     * {\bf k} = 2\pi\left( \frac{n_x}{L_x} , \frac{n_y}{L_y} ,\frac{n_z}{L_z} \right)  \;\;,\;\; {\bf n} \in \mathbb{Z}^3
     * @f]
     * 
     * @warning Polarization energy is yet to be implemented!
     * 
     * @todo documentation; optimization; automatic alpha calc.
     * 
     */
    template<
      class Tspace, \
      class Tpairpot, \
      bool useIonIon=true, bool useIonDipole=true, bool useDipoleDipole=true, \
      class Tbase=NonbondedVector<Tspace, \
      CombinedPairPotential<EwaldReal<useIonIon,useIonDipole,useDipoleDipole>, Tpairpot>>>

      class NonbondedEwald : public Tbase {
        private:
          using Tbase::spc;
          typedef typename Tbase::Tparticle Tparticle;
          typedef typename Tbase::Tpvec Tpvec;

          Average<double> timeReciprocal;

          int kVectorsLength, kVectorsInUse, N, updates, cnt_accepted, kcc_x, kcc_y, kcc_z, update_frequency;
          double kc_x, kc2_x, kc_y, kc2_y, kc_z, kc2_z, alpha, alpha2, V, eps_surf, eps_r, rc, lB, const_inf, check_k2_zero, delta, minL, maxL,Sq2, Smu2; 
          Point Li, L;                                                                                                     // Inverse box-length-vector, Box-length-vector
          bool user_alpha, user_kc, user_rc, spherical_sum, update_frequency_bool; 
          vector<double> values_alpha, values_rc, values_kcc_x, values_kcc_y, values_kcc_z, values_wavefunctions;
          vector<complex<double>> Q_ion_tot, Q_dip_tot;
	  typename Tspace::Change change;

          Eigen::MatrixXd kVectors;          // Matrices with k-vectors for ions and dipoles respectively
          Eigen::VectorXd k2s, Aks;  // Stores values based on k-vectors in order to minimize computational effort.

          /**
           * @brief Returns Ewald self energy for ions and dipoles.
           * @param p Vector with partciles
           * @param g Group to calculate from
           * @note Does only add to average value if all atoms are counted
           */
          template<class Tpvec, class Tgroup>
            double getSelfEnergy(const Tpvec &p, const Tgroup &g) {  // const
              double Eq = 0;
              double Emu = 0;
              for (auto i : g) {
                if (useIonIon || useIonDipole)
                  Eq += p[i].charge * p[i].charge;
#ifdef DIPOLEPARTICLE
                if (useIonDipole || useDipoleDipole)
                  Emu += p[i].muscalar * p[i].muscalar;
#endif
              }
              Sq2 = Eq;
              Smu2 = Emu;

              double selfEnergy = -alpha*( Sq2 + alpha2*(2.0/3.0)*Smu2 ) / sqrt(pc::pi);
              //selfEnergyAverage += selfEnergy;
              return selfEnergy;
            }

          /**
           * @brief Returns Ewald surface energy for ions and dipoles.
           * @param p Vector with partciles
           * @param g Group to calculate from
           * @warning Have trouble with few (2,...) particles!
           */
          template<class Tpvec, class Tgroup>
            double getSurfaceEnergy(const Tpvec &p, const Tgroup &g) { // const
              Point mus(0,0,0);
              Point qrs(0,0,0);
              for (auto i : g) {  
                if (useIonIon || useIonDipole)
                  qrs = qrs + p[i].charge*p[i];
#ifdef DIPOLEPARTICLE
                if (useIonDipole || useDipoleDipole)
                  mus = mus + p[i].mu*p[i].muscalar;
#endif
              }
              double surfaceEnergy = const_inf * (2*pc::pi/(( 2*eps_surf + 1)*V))*( qrs.dot(qrs) + 2*qrs.dot(mus) +  mus.dot(mus) );
              //surfaceEnergyAverage += surfaceEnergy;
              return surfaceEnergy;
            }

          /**
           * @brief Support-function to 'getReciEnergy'
           * @param p Vector with partciles
           * @param g Group to calculate from
           * @param kv Reciprocal space vector
           */
          template<class Tpvec, class Tgroup>
            double getQ2(const Tpvec &p, const Tgroup &g, const Point &kv, int index) const {
              double Q2 = 0.0;
	      
              complex<double> Q2_ion = Q_ion_tot.at(index);
              complex<double> Q2_dip = Q_dip_tot.at(index);
	      
              if (Tbase::isTrial(p)) {
		
		// Change k-vectors due to moved particles
                for (auto m : change.mvGroup) {
		  auto g = spc->groupList()[m.first];
		  for (auto i : m.second) {
		    double dotTrial = kv.dot(p[i]);
		    double dot = kv.dot(spc->p[i]);
		    if (useIonIon || useIonDipole) {
		      Q2_ion += p[i].charge * complex<double>(cos(dotTrial),sin(dotTrial));
		      Q2_ion -= spc->p[i].charge * complex<double>(cos(dot),sin(dot));
		    }
#ifdef DIPOLEPARTICLE
		    if (useDipoleDipole || useIonDipole) {
		      Q2_dip += kv.dot(p[i].mu) * p[i].muscalar * complex<double>(-sin(dotTrial),cos(dotTrial));
		      Q2_dip -= kv.dot(spc->p[i].mu) * spc->p[i].muscalar * complex<double>(-sin(dot),cos(dot));
		    }
#endif
		  }
                }
                
                // Change k-vectors due to removed particles
                for (auto m : change.rmGroup) {
		  auto g = spc->groupList()[m.first];
		  for (auto i : m.second) {
		    double dotTrial = kv.dot(p[i]);
		    double dot = kv.dot(spc->p[i]);
		    if (useIonIon || useIonDipole)
		      Q2_ion -= spc->p[i].charge * complex<double>(cos(dot),sin(dot));
#ifdef DIPOLEPARTICLE
		    if (useDipoleDipole || useIonDipole)
		      Q2_dip -= kv.dot(spc->p[i].mu) * spc->p[i].muscalar * complex<double>(-sin(dot),cos(dot));
#endif
		  }
                }
                
                // Change k-vectors due to inserted particles
                for (auto m : change.inGroup) {
		  auto g = spc->groupList()[m.first];
		  for (auto &i : m.second) {
		    double dot = kv.dot(i);
		    if (useIonIon || useIonDipole)
		      Q2_ion += i.charge * complex<double>(cos(dot),sin(dot));
#ifdef DIPOLEPARTICLE
		    if (useDipoleDipole || useIonDipole)
		      Q2_dip += kv.dot(i.mu) * i.muscalar * complex<double>(-sin(dot),cos(dot));
#endif
		  }
                }
                
                
                
              }
              
              
              if(useIonIon)
                Q2 += std::norm(Q2_ion);  // 'norm' returns squared magnitude
              if(useIonDipole)
                Q2 += 2.0*std::real(Q2_ion*Q2_dip);
              if(useDipoleDipole)
                Q2 += std::norm(Q2_dip);  // 'norm' returns squared magnitude
              return Q2;
            }

          /**
           * @brief Returns reciprocal space energy.
           * @param p Vector with partciles
           * @param g Group to calculate from
           */
          template<class Tpvec, class Tgroup>
            double getReciprocalEnergy(const Tpvec &p, const Tgroup &g) {
              double E = 0.0;
	      
              for (int k=0; k<kVectorsInUse; k++) {
                double Q2 = getQ2(p,g,kVectors.col(k),k);
                E += ( Aks[k] * Q2 );
              }

              double reciprocalEnergy = (2*pc::pi/V)*E;
              //reciprocalEnergyAverage += reciprocalEnergy;
              return reciprocalEnergy;
            }

          /**
           * @brief Updates all vectors and matrices which depends on the number of k-space vectors.
           * @note Needs to be called whenever 'kcc_x', 'kcc_y' or 'kcc_z' has been updated
           */
          void kVectorChange() {
            kVectorsLength = (2*kcc_x + 1)*(2*kcc_y + 1)*(2*kcc_z + 1) - 1;
            kVectors.resize(3, kVectorsLength); 
            k2s.resize(kVectorsLength);
            Aks.resize(kVectorsLength);
            kVectorsInUse = 0;
            kVectors.setZero();
            k2s.setZero();
            Aks.setZero();

            for (int kx = -kcc_x; kx <= kcc_x; kx++) {
              double dkx2 = double(kx*kx);
              for (int ky = -kcc_y; ky <= kcc_y; ky++) {
                double dky2 = double(ky*ky);
                for (int kz = -kcc_z; kz <= kcc_z; kz++) {
                  double dkz2 = double(kz*kz);
                  Point kv = 2*pc::pi*Point(kx*Li.x(),ky*Li.y(),kz*Li.z());
                  double k2 = kv.dot(kv);
                  if (k2 < check_k2_zero) { // Check if k2 != 0
                    continue;
                  }
                  if(spherical_sum) {
                    if( (dkx2/kc2_x) + (dky2/kc2_y) + (dkz2/kc2_z) > 1.0)
                      continue;
                  }
                  kVectors.col(kVectorsInUse) = kv; 
                  k2s[kVectorsInUse] = k2;
                  Aks[kVectorsInUse] = exp(-k2/(4*alpha2))/k2;
                  kVectorsInUse++;
                }
              }
            }
            values_wavefunctions.push_back(double(kVectorsInUse));
            Q_ion_tot.resize(kVectorsInUse);
            Q_dip_tot.resize(kVectorsInUse);
          }

          /**
           * @brief Solves the equation \f$ w\cdot e^{w} = x \f$
           * @note Implemented through DOI: 10.1145/361952.361970
           */
          double LambertW_solver(double x, double error_max=1e-14) const {
            double wn, wp;
            if(x < 0.7) {
              wp = (x + (4.0/3.0)*x*x ) / ( 1.0 + (7.0/3.0)*x + (5.0/6.0)*x*x );
            } else {
              wp = log(x) - (24.0*(log(x)*log(x) + 2.0*log(x) - 3.0)/(7.0*log(x)*log(x) + 58.0*log(x) + 127.0 ));
            }
            wn = wp + 2.0*error_max;
            double zn, en;
            while( abs(wp - wn) > error_max) {
              wn = wp;
              zn = log(x) - log(wn) - wn;
              en = ( zn / (1.0 + wn) ) * (2.0*(1.0 + wn)*(1.0 + wn + (2.0/3.0)*zn) - zn) / (2.0*(1.0 + wn)*(1.0 + wn + (2.0/3.0)*zn) - 2.0*zn);
              wp = wn*(1.0 + en);
            }
            return wp;
          }

          /**
           * @brief Calculates ionic damping parameter for the Ewald method.
           * @note According to DOI: 10.1080/08927029208049126  (DOI: 10.1007/b136795)
           * @warning 'tau_R' needs to be calculated
           */
          void calcAlphaIon() {
            double tau_R = 1.0; 
            double tau_F = timeReciprocal.avg();
            alpha = pow(tau_R*N/tau_F,1.0/6.0)*sqrt(pc::pi)/minL;
            values_alpha.push_back(alpha);
          }

          /**
           * @brief Calculates ionic real space cutoff parameter for the Ewald method.
           * @param alpha_ion_in Damping parameter (\f$ \AA^{-1} \f$)
           * @param delta_in Accuracy in ionic systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1080/08927029208049126  (DOI: 10.1007/b136795)
           */
          void calcRcIon(double alpha_ion_in, double delta_in = 5e-5) {
            rc = abs(0.5*sqrt(3.0*LambertW_solver( -(4.0/3.0)*pow(-pow(Sq2,4.0)/(alpha_ion_in*alpha_ion_in*pow(minL,6.0)*pow(delta_in,4.0)),1.0/3.0) )));
            values_rc.push_back(rc);
          }

          /**
           * @brief Calculates ionic reciprocal space cutoff parameter for the Ewald method.
           * @param alpha_ion_in Damping parameter (\f$ \AA^{-1} \f$)
           * @param delta_in Accuracy in ionic systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1080/08927029208049126  (DOI: 10.1007/b136795)
           */
          void calcKcIon(double alpha_ion_in, double delta_in = 5e-5) {
            double temp = (1.0/3.0)*pow(2.0,2.0/3.0)*pow(4.0,1.0/3.0)*pow(pow(Sq2,4.0)/(alpha_ion_in*alpha_ion_in*pow(delta_in,4.0)),1.0/3.0);
            kc_x = abs(0.5*sqrt(3.0*LambertW_solver( temp*pow(L.x(), -2.0) )));
            kc_x = abs(0.5*sqrt(3.0*LambertW_solver( temp*pow(L.y(), -2.0) )));
            kc_x = abs(0.5*sqrt(3.0*LambertW_solver( temp*pow(L.z(), -2.0) )));
            kc2_x = kc_x*kc_x;
            kc2_y = kc_y*kc_y;
            kc2_z = kc_z*kc_z;
            kcc_x = ceil(kc_x);
            kcc_y = ceil(kc_y);
            kcc_z = ceil(kc_z);
            values_kcc_x.push_back(kcc_x);
            values_kcc_y.push_back(kcc_y);
            values_kcc_z.push_back(kcc_z);
            kVectorChange();
          }


          /**
           * @brief Calculates dipolar damping parameter for the Ewald method.
           * @note According to DOI: 10.1063/1.1398588
           * @warning 'tau_R' needs to be calculated
           */
          void calcAlphaDipole() {
            double tau_R = 1.0; 
            double tau_F = timeReciprocal.avg();
            alpha = pow(tau_R*N/tau_F,1.0/6.0)*sqrt(pc::pi)/minL;
            values_alpha.push_back(alpha);
          }

          /**
           * @brief Calculates dipolar eciprocal space cut-off for the Ewald method.
           * @param alpha_ion_in Damping parameter
           * @param delta_in Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1063/1.1398588
           */
          void calcRcDipole(double alpha_dip_in, double delta_in = 5e-5) {
            double num = 0.5*delta_in*delta_in*N*pow(minL,3.0)/(Smu2*Smu2*pow(alpha_dip_in,4.0));
            double dennum = -pow(delta_in,4.0)*N*N*pow(minL,6.0) /( pow(alpha_dip_in,6.0)*pow(Smu2,4.0) );
            double denden = LambertW_solver(-3.515625*pow(delta_in,4.0)*N*N*pow(minL,6.0)/(pow(alpha_dip_in,6.0)*pow(Smu2,4.0)));
            rc = num/sqrt(dennum/denden);
            values_rc.push_back(rc);
          }

          /**
           * @brief Calculates dipolar reciprocal space cut-off for the Ewald method.
           * @param alpha_ion_in Damping parameter
           * @param delta Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$), force(\f$ e^2/\AA^2 \f$) and torque(\f$ e^2/\AA \f$).
           * @note According to DOI: 10.1063/1.1398588. Here we assume statistical error in reciprocal space is negligible with regard to the systematic error ( p.6355 )
           */
          void calcKcDipole(double alpha_dip_in, double delta_in = 5e-5) {
            double temp0 = exp(-0.5*LambertW_solver((-1.125*pow(pc::pi*delta/Smu2,2.0)*N*pow(alpha_dip_in,-6.0))));
            kc2_x = 0.75*delta*L.x()*sqrt(N)*temp0/(alpha_dip_in*alpha_dip_in*Smu2);
            kc2_y = 0.75*delta*L.y()*sqrt(N)*temp0/(alpha_dip_in*alpha_dip_in*Smu2);
            kc2_z = 0.75*delta*L.z()*sqrt(N)*temp0/(alpha_dip_in*alpha_dip_in*Smu2);
            kc2_x = kc_x*kc_x;
            kc2_y = kc_y*kc_y;
            kc2_z = kc_z*kc_z;
            kcc_x = ceil(kc_x);
            kcc_y = ceil(kc_y);
            kcc_z = ceil(kc_z);
            values_kcc_x.push_back(kcc_x);
            values_kcc_y.push_back(kcc_y);
            values_kcc_z.push_back(kcc_z);
            kVectorChange();
          }

          /**
           * @brief Calculates optimal parameters for the Ewald method. It will only calculate a value if the user has not defined it in the input.
           * @param delta_in Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$), force(\f$ e^2/\AA^2 \f$) and torque(\f$ e^2/\AA \f$). Should not be larger than \f$ 5\cdot 10^{-5} \f$.
           */
          void calcParameters(double delta_in=5e-5) {
            if(!user_alpha) {
              if(useIonIon || useIonDipole)
                calcAlphaIon();
              if(!useIonIon && !useIonDipole && useDipoleDipole)
                calcAlphaDipole();
            }
            if(!user_rc) {
              if(useIonIon || useIonDipole)
                calcRcIon(alpha,delta_in);
              if(!useIonIon && !useIonDipole && useDipoleDipole)
                calcRcDipole(alpha,delta_in);
            }
            if(!user_kc) {
              if(useIonIon || useIonDipole)
                calcKcIon(alpha,delta_in);
              if(!useIonIon && !useIonDipole && useDipoleDipole)
                calcKcDipole(alpha,delta_in);
            }
            Tbase::pairpot.first.updateAlpha(alpha);
            Tbase::pairpot.first.updateRcut(rc);
          }

          void updateBoxDimensions(Point L_in) {
            L = L_in;
            Li.x() = 1.0/L.x();
            Li.y() = 1.0/L.y();
            Li.z() = 1.0/L.z();
            minL = L.x();
            if(L.y() < minL)
              minL = L.y();
            if(L.z() < minL)
              minL = L.z();
            maxL = L.x();
            if(L.y() > maxL)
              maxL = L.y();
            if(L.z() > maxL)
              maxL = L.z();
            check_k2_zero = 0.1*(4*pc::pi*pc::pi)/(maxL*maxL);
          }

          string _info() {
            using namespace Faunus::textio;
            char w=25;
            std::ostringstream o;
            o << Tbase::_info();
            o << header("Ewald summation");
            o << pad(SUB,w, "Reciprocal cut-off") << "(" << kcc_x << "," << kcc_y << "," << kcc_z << ")" << endl;
            o << pad(SUB,w, "Wavefunctions") << values_wavefunctions.at(values_wavefunctions.size()-1) << endl;
            o << pad(SUB,w, "alpha") << alpha << endl;
            o << pad(SUB,w, "Real cut-off") << rc << endl;
            o << pad(SUB,w+1, epsilon_m+"(Surface)") << eps_surf << endl;
            o << pad(SUB,w, "updates") << updates << endl;

	    /*
            o << pad(SUB,w+4, bracket("Energies in kT")) << endl;
            o << pad(SUB,w+4, bracket("Self energy")) << selfEnergyAverage.mean()*lB << endl;
            o << pad(SUB,w+1, "    "+sigma) << selfEnergyAverage.std()*lB << endl;
            o << pad(SUB,w+4, bracket("Surf energy")) << surfaceEnergyAverage.mean()*lB << endl;
            o << pad(SUB,w+1, "    "+sigma) << surfaceEnergyAverage.std()*lB << endl;
            o << pad(SUB,w+4, bracket("Real energy")) << realEnergyAverage.mean()*lB << endl;
            o << pad(SUB,w+1, "    "+sigma) << realEnergyAverage.std()*lB << endl;
            o << pad(SUB,w+4, bracket("Reci energy")) << reciprocalEnergyAverage.mean()*lB << endl;
            o << pad(SUB,w+1, "    "+sigma) << reciprocalEnergyAverage.std()*lB << endl;
	    */
	    
            return o.str();
          }

        public:
          //MeanValue<double> selfEnergyAverage, surfaceEnergyAverage, realEnergyAverage, reciprocalEnergyAverage;

          NonbondedEwald(Tmjson &j, const string &sec="nonbonded") : Tbase(j,sec)  {
	    // , selfEnergyAverage(10) ,surfaceEnergyAverage(10) , realEnergyAverage(10) , reciprocalEnergyAverage(10)
            Tbase::name += " (Ewald)";
            updateBoxDimensions(typename Tspace::GeometryType(j).len);
	    auto _j = j["energy"]["nonbonded"]["ewald"];
	    
            cnt_accepted = 0;
            updates = 0;
            values_alpha.clear();
            values_rc.clear();
            values_kcc_x.clear();
            values_kcc_y.clear();
            values_kcc_z.clear();
            values_wavefunctions.clear();
            lB = Tbase::pairpot.first.bjerrumLength();
	    eps_surf = ( _j["eps_surf"] | pc::infty );
            const_inf = 1.0;    
            if(eps_surf < 1)  
              const_inf = 0.0; // \varepsilon_{Surface} = \infty
	    
	    spherical_sum = _j["spherical_sum"] | true;
	    delta = ( _j["delta"] | 5e-5 );
            if( ( _j["update_frequency"] | -1 )  > 0 ) {
              update_frequency_bool = true;
              update_frequency = _j["update_frequency"] | -1;
            } else {
              update_frequency_bool = false;
            }

            // Initiate alpha-values
            user_alpha = false;
            if ( ( _j["alpha"] | -1.0 ) > 0.0 ) {
              user_alpha = true;
              alpha = ( _j["alpha"] | -1.0 );
              alpha2 = alpha*alpha;
            } else {
              if(useIonIon || useIonDipole) {
                calcAlphaIon();
              }else {
                calcAlphaDipole();
              }
            }

            // Initiate real cut-off
            user_rc = false;
            if ( ( _j["cutoff"] | -1.0 ) > 0.0 ) {
              user_rc = true;
              rc = ( _j["cutoff"] | -1.0 );
              if (rc > minL/2.0)
                rc = minL/2.0;
            } else {
              if(useIonIon || useIonDipole) {
                calcRcIon(alpha);
              } else {
                calcRcDipole(alpha);
              }
            }

            // Initiate reciprocal cut-off
            user_kc = false;
            if ( ( _j["cutoffK"] | -1.0 ) > 0.99 ) {
              user_kc = true;
              kc_x = ( _j["cutoffK"] | -1.0 );
              kc_y = kc_x;
              kc_z = kc_x;
            } else {
              if( ( _j["cutoffK_x"] | -1.0 ) > 0.99 || ( _j["cutoffK_y"] | -1.0 ) > 0.99 || ( _j["cutoffK_z"] | -1.0 ) > 0.99) {
                user_kc = true;
                kc_x = ( _j["cutoffK_x"] | -1.0 );
                kc_y = ( _j["cutoffK_y"] | kc_x );
                kc_z = ( _j["cutoffK_z"] | kc_y );

                // Assures that kc_x_ion, kc_y_ion and kc_z_ion willhave non-negative values
                if(kc_x < 0.0) {
                  if(kc_y < 0.0) {
                    kc_x = kc_z;
                    kc_y = kc_z;
                  } else {
                    kc_x = kc_y;
                    if(kc_z < 0.0)
                      kc_z = kc_y;
                  }
                } else {
                  if(kc_y < 0.0)
                    kc_y = kc_x;
                  if(kc_z < 0.0)
                    kc_z = kc_x;
                }
              }
            }
            kc2_x = kc_x*kc_x;
            kc2_y = kc_y*kc_y;
            kc2_z = kc_z*kc_z;
            kcc_x = ceil(kc_x);
            kcc_y = ceil(kc_y);
            kcc_z = ceil(kc_z);
	    
            calcParameters(delta);
            if(user_kc)
              kVectorChange();
          }

          double i_external(const Tpvec &p, int i) override {
            Group g = Group(i,i);
            return g_external(p,g);
          }

          double g_external(const Tpvec &p, Group &g) override {
            return lB * getSelfEnergy(p,g);
          }
          
          /**
           * @brief Updates vectors of complex numbers if a move has been accepted. These vectors are implemented in order to increase computational speed. 
           * After 'update_frequency' number of updates the entirety of the complex vectors are recalculated in order to avoid numerical errors.
           * @param move_accepted True if a move is accepted, otherwise false
           */
          virtual double update(bool move_accepted) {
	    if(move_accepted) {
	      cnt_accepted++;
              if(cnt_accepted > update_frequency - 1) {
                updateAllComplexNumbers(spc->trial);
                cnt_accepted = 0;
	      } else {
                for (int k=0; k<kVectorsInUse; k++) {
                  Point kv = kVectors.col(k);
		  
		  // Change k-vectors due to moved particles
                for (auto m : change.mvGroup) {
		  auto g = spc->groupList()[m.first];
		  for (auto i : m.second) {
                    double dotTrial = kv.dot(spc->trial[i]);
                    double dot = kv.dot(spc->p[i]);
                    if(useIonIon || useIonDipole) {
                      Q_ion_tot.at(k) += spc->trial[i].charge * complex<double>(cos(dotTrial),sin(dotTrial));
                      Q_ion_tot.at(k) -= spc->p[i].charge * complex<double>(cos(dot),sin(dot));
                    }
#ifdef DIPOLEPARTICLE 
                    if(useDipoleDipole) {
                      Q_dip_tot.at(k) += kv.dot(spc->trial[i].mu) * spc->trial[i].muscalar * complex<double>(-sin(dotTrial),cos(dotTrial));
                      Q_dip_tot.at(k) -= kv.dot(spc->p[i].mu) * spc->p[i].muscalar * complex<double>(-sin(dot),cos(dot));
                    }
#endif
                  }
		}
                  
                  // Change k-vectors due to removed particles
                for (auto m : change.rmGroup) {
		  auto g = spc->groupList()[m.first];
		  for (auto i : m.second) {
                    double dotTrial = kv.dot(spc->trial[i]);
                    double dot = kv.dot(spc->p[i]);
                    if(useIonIon || useIonDipole)
                      Q_ion_tot.at(k) -= spc->p[i].charge * complex<double>(cos(dot),sin(dot));
#ifdef DIPOLEPARTICLE 
                    if(useDipoleDipole)
                      Q_dip_tot.at(k) -= kv.dot(spc->p[i].mu) * spc->p[i].muscalar * complex<double>(-sin(dot),cos(dot));
#endif
                  }
		}
                  // Change k-vectors due to inserted particles
                for (auto m : change.inGroup) {
		  auto g = spc->groupList()[m.first];
		  for (auto &i : m.second) {
                    double dot = kv.dot(i);
                    if(useIonIon || useIonDipole)
                      Q_ion_tot.at(k) += i.charge * complex<double>(cos(dot),sin(dot));
#ifdef DIPOLEPARTICLE 
                    if(useDipoleDipole)
                      Q_dip_tot.at(k) += kv.dot(i.mu) * i.muscalar * complex<double>(-sin(dot),cos(dot));
#endif
                  }
		}
                }
	      }
	    }
	    return 0.0; 
	  }
	  
          /** @brief Update energy function due to Change */
          virtual double updateChange(const typename Tspace::Change &c) {
	    change = c;
	    return 0.0;
	  }

          double external(const Tpvec &p) override {
            Group g(0, p.size()-1);
            std::chrono::steady_clock::time_point init, final;
            init = std::chrono::steady_clock::now();
            double reciprocalEnergy = getReciprocalEnergy(p,g);
            final = std::chrono::steady_clock::now();
            timeReciprocal += double(std::chrono::duration_cast<std::chrono::milliseconds >(final-init).count());
            return lB * ( getSurfaceEnergy(p,g) + reciprocalEnergy   );
          }

          /**
           * @brief Re-calculates the vectors of complex numbers used in getQ2_ion and getQ2_dip.
           * @param p Particle vector
           */
          void updateAllComplexNumbers(const Tpvec &p) {
              for (int k=0; k<kVectorsInUse; k++) {
                Point kv = kVectors.col(k);
                complex<double> Q_temp_ion(0.0,0.0);
                complex<double> Q_temp_dip(0.0,0.0);
                for (int i = 0; i < N; i++) {
                  double dot = kv.dot(p[i]);
                  if(useIonIon || useIonDipole)
                    Q_temp_ion += p[i].charge * complex<double>(cos(dot),sin(dot));
#ifdef DIPOLEPARTICLE 
                  if(useDipoleDipole)
                    Q_temp_dip += kv.dot(p[i].mu) * p[i].muscalar * complex<double>(-sin(dot),cos(dot));
#endif
                }
                Q_ion_tot.at(k) = Q_temp_ion;
                Q_dip_tot.at(k) = Q_temp_dip;
            }
          }

          /**
           * @brief Set space and update k-vectors and alpha
           */
          void setSpace(Tspace &s) override {
            Tbase::setSpace(s);
            N = s.p.size();
            V = spc->geo.getVolume();
            updateBoxDimensions(spc->geo.len);
            calcParameters();
            updateAllComplexNumbers(s.p);
            if(!update_frequency_bool)
              update_frequency = N;
            updates++;
          }
      };

  }//namespace
}//namespace
#endif

