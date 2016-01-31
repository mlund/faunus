#ifndef FAUNUS_EWALD_H
#define FAUNUS_EWALD_H

#include <faunus/auxiliary.h>
#include "faunus/inputfile.h"
#include <complex>

namespace Faunus {

  namespace Potential {

    /**
     * @brief Pair potential for real-space part of Ewald.
     * 
     * @f[
     * \boldsymbol{{\rm T}}_0(\boldsymbol{r}) = \frac{{\rm erfc}(\alpha |\boldsymbol{r}|)}{|\boldsymbol{r}|}
     * @f]
     * 
     * @f[
     * \boldsymbol{{\rm T}}_1(\boldsymbol{r}) = \nabla\left(\frac{{\rm erfc}(\alpha |\boldsymbol{r}|)}{|\boldsymbol{r}|}\right)
     * @f]
     * 
     * @f[
     * \boldsymbol{{\rm T}}_2(\boldsymbol{r}) = \nabla^2\left(\frac{{\rm erfc}(\alpha |\boldsymbol{r}|)}{|\boldsymbol{r}|}\right)
     * @f]
     * 
     * @f[
     * E_{Real} = \sum_{i=1}^{N-1}\sum_{j=i+1}^N \left( q_i\boldsymbol{{\rm T}}_0(\boldsymbol{r}_{ij})q_j   +      q_i\boldsymbol{{\rm T}}_1(\boldsymbol{r}_{ij})\cdot \boldsymbol{\mu}_j  - q_j\boldsymbol{{\rm T}}_1(\boldsymbol{r}_{ij})\cdot \boldsymbol{\mu}_i   -       \boldsymbol{\mu}_i^T\boldsymbol{{\rm T}}_2(\boldsymbol{r}_{ij})\boldsymbol{\mu}_j    \right)
     * @f]
     * 
     * @note Optimal parameters for ion-dipole interactions are assumed to be the same as for ion-ion interactions.
     *
     */
    template<bool useIonIon=true, bool useIonDipole=false, bool useDipoleDipole=false> 
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
      
    template<bool useIonIon=true, bool useIonDipole=false, bool useDipoleDipole=false>
      struct EwaldParameters {

	bool user_alpha, user_rc, user_kc;
	int N;
	double tau_R, tau_F, minL, maxL, Sq2, Smu2, delta, check_k2_zero;
	Point L, Li;
	double alpha, alpha2, rc, kc_x, kc_y, kc_z, kc2_x, kc2_y, kc2_z, kcc_x, kcc_y, kcc_z;
	vector<double> values_alpha, values_rc, values_kcc_x, values_kcc_y, values_kcc_z;
	
        EwaldParameters() { }

          /**
           * @brief Calculates ionic damping parameter for the Ewald method.
	   * @param tau_R_in Time needed in real space to execute one cycle
	   * @param tau_F_in Time needed in reciprocal space to execute one cycle
	   * @param N_in Number of particles
	   * @param minL_in Shortest side of the box
           * @note According to DOI: 10.1080/08927029208049126  (DOI: 10.1007/b136795)
           * @warning 'tau_R' needs to be calculated
           */
          double calcAlphaIon(double tau_R_in, double tau_F_in, int N_in, double minL_in) {
            alpha = pow(tau_R_in*double(N_in)/tau_F_in,1.0/6.0)*sqrt(pc::pi)/minL_in;
            values_alpha.push_back(alpha);
	    alpha2 = alpha*alpha;
	    return alpha;
          }
          
          double calcAlphaIon() {
	    return calcAlphaIon(tau_R,tau_F,N,minL);
	  }

          /**
           * @brief Calculates ionic real space cutoff parameter for the Ewald method.
           * @param alpha_ion_in Damping parameter (\f$ \AA^{-1} \f$)
	   * @param Sq2_in Sum of the squared charges 
	   * @param minL_in Shortest side of the box
           * @param delta_in Accuracy in ionic systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1080/08927029208049126  (DOI: 10.1007/b136795)
           */
          double calcRcIon(double alpha_ion_in, double Sq2_in, double minL_in, double delta_in = 5e-5) {
            rc = abs(0.5*sqrt(3.0*LambertW( -(4.0/3.0)*pow(-pow(Sq2_in,4.0)/(alpha_ion_in*alpha_ion_in*pow(minL_in,6.0)*pow(delta_in,4.0)),1.0/3.0) )));
            values_rc.push_back(rc);
	    return rc;
          }
          
          double calcRcIon() {
            return calcRcIon(alpha,Sq2,minL,delta);
          }

          /**
           * @brief Calculates ionic reciprocal space cutoff parameter for the Ewald method.
           * @param alpha_ion_in Damping parameter (\f$ \AA^{-1} \f$)
	   * @param Sq2_in Sum of the squared charges 
	   * @param L_in Vector with box-dimansions
           * @param delta_in Accuracy in ionic systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1080/08927029208049126  (DOI: 10.1007/b136795)
           * @warning Need to call 'kVectorChange' after change of number of k-vectors
           */
          Point calcKcIon(double alpha_ion_in, double Sq2_in, Point L_in, double delta_in = 5e-5) {
            double temp = (1.0/3.0)*pow(2.0,2.0/3.0)*pow(4.0,1.0/3.0)*pow(pow(Sq2_in,4.0)/(alpha_ion_in*alpha_ion_in*pow(delta_in,4.0)),1.0/3.0);
            kc_x = abs(0.5*sqrt(3.0*LambertW( temp*pow(L_in.x(), -2.0) )));
            kc_y = abs(0.5*sqrt(3.0*LambertW( temp*pow(L_in.y(), -2.0) )));
            kc_z = abs(0.5*sqrt(3.0*LambertW( temp*pow(L_in.z(), -2.0) )));
            kc2_x = kc_x*kc_x;
            kc2_y = kc_y*kc_y;
            kc2_z = kc_z*kc_z;
            kcc_x = ceil(kc_x);
            kcc_y = ceil(kc_y);
            kcc_z = ceil(kc_z);
            values_kcc_x.push_back(kcc_x);
            values_kcc_y.push_back(kcc_y);
            values_kcc_z.push_back(kcc_z);
	    return Point(kcc_x,kcc_y,kcc_z);
          }
          
          Point calcKcIon() {
            return calcKcIon(alpha,Sq2,L,delta);
          }
        
          /**
           * @brief Calculates dipolar damping parameter for the Ewald method.
	   * @param tau_R_in Time needed in real space to execute one cycle
	   * @param tau_F_in Time needed in reciprocal space to execute one cycle
	   * @param N_in Number of particles
	   * @param minL_in Shortest side of the box
           * @note According to DOI: 10.1063/1.1398588
           * @warning 'tau_R' needs to be calculated
           */
          double calcAlphaDipole(double tau_R_in, double tau_F_in, int N_in, double minL_in) {
            alpha = pow(tau_R_in*double(N_in)/tau_F_in,1.0/6.0)*sqrt(pc::pi)/minL_in;
            values_alpha.push_back(alpha);
	    alpha2 = alpha*alpha;
	    return alpha;
          }
          
          double calcAlphaDipole() {
	    return calcAlphaDipole(tau_R,tau_F,N,minL);
	  }

          /**
           * @brief Calculates dipolar reciprocal space cut-off for the Ewald method.
           * @param alpha_ion_in Damping parameter
	   * @param Sq2_in Sum of the squared abolute dipole moments 
	   * @param N_in Number of particles
	   * @param minL_in Shortest side of the box
           * @param delta_in Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1063/1.1398588
           */
          double calcRcDipole(double alpha_dip_in, double Smu2_in, int N_in, double minL_in, double delta_in = 5e-5) {
            double num = 0.5*delta_in*delta_in*double(N_in)*pow(minL_in,3.0)/(Smu2_in*Smu2_in*pow(alpha_dip_in,4.0));
            double dennum = -pow(delta_in,4.0)*double(N_in*N_in)*pow(minL_in,6.0) /( pow(alpha_dip_in,6.0)*pow(Smu2_in,4.0) );
            double denden = LambertW(-3.515625*pow(delta_in,4.0)*double(N_in*N_in)*pow(minL_in,6.0)/(pow(alpha_dip_in,6.0)*pow(Smu2_in,4.0)));
            rc = num/sqrt(dennum/denden);
            values_rc.push_back(rc);
	    return rc;
          }
          
          double calcRcDipole() {
            return calcRcDipole(alpha,Smu2,minL,delta);
          }

          /**
           * @brief Calculates dipolar reciprocal space cut-off for the Ewald method.
           * @param alpha_ion_in Damping parameter
	   * @param Sq2_in Sum of the squared abolute dipole moments 
	   * @param L_in Vector with box-dimansions
	   * @param N_in Number of particles
           * @param delta Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$), force(\f$ e^2/\AA^2 \f$) and torque(\f$ e^2/\AA \f$).
           * @note According to DOI: 10.1063/1.1398588. Here we assume statistical error in reciprocal space is negligible with regard to the systematic error ( p.6355 )
           * @warning Need to call 'kVectorChange' after change of number of k-vectors
           */
          Point calcKcDipole(double alpha_dip_in, double Smu2_in, Point L_in, int N_in, double delta_in = 5e-5) {
            double temp0 = exp(-0.5*LambertW((-1.125*pow(pc::pi*delta_in/Smu2_in,2.0)*double(N_in)*pow(alpha_dip_in,-6.0))));
            kc2_x = 0.75*delta_in*L_in.x()*sqrt(double(N_in))*temp0/(alpha_dip_in*alpha_dip_in*Smu2_in);
            kc2_y = 0.75*delta_in*L_in.y()*sqrt(double(N_in))*temp0/(alpha_dip_in*alpha_dip_in*Smu2_in);
            kc2_z = 0.75*delta_in*L_in.z()*sqrt(double(N_in))*temp0/(alpha_dip_in*alpha_dip_in*Smu2_in);
            kc2_x = kc_x*kc_x;
            kc2_y = kc_y*kc_y;
            kc2_z = kc_z*kc_z;
            kcc_x = ceil(kc_x);
            kcc_y = ceil(kc_y);
            kcc_z = ceil(kc_z);
            values_kcc_x.push_back(kcc_x);
            values_kcc_y.push_back(kcc_y);
            values_kcc_z.push_back(kcc_z);
	    return Point(kcc_x,kcc_y,kcc_z);
          }
          
          Point calcKcDipole() {
            return calcKcDipole(alpha,Smu2,L,N,delta);
          }
          
          void update(double tau_R_in, double tau_F_in, int N_in, double Sq2_in, double Smu2_in, Point L_in) {
	    tau_R = tau_R_in;
	    tau_F = tau_F_in;
	    N = N_in;
	    Sq2 = Sq2_in;
	    Smu2 = Smu2_in;
	    L = L_in;
	    Li.x() = 1.0/L.x();
	    Li.y() = 1.0/L.y();
	    Li.z() = 1.0/L.z();
	    
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
          
          /**
           * @brief Calculates optimal parameters for the Ewald method. It will only calculate a value if the user has not defined it in the input.
           * @param delta_in Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$)). Should not be larger than \f$ 5\cdot 10^{-5} \f$.
           * @warning Need to call 'kVectorChange' after change of number of k-vectors
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
                calcRcIon();
              if(!useIonIon && !useIonDipole && useDipoleDipole)
                calcRcDipole();
            }
            if(!user_kc) {
              if(useIonIon || useIonDipole)
                calcKcIon();
              if(!useIonIon && !useIonDipole && useDipoleDipole)
                calcKcDipole();
            }
          }
      };
      
  }//namespace

  namespace Energy {

    using namespace Faunus::Potential;

    /**
     * @brief Ewald summation for electrostatic interactions
     * @date Lund 2014
     *
     * This will take care of long-ranged electrostatic interactions using the Ewald summation scheme Currently it is possible to include ions and dipoles. 
     * 
     *  Keyword          |  Description
     * :--------------   | :---------------
     * `delta`           |  Accuracy in system for energy (\f$ e^2/\AA \f$).                                                                 (Default: 5e-5)
     * `eps_surf`        |  Dielectric constant of the surroundings                                                                          (Default: \f$ \varepsilon_r = \infty \f$)
     * `spherical_sum`   |  Spherical (or ellipsoid) summation in reciprocal space if set to be true. Cubic summation is applied if false.   (Default: true)
     * `cutoff`          |  Real space cut-off.                                                                                              (Default: Half the minimum box-length)              
     * `alpha`           |  Damping parameter.                                                                                               (Default: According to DOI: (Ions) 10.1080/08927029208049126 or (Dipoles) 10.1063/1.1398588 )  
     * `cutoffK`         |  Maximum number of vectors in any axis in k-space.                                                                (Default: According to DOI: (Ions) 10.1080/08927029208049126 or (Dipoles) 10.1063/1.1398588 )  
     * `cutoffK_x`       |  Maximum number of vectors in x-axis in k-space for ions. Is overridden if 'cutoffK' is set.                      (Default: According to DOI: (Ions) 10.1080/08927029208049126 or (Dipoles) 10.1063/1.1398588 )  
     * `cutoffK_y`       |  Maximum number of vectors in y-axis in k-space for ions. Is overridden if 'cutoffK' is set.                      (Default: According to DOI: (Ions) 10.1080/08927029208049126 or (Dipoles) 10.1063/1.1398588 )  
     * `cutoffK_z`       |  Maximum number of vectors in z-axis in k-space for ions. Is overridden if 'cutoffK' is set.                      (Default: According to DOI: (Ions) 10.1080/08927029208049126 or (Dipoles) 10.1063/1.1398588 )  
     * `update_frequency`|  The frequency of how often the total sum of all complex numbers are updated (an optimization optin).             (Default: Number of particles in system)
     * 
     * @note Tested and implemented through DOI: 10.1063/1.481216
     * @note If parameters are not set; optimized parameters for ion-ion-interactions will be used for all interactions save the case when only dipoles are present, then optimized parameters for dipole-dipole-interactions will be used. This is done since the wave-functions for ions and dipoles needs to be compatible with each other in order to get ion-dipole interactions.
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
     * where 
     * 
     * @f[
     * {\bf k} = 2\pi\left( \frac{n_x}{L_x} , \frac{n_y}{L_y} ,\frac{n_z}{L_z} \right)  \;\;,\;\; {\bf n} \in \mathbb{Z}^3
     * @f]
     * 
     * @todo Need to implement measurement of the average evaluation time for real part of Ewald, this in order to get optimal parametrization to work properly.
     * 
     */
    template<
      class Tspace, \
      class Tpairpot, \
      bool useIonIon=true, bool useIonDipole=false, bool useDipoleDipole=false, \
      class Tbase=NonbondedVector<Tspace, \
      CombinedPairPotential<EwaldReal<useIonIon,useIonDipole,useDipoleDipole>, Tpairpot>>>

      class NonbondedEwald : public Tbase {
        private:
          using Tbase::spc;
          typedef typename Tbase::Tparticle Tparticle;
          typedef typename Tbase::Tpvec Tpvec;

	  EwaldParameters<useIonIon,useIonDipole,useDipoleDipole> parameters, parameters_trial;
          Average<double> timeReciprocal;
          int kVectorsInUse, kVectorsInUse_trial, N, updates, cnt_accepted, update_frequency;
          double V, V_trial, eps_surf, eps_r, lB, const_inf, delta, Sq2, Smu2, selfEnergy, selfEnergyTrial, surfaceEnergy, surfaceEnergyTrial, reciprocalEnergy, reciprocalEnergyTrial; 
          bool spherical_sum, update_frequency_bool;
          vector<complex<double>> Q_ion_tot, Q_dip_tot, Q_ion_tot_trial, Q_dip_tot_trial;
          typename Tspace::Change change;

          Eigen::MatrixXd kVectors, kVectors_trial;  // Matrices with k-vectors
          Eigen::VectorXd Aks, Aks_trial;  // Stores values based on k-vectors in order to minimize computational effort. (See Eq.24 in DOI: 10.1063/1.481216)

          /**
           * @brief Returns Ewald self energy (kT) for ions and dipoles.
           * @param p Vector with partciles
           * @param g Group to calculate from
           * @note Does only add to average value if all atoms are counted
           */
          template<class Tpvec, class Tgroup>
            double getSelfEnergy(const Tpvec &p, const Tgroup &g, const EwaldParameters<useIonIon,useIonDipole,useDipoleDipole> &parameters_in) { 
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
              return (-parameters_in.alpha*( Sq2 + parameters_in.alpha2*(2.0/3.0)*Smu2 ) / sqrt(pc::pi))*lB;
            }

          /**
           * @brief Returns Ewald surface energy (kT) for ions and dipoles.
           * @param p Vector with partciles
           * @param g Group to calculate from
           */
          template<class Tpvec, class Tgroup>
            double getSurfaceEnergy(const Tpvec &p, const Tgroup &g, double V_in) { 
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
              return const_inf * (2*pc::pi/(( 2*eps_surf + 1)*V_in))*( qrs.dot(qrs) + 2*qrs.dot(mus) +  mus.dot(mus) )*lB;
            }

          /**
           * @brief Returns reciprocal space energy (kT).
           * @param p Vector with partciles
           * @param g Group to calculate from
           */
          double getReciprocalEnergy(const vector<complex<double>> &Q_ion_tot_in, const vector<complex<double>> &Q_dip_tot_in, const Eigen::VectorXd &Aks_in, double V_in) const {
            double E = 0.0;
            for (size_t k=0; k < Q_ion_tot_in.size(); k++) {
              double Q2 = 0.0;
              if(useIonIon)
                Q2 += std::norm(Q_ion_tot_in.at(k));  // 'norm' returns squared magnitude
              if(useIonDipole)
                Q2 += 2.0*std::real(Q_ion_tot_in.at(k)*Q_dip_tot_in.at(k));
              if(useDipoleDipole)
                Q2 += std::norm(Q_dip_tot_in.at(k));  // 'norm' returns squared magnitude
              E += ( Aks_in[k] * Q2 );
            }
            return (2*pc::pi/V_in)*E*lB;
          }

          /**
           * @brief Updates all vectors and matrices which depends on the number of k-space vectors.
           * @note Needs to be called whenever 'kcc_x', 'kcc_y' or 'kcc_z' has been updated
           */
          void kVectorChange(const EwaldParameters<useIonIon,useIonDipole,useDipoleDipole> &parameters_in, Eigen::MatrixXd &kVectors_in, Eigen::VectorXd &Aks_in, vector<complex<double>> &Q_ion_tot_in, vector<complex<double>> &Q_dip_tot_in, int &kVectorsInUse_in ) {
            int kVectorsLength = (2*parameters_in.kcc_x + 1)*(2*parameters_in.kcc_y + 1)*(2*parameters_in.kcc_z + 1) - 1;
	    if(kVectorsLength == 0) {
	      kVectors_in.resize(3, 1); 
	      Aks_in.resize(1);
	      kVectors_in.col(0) = Point(1.0,0.0,0.0); // Just so it is not the zero-vector
	      Aks_in[0] = 0.0;
	      kVectorsInUse_in = 1;
	      Q_ion_tot_in.resize(1);
	      Q_dip_tot_in.resize(1);
	      return;
	    }
            kVectors_in.resize(3, kVectorsLength); 
            Aks_in.resize(kVectorsLength);
            kVectorsInUse_in = 0;
            kVectors_in.setZero();
            Aks_in.setZero();

            double factor = 1.0;
            for (int kx = 0; kx <= parameters_in.kcc_x; kx++) {
              if(kx > 0)
                factor = 2.0;
              double dkx2 = double(kx*kx);
              for (int ky = -parameters_in.kcc_y; ky <= parameters_in.kcc_y; ky++) {
                double dky2 = double(ky*ky);
                for (int kz = -parameters_in.kcc_z; kz <= parameters_in.kcc_z; kz++) {
                  double dkz2 = double(kz*kz);
                  Point kv = 2*pc::pi*Point(kx*parameters_in.Li.x(),ky*parameters_in.Li.y(),kz*parameters_in.Li.z());
                  double k2 = kv.dot(kv);
                  if (k2 < parameters_in.check_k2_zero) { // Check if k2 != 0
                    continue;
                  }
                  if(spherical_sum) {
                    if( (dkx2/parameters_in.kc2_x) + (dky2/parameters_in.kc2_y) + (dkz2/parameters_in.kc2_z) > 1.0)
                      continue;
                  }
                  kVectors_in.col(kVectorsInUse_in) = kv; 
                  Aks_in[kVectorsInUse_in] = factor*exp(-k2/(4.0*parameters_in.alpha2))/k2;
                  kVectorsInUse_in++;
                }
              }
            }
            Q_ion_tot_in.resize(kVectorsInUse_in);
            Q_dip_tot_in.resize(kVectorsInUse_in);
          }

          /**
           * @brief Calculates optimal parameters for the Ewald method. It will only calculate a value if the user has not defined it in the input.
           * @param delta_in Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$)). Should not be larger than \f$ 5\cdot 10^{-5} \f$.
           * @warning Need to call 'kVectorChange' after change of number of k-vectors
           */
          void calcParameters(double delta_in=5e-5) {
            parameters.calcParameters(delta_in);
            Tbase::pairpot.first.updateAlpha(parameters.alpha);
            Tbase::pairpot.first.updateRcut(parameters.rc);
          }

          string _info() override {
            using namespace Faunus::textio;
            char w=25;
            std::ostringstream o;
            o << Tbase::_info();
            o << header("Ewald summation");
            string text = "";
            if(useIonIon)
              text += "Ion-Ion, ";
            if(useIonDipole)
              text += "Ion-dipole, ";
            if(useDipoleDipole)
              text += "Dipole-Dipole, ";
            if(useIonIon || useIonDipole || useDipoleDipole)
              text.erase (text.end()-2, text.end());
            o << pad(SUB,w, "Interactions") << text << endl;
            o << pad(SUB,w, "Reciprocal cut-off") << "(" << parameters.kc_x << "," << parameters.kc_y << "," << parameters.kc_z << ")" << endl;
	    if(parameters.kcc_x == 0 && parameters.kcc_y == 0 && parameters.kcc_z == 0) {
	      o << pad(SUB,w, "Wavefunctions") << 0 << endl;
	    } else {
	      o << pad(SUB,w, "Wavefunctions") << kVectorsInUse << endl;
	    }
            o << pad(SUB,w, "alpha") << parameters.alpha << endl;
            o << pad(SUB,w, "Real cut-off") << parameters.rc << endl;
            if(const_inf < 0.5) {
              o << pad(SUB,w+1, epsilon_m+"(Surface)") << infinity << endl;
            } else {
              o << pad(SUB,w+1, epsilon_m+"(Surface)") << eps_surf << endl;
            }
            o << pad(SUB,w, "updates") << updates << endl;
	    
            o << pad(SUB,w+4, bracket("Energies / kT")) << endl;
            o << pad(SUB,w+4, bracket("Self energy")) << selfEnergyAverage.avg() << endl;
            o << pad(SUB,w+1, "    "+sigma) << selfEnergyAverage.std() << endl;
            o << pad(SUB,w+4, bracket("Surf energy")) << surfaceEnergyAverage.avg() << endl;
            o << pad(SUB,w+1, "    "+sigma) << surfaceEnergyAverage.std() << endl;
            o << pad(SUB,w+4, bracket("Real energy")) << realEnergyAverage.avg() << endl;
            o << pad(SUB,w+1, "    "+sigma) << realEnergyAverage.std() << endl;
            o << pad(SUB,w+4, bracket("Reci energy")) << reciprocalEnergyAverage.avg() << endl;
            o << pad(SUB,w+1, "    "+sigma) << reciprocalEnergyAverage.std() << endl;
            
            return o.str();
          }

        public:
          MeanValue<double> selfEnergyAverage, surfaceEnergyAverage, realEnergyAverage, reciprocalEnergyAverage;

          NonbondedEwald(Tmjson &j, const string &sec="nonbonded") : Tbase(j,sec) , selfEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) ,surfaceEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) , realEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) , reciprocalEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100)  {
            Tbase::name += " (Ewald)";
	    parameters.update(1.0,timeReciprocal.avg(),N,Sq2,Smu2,typename Tspace::GeometryType(j).len);
            auto _j = j["energy"]["nonbonded"]["ewald"];
            cnt_accepted = 0;
            updates = 0;           // will finally have the value of how many time setSpace() has been used
            lB = Tbase::pairpot.first.bjerrumLength();
            eps_surf = ( _j["eps_surf"] | 0.0 );
            const_inf = 1.0;    
            if(eps_surf < 1)  // if the value is unphysical (< 1) then we set infinity as the dielectric sonatant of the surronding medium
              const_inf = 0.0; // \varepsilon_{Surface} = \infty

            spherical_sum = _j["spherical_sum"] | true; // specifies if Spherical or Cubical summation should be used in reciprocal space
            parameters.delta = ( _j["delta"] | 5e-5 );
            if( ( _j["update_frequency"] | -1 )  > 0 ) {  // specifies after how many successful moves all complex numbers should the updated
              update_frequency_bool = true;
              update_frequency = _j["update_frequency"] | -1;
            } else {
              update_frequency_bool = false;
            }

            // Initiate alpha-values
            parameters.user_alpha = false;
            if ( ( _j["alpha"] | -1.0 ) > -1e-9 ) {  // set damping-parameter
              parameters.user_alpha = true;
              parameters.alpha = ( _j["alpha"] | -1.0 );
	      assert(!(parameters.alpha < 0) && "alpha parameter cannot be negative!");
              parameters.alpha2 = parameters.alpha*parameters.alpha;
            }

            // Initiate real cut-off
            parameters.user_rc = false;
            if ( ( _j["cutoff"] | -1.0 ) > 0.0 ) {  // set cutoff
              parameters.user_rc = true;
              parameters.rc = ( _j["cutoff"] | -1.0 );
              if (parameters.rc > parameters.minL/2.0)
                parameters.rc = parameters.minL/2.0;
            }
            // Initiate reciprocal cut-off
            parameters.user_kc = false;
            if ( ( _j["cutoffK"] | -1.0 ) > -0.5 ) {  // set reciprocal cutoff
              parameters.user_kc = true;
              parameters.kc_x = ( _j["cutoffK"] | -1.0 );
              parameters.kc_y = parameters.kc_x;
              parameters.kc_z = parameters.kc_x;
            } else {
              if( ( _j["cutoffK_x"] | -1.0 ) > 0.0 || ( _j["cutoffK_y"] | -1.0 ) > 0.0 || ( _j["cutoffK_z"] | -1.0 ) > 0.0) {
                parameters.user_kc = true;
                parameters.kc_x = ( _j["cutoffK_x"] | -1.0 );
                parameters.kc_y = ( _j["cutoffK_y"] | parameters.kc_x );
                parameters.kc_z = ( _j["cutoffK_z"] | parameters.kc_y );

                // Assures that kc_x_ion, kc_y_ion and kc_z_ion willhave non-negative values
                if(parameters.kc_x < 0.0) {
                  if(parameters.kc_y < 0.0) {
                    parameters.kc_x = parameters.kc_z;
                    parameters.kc_y = parameters.kc_z;
                  } else {
                    parameters.kc_x = parameters.kc_y;
                    if(parameters.kc_z < 0.0)
                      parameters.kc_z = parameters.kc_y;
                  }
                } else {
                  if(parameters.kc_y < 0.0)
                    parameters.kc_y = parameters.kc_x;
                  if(parameters.kc_z < 0.0)
                    parameters.kc_z = parameters.kc_x;
                }
              }
            }
            parameters.kc2_x = parameters.kc_x*parameters.kc_x;
            parameters.kc2_y = parameters.kc_y*parameters.kc_y;
            parameters.kc2_z = parameters.kc_z*parameters.kc_z;
            parameters.kcc_x = ceil(parameters.kc_x);
            parameters.kcc_y = ceil(parameters.kc_y);
            parameters.kcc_z = ceil(parameters.kc_z);
	    
	    calcParameters(parameters.delta);
	    parameters_trial = parameters;
            kVectorChange(parameters,kVectors,Aks,Q_ion_tot,Q_dip_tot,kVectorsInUse);
	    kVectorChange(parameters_trial,kVectors_trial,Aks_trial,Q_ion_tot_trial,Q_dip_tot_trial,kVectorsInUse_trial);
            change.clear();

          }

          double i_external(const Tpvec &p, int i) override {
            Group g = Group(i,i);
            return g_external(p,g);
          }

          double g_external(const Tpvec &p, Group &g) override {
	    if (Tbase::isTrial(p))
	      return getSelfEnergy(p,g,parameters_trial);
	    return getSelfEnergy(p,g,parameters);
          }

          /**
           * @brief Updates vectors of complex numbers if a move has been accepted. These vectors are implemented in order to increase computational speed. 
           * After 'update_frequency' number of non-isobaric updates the entirety of the complex vectors are recalculated in order to avoid numerical errors.
           * @param move_accepted True if a move is accepted, otherwise false
           * 
           * @note Needs to be called before particle-vectors are updated after acception/rejection.
           */
          double update(bool move_accepted) override {
            if (!move_accepted ) {
	      selfEnergyAverage += selfEnergy;
	      surfaceEnergyAverage += surfaceEnergy;
	      reciprocalEnergyAverage += reciprocalEnergy;
              return 0.0;
	    }

            assert(!change.empty());
	    
	    selfEnergyAverage += selfEnergyTrial;
	    surfaceEnergyAverage += surfaceEnergyTrial;
	    reciprocalEnergyAverage += reciprocalEnergyTrial;
	    
	    double du = 0.0;
	    
	    if(++cnt_accepted > update_frequency) {
	      du -= getReciprocalEnergy(Q_ion_tot_trial,Q_dip_tot_trial,Aks_trial,V_trial);
	      updateAllComplexNumbers(spc->trial, Q_ion_tot_trial, Q_dip_tot_trial,kVectors_trial,kVectorsInUse_trial);
	      du += getReciprocalEnergy(Q_ion_tot_trial,Q_dip_tot_trial,Aks_trial,V_trial);
	      cnt_accepted = 0;
	    }
	    
	    V = V_trial;
            kVectorsInUse = kVectorsInUse_trial;
            Q_ion_tot.resize(kVectorsInUse);
            Q_dip_tot.resize(kVectorsInUse);
            kVectors.resize(3, kVectorsInUse); 
            for (int k=0; k < kVectorsInUse; k++) {
              Q_ion_tot.at(k) = Q_ion_tot_trial.at(k);
              Q_dip_tot.at(k) = Q_dip_tot_trial.at(k);
              kVectors.col(k) = kVectors_trial.col(k);
              Aks[k] = Aks_trial[k];
            }
            return du;
          }

          /** @brief Update energy function due to Change 
           */
          void setChange(const typename Tspace::Change &c) override {
            change = c;

            if(c.geometryChange) {
              updateAllComplexNumbers(spc->trial, Q_ion_tot_trial, Q_dip_tot_trial,kVectors_trial,kVectorsInUse_trial);
              V_trial = V + c.dV;
              return;
            }

            for (int k=0; k<kVectorsInUse_trial; k++) {
              complex<double> Q2_ion = Q_ion_tot.at(k);
              complex<double> Q2_dip = Q_dip_tot.at(k);
              for (auto m : change.mvGroup) {
                for (auto i : m.second) {
                  double dotTrial = kVectors_trial.col(k).dot(spc->trial[i]);
                  double dot = kVectors_trial.col(k).dot(spc->p[i]);
                  if (useIonIon || useIonDipole) {
                    Q2_ion += spc->trial[i].charge * complex<double>(cos(dotTrial),sin(dotTrial));
                    Q2_ion -= spc->p[i].charge * complex<double>(cos(dot),sin(dot));
                  }
#ifdef DIPOLEPARTICLE
                  if (useDipoleDipole || useIonDipole) {
                    Q2_dip += kVectors_trial.col(k).dot(spc->trial[i].mu) * spc->trial[i].muscalar * complex<double>(-sin(dotTrial),cos(dotTrial));
                    Q2_dip -= kVectors_trial.col(k).dot(spc->p[i].mu) * spc->p[i].muscalar * complex<double>(-sin(dot),cos(dot));
                  }
#endif
                }
              }
              Q_ion_tot_trial.at(k) = Q2_ion;
              Q_dip_tot_trial.at(k) = Q2_dip;
            }
          }

          double external(const Tpvec &p) override {
            Group g(0, p.size()-1);
            std::chrono::steady_clock::time_point init, final;
            init = std::chrono::steady_clock::now();
	    double total = 0.0;
            if (Tbase::isTrial(p)) {
              surfaceEnergyTrial = getSurfaceEnergy(p,g,V_trial);
              reciprocalEnergyTrial = getReciprocalEnergy(Q_ion_tot_trial,Q_dip_tot_trial,Aks_trial,V_trial);
	      total = surfaceEnergyTrial + reciprocalEnergyTrial;
            } else {
              surfaceEnergy = getSurfaceEnergy(p,g,V);
              reciprocalEnergy = getReciprocalEnergy(Q_ion_tot,Q_dip_tot,Aks,V);
	      total = surfaceEnergy + reciprocalEnergy;
            }
            final = std::chrono::steady_clock::now();
            timeReciprocal += double(std::chrono::duration_cast<std::chrono::milliseconds >(final-init).count());
            return total;
          }

          /**
           * @brief Re-calculates the vectors of complex numbers used in getQ2
           * @param p Particle vector
           */
          void updateAllComplexNumbers(const Tpvec &p, vector<complex<double>> &Q_ion_tot_in, vector<complex<double>> &Q_dip_tot_in, const Eigen::MatrixXd &kVectors_in, int kVectorsInUse_in) {
            for (int k=0; k<kVectorsInUse_in; k++) {
              Point kv = kVectors_in.col(k);
              complex<double> Q_temp_ion(0.0,0.0);
              complex<double> Q_temp_dip(0.0,0.0);
              for (size_t i = 0; i < p.size(); i++) {
                double dot = kv.dot(p[i]);
                if(useIonIon || useIonDipole)
                  Q_temp_ion += p[i].charge * complex<double>(cos(dot),sin(dot));
#ifdef DIPOLEPARTICLE 
                if(useDipoleDipole)
                  Q_temp_dip += kv.dot(p[i].mu) * p[i].muscalar * complex<double>(-sin(dot),cos(dot));
#endif
              }
              Q_ion_tot_in.at(k) = Q_temp_ion;
              Q_dip_tot_in.at(k) = Q_temp_dip;
            }
          }

          /**
           * @brief Set space and updates parameters (if not set by user)
	   * @warning Sets all trial-entities as the original ones 
           */
          void setSpace(Tspace &s) override {
            Tbase::setSpace(s);
            N = s.p.size();
            V = spc->geo.getVolume();
	    V_trial = V;
            calcParameters(delta);
            if(change.geometryChange || !parameters.user_kc)
              kVectorChange(parameters,kVectors,Aks,Q_ion_tot,Q_dip_tot,kVectorsInUse);
            updateAllComplexNumbers(s.p,Q_ion_tot,Q_dip_tot,kVectors,kVectorsInUse);
            if(!update_frequency_bool)
              update_frequency = N;
            updates++;
          }
      };

  }//namespace
}//namespace
#endif

