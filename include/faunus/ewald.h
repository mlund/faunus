#ifndef FAUNUS_EWALD_H
#define FAUNUS_EWALD_H

#include "faunus/inputfile.h"
#include <faunus/ewald.h>
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
    template<bool useIonIon=true, bool useIonDipole=true, bool useDipoleDipole=true> 
      struct EwaldReal : public Potential::Coulomb {

        typedef Potential::Coulomb Tbase;
        double alpha_ion, alpha2_ion, alpha_dip, alpha2_dip, constant_ion, constant_dip, rc_ion, rc2i_ion, rc_dip, rc2i_dip;

        EwaldReal(InputMap &in, string pfx="coulomb") : Tbase(in,pfx), alpha_ion(0), alpha_dip(0), rc2i_ion(0), rc2i_dip(0) {
          Tbase::name="Ewald Real";
        }

        void updateAlpha(double alpha_ion_in, double alpha_dip_in=-1.0) {
	  if(alpha_dip_in < 0)
	    alpha_dip_in = alpha_ion_in;
          if(useIonIon || useIonDipole) {
            alpha_ion = alpha_ion_in;
            alpha2_ion = alpha_ion*alpha_ion;
            constant_ion = 2*alpha_ion/sqrt(pc::pi);
          }
          if(useDipoleDipole) {
            alpha_dip = alpha_dip_in;
            alpha2_dip = alpha_dip*alpha_dip;
            constant_dip = 2*alpha_dip/sqrt(pc::pi);
          }
        }
        
        void updateRcut(double rc_ion_in, double rc_dip_in=-1.0) {
	  if(rc_dip_in < 0)
	    rc_dip_in = rc_ion_in;
          if(useIonIon || useIonDipole) {
            rc2i_ion = 1.0/(rc_ion_in*rc_ion_in);
	    rc_ion = rc_ion_in;
	  }
          if(useDipoleDipole) {
            rc2i_dip = 1.0/(rc_dip_in*rc_dip_in);
	    rc_dip = rc_dip_in;
	  }
        }

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
	    
	    // No interaction
	    if(!useIonIon && !useIonDipole && !useDipoleDipole) {
	      return 0.0;
	    }

	    // Ion-Ion interaction
	    if(useIonIon && !useIonDipole && !useDipoleDipole) {
	      double r2i = 1.0/r.squaredNorm();
	      if (r2i < rc2i_ion)
		return 0.0;
	      double r1i = sqrt(r2i);
              return erfc_x(alpha_ion/r1i)*r1i*a.charge*b.charge * Tbase::bjerrumLength();
	    }
	    
#ifdef DIPOLEPARTICLE   

	    // Ion-Ion and Ion-Dipole interaction
	    if(useIonIon && useIonDipole && !useDipoleDipole) {
	      double r2i = 1.0/r.squaredNorm();
	      if (r2i < rc2i_ion)
		return 0.0;
	      double r1i = sqrt(r2i);
	      double r1i_d_ion = erfc_x(alpha_ion/r1i)*r1i;
	      double E = r1i_d_ion*a.charge*b.charge;
              double T1_ion = (constant_ion*exp(-alpha2_ion/r2i) + r1i_d_ion)*r2i;
              E += a.charge*r.dot(b.mu)*b.muscalar*T1_ion;
              E -= b.charge*r.dot(a.mu)*a.muscalar*T1_ion;
	      return E * Tbase::bjerrumLength();
	    }
	    
	    // Ion-Ion and Dipole-Dipole interaction
	    if(useIonIon && !useIonDipole && useDipoleDipole) {
	      double r1 = r.norm();
	      double E = 0.0;
	      if (r1 < rc_ion) {
		E += erfc_x(alpha_ion*r1)*a.charge*b.charge/r1; // Add ion-ion
	      }
	      if (r1 < rc_dip) {
		double r2i = 1.0/(r1*r1);
		double r1i_d_dip = erfc_x(alpha_dip*r1)/r1;
		double expK_dip = constant_dip*exp(-alpha2_dip/r2i);
		double t3 = -b.mu.dot(a.mu)*(expK_dip + r1i_d_dip)*r2i;
		double t5 = b.mu.dot(r)*a.mu.dot(r)*(3.*r1i_d_dip*r2i + (3.*r2i + 2.*alpha2_dip)*expK_dip)*r2i;
		E += -(t5 + t3)*b.muscalar*a.muscalar; // Add dipole-dipole
	      }
	      return E * Tbase::bjerrumLength();
	    }
	    
	    // Ion-Ion, Ion-Dipole and Dipole-Dipole interaction
	    if(useIonIon && useIonDipole && useDipoleDipole) {
	      double r2i = 1.0/r.squaredNorm();
	      double E = 0.0;
	      double r1i = sqrt(r2i);
	      if (r2i > rc2i_ion) {
		double r1i_d_ion = erfc_x(alpha_ion/r1i)*r1i;
		E += r1i_d_ion*a.charge*b.charge; // Add ion-ion
		double T1_ion = (constant_ion*exp(-alpha2_ion/r2i) + r1i_d_ion)*r2i;
		E += a.charge*r.dot(b.mu)*b.muscalar*T1_ion;
		E -= b.charge*r.dot(a.mu)*a.muscalar*T1_ion; // Add ion-dipole
	      }
	      if (r2i > rc2i_dip) {
		double r1i_d_dip = erfc_x(alpha_dip/r1i)*r1i;
		double expK_dip = constant_dip*exp(-alpha2_dip/r2i);
		double t3 = -b.mu.dot(a.mu)*(expK_dip + r1i_d_dip)*r2i;
		double t5 = b.mu.dot(r)*a.mu.dot(r)*(3.*r1i_d_dip*r2i + (3.*r2i + 2.*alpha2_dip)*expK_dip)*r2i;
		E += -(t5 + t3)*b.muscalar*a.muscalar; // Add dipole-dipole
	      }
	      return E * Tbase::bjerrumLength();
	    }
	    
	    // Ion-Dipole interaction
	    if(!useIonIon && useIonDipole && !useDipoleDipole) {
	      double r2i = 1.0/r.squaredNorm();
	      if (r2i < rc2i_ion)
		return 0.0;
	      double r1i = sqrt(r2i);
	      double T1_ion = (constant_ion*exp(-alpha2_ion/r2i) + erfc_x(alpha_ion/r1i)*r1i)*r2i;
	      double E = a.charge*r.dot(b.mu)*b.muscalar*T1_ion;
	      E -= b.charge*r.dot(a.mu)*a.muscalar*T1_ion; // Add ion-dipole
	      return E * Tbase::bjerrumLength();
	    }
	    
	    // Dipole-Dipole interaction
	    if(!useIonIon && !useIonDipole && useDipoleDipole) {
	      double r2i = 1.0/r.squaredNorm();
	      if (r2i < rc2i_dip)
		return 0.0;
	      double r1i = sqrt(r2i);
	      double r1i_d_dip = erfc_x(alpha_dip/r1i)*r1i;
	      double expK_dip = constant_dip*exp(-alpha2_dip/r2i);
	      double t3 = -b.mu.dot(a.mu)*(expK_dip + r1i_d_dip)*r2i;
	      double t5 = b.mu.dot(r)*a.mu.dot(r)*(3.*r1i_d_dip*r2i + (3.*r2i + 2.*alpha2_dip)*expK_dip)*r2i;
	      return -(t5 + t3)*b.muscalar*a.muscalar * Tbase::bjerrumLength(); // Add dipole-dipole
	    }
	    
	    // Ion-Dipole and Dipole-Dipole interaction
	    if(!useIonIon && useIonDipole && useDipoleDipole) {
	      double r2i = 1.0/r.squaredNorm();
	      double E = 0.0;
	      double r1i = sqrt(r2i);
	      if (r2i > rc2i_ion) {
		double T1_ion = (constant_ion*exp(-alpha2_ion/r2i) + erfc_x(alpha_ion/r1i)*r1i)*r2i;
		double E = a.charge*r.dot(b.mu)*b.muscalar*T1_ion;
		E -= b.charge*r.dot(a.mu)*a.muscalar*T1_ion; // Add ion-dipole
	      }
	      if (r2i > rc2i_dip) {
		double r1i_d_dip = erfc_x(alpha_dip/r1i)*r1i;
		double expK_dip = constant_dip*exp(-alpha2_dip/r2i);
		double t3 = -b.mu.dot(a.mu)*(expK_dip + r1i_d_dip)*r2i;
		double t5 = b.mu.dot(r)*a.mu.dot(r)*(3.*r1i_d_dip*r2i + (3.*r2i + 2.*alpha2_dip)*expK_dip)*r2i;
		E += -(t5 + t3)*b.muscalar*a.muscalar; // Add dipole-dipole
	      }
	      return E * Tbase::bjerrumLength();
	    }
#endif
	  return 0.0;
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
     * `cutoff_ion`      |  Real space cut-off for ions. Is overridden if 'cutoff' is set.                                                   (Default: Half the minimum box-length)       
     * `cutoff_dip`      |  Real space cut-off for dipoles. Is overridden if 'cutoff' is set.                                                (Default: Half the minimum box-length)       
     * `alpha`           |  Damping parameter.                                                                                               (Default: According to DOI: 10.1080/08927029208049126 )  
     * `alpha_ion`       |  Damping parameter for ions. Is overridden if 'alpha' is set.                                                     (Default: According to DOI: 10.1080/08927029208049126 )
     * `alpha_dip`       |  Damping parameter for dipoles. Is overridden if 'alpha' is set.                                                  (Default: According to DOI: 10.1063/1.1398588 )
     * `cutoffK`         |  Maximum number of vectors in any axis in k-space.                                                                (Default: According to DOI: 10.1080/08927029208049126 )
     * `cutoffK_ion`     |  Maximum number of vectors in any axis in k-space for ions. Is overridden if 'cutoffK' is set.                    (Default: According to DOI: 10.1080/08927029208049126 )
     * `cutoffK_x_ion`   |  Maximum number of vectors in x-axis in k-space for ions. Is overridden if 'cutoffK' or 'cutoffK_ion' is set.     (Default: According to DOI: 10.1080/08927029208049126 )
     * `cutoffK_y_ion`   |  Maximum number of vectors in y-axis in k-space for ions. Is overridden if 'cutoffK' or 'cutoffK_ion' is set.     (Default: According to DOI: 10.1080/08927029208049126 )
     * `cutoffK_z_ion`   |  Maximum number of vectors in z-axis in k-space for ions. Is overridden if 'cutoffK' or 'cutoffK_ion' is set.     (Default: According to DOI: 10.1080/08927029208049126 )
     * `cutoffK_dip`     |  Maximum number of vectors in any axis in k-space for dipoles. Is overridden if 'cutoffK' is set.                 (Default: According to DOI: 10.1063/1.1398588 )
     * `cutoffK_x_dip`   |  Maximum number of vectors in x-axis in k-space for dipoles. Is overridden if 'cutoffK' or 'cutoffK_dip' is set.  (Default: According to DOI: 10.1063/1.1398588 )
     * `cutoffK_y_dip`   |  Maximum number of vectors in y-axis in k-space for dipoles. Is overridden if 'cutoffK' or 'cutoffK_dip' is set.  (Default: According to DOI: 10.1063/1.1398588 )
     * `cutoffK_z_dip`   |  Maximum number of vectors in z-axis in k-space for dipoles. Is overridden if 'cutoffK' or 'cutoffK_dip' is set.  (Default: According to DOI: 10.1063/1.1398588 )
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

          Average<double> timeReciIonIon, timeRealIonIon, timeReciDipDip, timeRealDipDip;

          int kcc_x_ion, kcc_y_ion, kcc_z_ion, kcc_x_dip, kcc_y_dip, kcc_z_dip;
	  int kVectorsLength_ion, kVectorsLength_dip, wavefunctions_ion, wavefunctions_dip, N, updates;
          double kc_x_ion, kc2_x_ion, kc_y_ion, kc2_y_ion, kc_z_ion, kc2_z_ion, kc_x_dip, kc2_x_dip, kc_y_dip, kc2_y_dip, kc_z_dip, kc2_z_dip; 
	  double alpha_ion, alpha2_ion, alpha_dip, alpha2_dip, V, eps_surf, eps_r, rc_ion, rc_dip, lB, const_inf, check_k2_zero, delta, minL, maxL,Sq2, Smu2;
          Point Li, L;                                                                                                     // Inverse box-length-vector, Box-length-vector
          bool user_alpha_ion, user_alpha_dip, user_kc_ion, user_kc_dip, user_rc_ion, user_rc_dip, spherical_sum; 
          vector<double> values_alpha_ion, values_alpha_dip, values_kcc_x_ion, values_kcc_y_ion, values_kcc_z_ion, values_kcc_x_dip, values_kcc_y_dip, values_kcc_z_dip, values_rc_ion, values_rc_dip, values_wavefunctions_ion, values_wavefunctions_dip;

          Eigen::MatrixXd kVectors_ion, kVectors_dip;          // Matrices with k-vectors for ions and dipoles respectively
          Eigen::VectorXd k2s_ion, k2s_dip, Aks_ion, Aks_dip;  // Stores values based on k-vectors in order to minimize computational effort.

          /**
           * @brief Returns Ewald self energy for ions and dipoles.
           * @param p Vector with partciles
           * @param g Group to calculate from
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
	      
              selfEA += -( alpha_ion * Eq + alpha_dip*alpha2_dip*(2.0/3.0)*Emu ) / sqrt(pc::pi);
              return -( alpha_ion * Eq + 2*alpha_dip*alpha2_dip/3*Emu ) / sqrt(pc::pi);
            }

          /**
           * @brief Returns Ewald surface energy for ions and dipoles.
           * @param p Vector with partciles
           * @param g Group to calculate from
           * @warning Have trouble with few (2,...) particles!
           */
          template<class Tpvec, class Tgroup>
            double getSurfEnergy(const Tpvec &p, const Tgroup &g) { // const
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
              surfEA += const_inf * (2*pc::pi/(( 2*eps_surf + 1)*V))*( qrs.dot(qrs) + 2*qrs.dot(mus) +  mus.dot(mus) );
              return const_inf * (2*pc::pi/(( 2*eps_surf + 1)*V)) * ( qrs.dot(qrs) + 2*qrs.dot(mus) +  mus.dot(mus) );
            }

          /**
           * @brief Support-function to 'getReciEnergy' (for ions)
           * @param p Vector with partciles
           * @param g Group to calculate from
           * @param kv Reciprocal space vector
           */
          template<class Tpvec, class Tgroup>
            double getQ2_ion(const Tpvec &p, const Tgroup &g, const Point &kv) {
              double dQ2 = 0.0;
              complex<double> Q2(0,0);
              for (auto i : g) {
                double dot = kv.dot(p[i]);
                Q2 += p[i].charge * complex<double>(cos(dot),sin(dot));
              }
              dQ2 = std::abs(Q2);
              return dQ2*dQ2;
            }

          /**
           * @brief Support-function to 'getReciEnergy' (for dipoles)
           * @param p Vector with partciles
           * @param g Group to calculate from
           * @param kv Reciprocal space vector
           */
          template<class Tpvec, class Tgroup>
            double getQ2_dip(const Tpvec &p, const Tgroup &g, const Point &kv) {
              double dQ2 = 0.0;
              complex<double> Q2(0,0);
              for (auto i : g) {
                double dot = kv.dot(p[i]);
                Q2 += kv.dot(p[i].mu) * p[i].muscalar * complex<double>(-sin(dot),cos(dot));
              }
              dQ2 = std::abs(Q2);
              return dQ2*dQ2;
            }

          /**
           * @brief Returns reciprocal space energy for ion-ion interactions.
           * @param p Vector with partciles
           * @param g Group to calculate from
           */
          template<class Tpvec, class Tgroup>
            double getReciEnergyIonIon(const Tpvec &p, const Tgroup &g) { // const
              double dQ2 = 0.0;
              double E = 0.0;

              if (useIonIon) {
                for (int i=0; i<wavefunctions_ion; i++) {
                  dQ2 = getQ2_ion(p,g,kVectors_ion.col(i));
                  E += ( Aks_ion[i] * dQ2 );
                }
	      } else {
		return 0.0;
	      }

              reciEA += (2*pc::pi/V)*E;
              return (2*pc::pi/V)*E;
            }
            
          /**
           * @brief Returns reciprocal space energy for ion-dipole interactions.
           * @param p Vector with partciles
           * @param g Group to calculate from
           * @warning Have to be fixed!
           */
          template<class Tpvec, class Tgroup>
            double getReciEnergyIonDip(const Tpvec &p, const Tgroup &g) { // const
              double dQ2 = 0.0;
              double E = 0.0;

              if (useIonDipole) {
                for (int i=0; i<wavefunctions_ion; i++) {
                  dQ2 = getQ2_ion(p,g,kVectors_ion.col(i));
                  E += ( Aks_ion[i] * dQ2 );
                }
	      } else {
		return 0.0;
	      }

              reciEA += (2*pc::pi/V)*E;
              return (2*pc::pi/V)*E;
            }
            
          /**
           * @brief Returns reciprocal space energy for dipole-dipole interactions.
           * @param p Vector with partciles
           * @param g Group to calculate from
           */
          template<class Tpvec, class Tgroup>
            double getReciEnergyDipDip(const Tpvec &p, const Tgroup &g) { // const
              double dQ2 = 0.0;
              double E = 0.0;

              if (useDipoleDipole) {
                for (int i=0; i<wavefunctions_dip; i++) {
                  dQ2 = getQ2_dip(p,g,kVectors_dip.col(i));
                  E += ( Aks_dip[i] * dQ2 );
                }
	      } else {
		return 0.0;
	      }

              reciEA += (2*pc::pi/V)*E;
              return (2*pc::pi/V)*E;
            }
            

          /**
           * @brief Updates all vectors and matrices which depends on the number of k-space vectors for ions
           * @note Needs to be called whenever 'kcc_x_ion', 'kcc_y_ion' or 'kcc_z_ion' has been updated
           */
          void kVectorChangeIon() {
            kVectorsLength_ion = (2*kcc_x_ion + 1)*(2*kcc_y_ion + 1)*(2*kcc_z_ion + 1) - 1;
            kVectors_ion.resize(3, kVectorsLength_ion); 
            k2s_ion.resize(kVectorsLength_ion);
            Aks_ion.resize(kVectorsLength_ion);
            wavefunctions_ion = 0;
            kVectors_ion.setZero();
            k2s_ion.setZero();
            Aks_ion.setZero();

            for (int kx = -kcc_x_ion; kx <= kcc_x_ion; kx++) {
	      double dkx2 = double(kx*kx);
              for (int ky = -kcc_y_ion; ky <= kcc_y_ion; ky++) {
		double dky2 = double(ky*ky);
                for (int kz = -kcc_z_ion; kz <= kcc_z_ion; kz++) {
		  double dkz2 = double(kz*kz);
                  Point kv = 2*pc::pi*Point(kx*Li.x(),ky*Li.y(),kz*Li.z());
                  double k2 = kv.dot(kv);
                  if (k2 < check_k2_zero) { // Check if k2 != 0
                    continue;
                  }
                  if(spherical_sum) {
		    if( (dkx2/kc2_x_ion) + (dky2/kc2_y_ion) + (dkz2/kc2_z_ion) > 1.0)
		      continue;
                  }
                  kVectors_ion.col(wavefunctions_ion) = kv; 
                  k2s_ion[wavefunctions_ion] = k2;
                  Aks_ion[wavefunctions_ion] = exp(-k2/(4*alpha2_ion))/k2;
                  wavefunctions_ion++;
                }
	      }
	    }
            values_wavefunctions_ion.push_back(double(wavefunctions_ion));
          }

          /**
           * @brief Updates all vectors and matrices which depends on the number of k-space vectors for dipoles
           * @note Needs to be called whenever 'kcc_x_dip', 'kcc_y_dip' or 'kcc_z_dip' has been updated
           */
          void kVectorChangeDip() {
            kVectorsLength_dip = (2*kcc_x_dip + 1)*(2*kcc_y_dip + 1)*(2*kcc_z_dip + 1) - 1;
            kVectors_dip.resize(3, kVectorsLength_dip); 
            k2s_dip.resize(kVectorsLength_dip);
            Aks_dip.resize(kVectorsLength_dip);
            wavefunctions_dip = 0;
            kVectors_dip.setZero();
            k2s_dip.setZero();
            Aks_dip.setZero();

            for (int kx = -kcc_x_dip; kx <= kcc_x_dip; kx++) {
	      double dkx2 = double(kx*kx);
              for (int ky = -kcc_y_dip; ky <= kcc_y_dip; ky++) {
		double dky2 = double(ky*ky);
                for (int kz = -kcc_z_dip; kz <= kcc_z_dip; kz++) {
		  double dkz2 = double(kz*kz);
                  Point kv = 2*pc::pi*Point(kx*Li.x(),ky*Li.y(),kz*Li.z());
                  double k2 = kv.dot(kv);
                  if (k2 < check_k2_zero) { // Check if k2 != 0
                    continue;
                  }
                  if(spherical_sum) {
		    if( (dkx2/kc2_x_dip) + (dky2/kc2_y_dip) + (dkz2/kc2_z_dip) > 1.0 )
		      continue;
                  }
                  kVectors_dip.col(wavefunctions_dip) = kv; 
                  k2s_dip[wavefunctions_dip] = k2;
                  Aks_dip[wavefunctions_dip] = exp(-k2/(4*alpha2_dip))/k2;
                  wavefunctions_dip++;
                }
	      }
	    }
            values_wavefunctions_dip.push_back(double(wavefunctions_dip));
          }

          /**
           * @brief Solves the equation \f$ w\cdot e^{w} = x \f$
           * @note Implemented through DOI: 10.1145/361952.361970
           */
          double LambertW_solver(double x, double error_max=1e-14) {
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
            double tau_F = timeReciIonIon.avg();
            alpha_ion = pow(tau_R*N/tau_F,1.0/6.0)*sqrt(pc::pi)/minL;
            values_alpha_ion.push_back(alpha_ion);
          }
          
          /**
           * @brief Calculates ionic real space cutoff parameter for the Ewald method.
           * @param alpha_ion_in Damping parameter (\f$ \AA^{-1} \f$)
           * @param delta_in Accuracy in ionic systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1080/08927029208049126  (DOI: 10.1007/b136795)
           */
          void calcRcIon(double alpha_ion_in, double delta_in = 5e-5) {
	    rc_ion = abs(0.5*sqrt(3.0*LambertW_solver( -(4.0/3.0)*pow(-pow(Sq2,4.0)/(alpha_ion_in*alpha_ion_in*pow(minL,6.0)*pow(delta_in,4.0)),1.0/3.0) )));
	    values_rc_ion.push_back(rc_ion);
	  }
          
          /**
           * @brief Calculates ionic reciprocal space cutoff parameter for the Ewald method.
           * @param alpha_ion_in Damping parameter (\f$ \AA^{-1} \f$)
           * @param delta_in Accuracy in ionic systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1080/08927029208049126  (DOI: 10.1007/b136795)
           */
          void calcKcIon(double alpha_ion_in, double delta_in = 5e-5) {
	    double temp = (1.0/3.0)*pow(2.0,2.0/3.0)*pow(4.0,1.0/3.0)*pow(pow(Sq2,4.0)/(alpha_ion_in*alpha_ion_in*pow(delta_in,4.0)),1.0/3.0);
	    kc_x_ion = abs(0.5*sqrt(3.0*LambertW_solver( temp*pow(L.x(), -2.0) )));
	    kc_x_ion = abs(0.5*sqrt(3.0*LambertW_solver( temp*pow(L.y(), -2.0) )));
	    kc_x_ion = abs(0.5*sqrt(3.0*LambertW_solver( temp*pow(L.z(), -2.0) )));
            kc2_x_ion = kc_x_ion*kc_x_ion;
	    kc2_y_ion = kc_y_ion*kc_y_ion;
	    kc2_z_ion = kc_z_ion*kc_z_ion;
            kcc_x_ion = ceil(kc_x_ion);
	    kcc_y_ion = ceil(kc_y_ion);
	    kcc_z_ion = ceil(kc_z_ion);
	    values_kcc_x_ion.push_back(kcc_x_ion);
	    values_kcc_y_ion.push_back(kcc_y_ion);
	    values_kcc_z_ion.push_back(kcc_z_ion);
	    kVectorChangeIon();
	  }


          /**
           * @brief Calculates dipolar damping parameter for the Ewald method.
           * @note According to DOI: 10.1063/1.1398588
	   * @warning 'tau_R' needs to be calculated
           */
          void calcAlphaDip() {
            double tau_R = 1.0; 
            double tau_F = timeReciDipDip.avg();
	    alpha_dip = pow(tau_R*N/tau_F,1.0/6.0)*sqrt(pc::pi)/minL;
	    values_alpha_dip.push_back(alpha_dip);
          }
          
          /**
           * @brief Calculates dipolar eciprocal space cut-off for the Ewald method.
           * @param alpha_ion_in Damping parameter
           * @param delta_in Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$)
           * @note According to DOI: 10.1063/1.1398588
           */
          void calcRcDip(double alpha_dip_in, double delta_in = 5e-5) {
	    double num = 0.5*delta_in*delta_in*N*pow(minL,3.0)/(Smu2*Smu2*pow(alpha_dip_in,4.0));
	    double dennum = -pow(delta_in,4.0)*N*N*pow(minL,6.0) /( pow(alpha_dip_in,6.0)*pow(Smu2,4.0) );
	    double denden = LambertW_solver(-3.515625*pow(delta_in,4.0)*N*N*pow(minL,6.0)/(pow(alpha_dip_in,6.0)*pow(Smu2,4.0)));
	    rc_dip = num/sqrt(dennum/denden);
	    values_rc_dip.push_back(rc_dip);
          }

          /**
           * @brief Calculates dipolar reciprocal space cut-off for the Ewald method.
           * @param alpha_ion_in Damping parameter
           * @param delta Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$), force(\f$ e^2/\AA^2 \f$) and torque(\f$ e^2/\AA \f$).
           * @note According to DOI: 10.1063/1.1398588. Here we assume statistical error in reciprocal space is negligible with regard to the systematic error ( p.6355 )
           */
          void calcKcDip(double alpha_dip_in, double delta_in = 5e-5) {
	    double temp0 = exp(-0.5*LambertW_solver((-1.125*pow(pc::pi*delta/Smu2,2.0)*N*pow(alpha_dip_in,-6.0))));
	    kc2_x_dip = 0.75*delta*L.x()*sqrt(N)*temp0/(alpha_dip_in*alpha_dip_in*Smu2);
	    kc2_y_dip = 0.75*delta*L.y()*sqrt(N)*temp0/(alpha_dip_in*alpha_dip_in*Smu2);
	    kc2_z_dip = 0.75*delta*L.z()*sqrt(N)*temp0/(alpha_dip_in*alpha_dip_in*Smu2);
            kc2_x_dip = kc_x_dip*kc_x_dip;
	    kc2_y_dip = kc_y_dip*kc_y_dip;
	    kc2_z_dip = kc_z_dip*kc_z_dip;
            kcc_x_dip = ceil(kc_x_dip);
	    kcc_y_dip = ceil(kc_y_dip);
	    kcc_z_dip = ceil(kc_z_dip);
	    values_kcc_x_dip.push_back(kcc_x_dip);
	    values_kcc_y_dip.push_back(kcc_y_dip);
	    values_kcc_z_dip.push_back(kcc_z_dip);
	    kVectorChangeDip();
          }

          /**
           * @brief Calculates optimal parameters for the Ewald method. It will only calculate a value if the user has not defined it in the input.
           * @param delta_in Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$), force(\f$ e^2/\AA^2 \f$) and torque(\f$ e^2/\AA \f$). Should not be larger than \f$ 5\cdot 10^{-5} \f$.
           * @warning Not fully implemented of tried!
           */
          void calcParameters(double delta_in=5e-5) {
            if(!user_alpha_ion)
              calcAlphaIon();
	    if(!user_rc_ion)
	      calcRcDip(alpha_ion);
            if(!user_kc_ion)
              calcKcIon(alpha_ion,delta_in);

	    if(!user_alpha_dip)
              calcAlphaDip();
	    if(!user_rc_dip)
	      calcRcDip(alpha_dip);
            if(!user_kc_dip)
              calcKcDip(alpha_dip,delta_in);
	    
	    Tbase::pairpot.first.updateAlpha(alpha_ion, alpha_dip);
	    Tbase::pairpot.first.updateRcut(rc_ion, rc_dip);
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
            o << pad(SUB,w, "Reciprocal cut-off (ion)") << "(" << kcc_x_ion << "," << kcc_y_ion << "," << kcc_z_ion << ")" << endl;
            o << pad(SUB,w, "Reciprocal cut-off (dip)") << "(" << kcc_x_dip << "," << kcc_y_dip << "," << kcc_z_dip << ")" << endl;
            o << pad(SUB,w, "Wavefunctions (ion)") << values_wavefunctions_ion.at(0) << endl;
            o << pad(SUB,w, "Wavefunctions (dip)") << values_wavefunctions_dip.at(0) << endl;
            o << pad(SUB,w, "alpha (ion)") << alpha_ion << endl;
            o << pad(SUB,w, "alpha (dip)") << alpha_dip << endl;
            o << pad(SUB,w, "Real cut-off (ion)") << rc_ion << endl;
	    o << pad(SUB,w, "Real cut-off (dip)") << rc_dip << endl;
            o << pad(SUB,w+1, epsilon_m+"(Surface)") << eps_surf << endl;
            o << pad(SUB,w, "updates") << updates << endl;
            /*
               o << pad(SUB,w+4, bracket("Self energy")) << selfEA.avg() << endl;
               o << pad(SUB,w+1, "    "+sigma) << selfEA.stdev() << endl;
               o << pad(SUB,w+4, bracket("Surf energy")) << surfEA.avg() << endl;
               o << pad(SUB,w+1, "    "+sigma) << surfEA.stdev() << endl;
               o << pad(SUB,w+4, bracket("Real energy")) << realEA.avg() << endl;
               o << pad(SUB,w+1, "    "+sigma) << realEA.stdev() << endl;
               o << pad(SUB,w+4, bracket("Reci energy")) << reciEA.avg() << endl;
               o << pad(SUB,w+1, "    "+sigma) << reciEA.stdev() << endl;
               */
            return o.str();
          }

        public:
          Average<double> selfEA, surfEA, realEA, reciEA;

          NonbondedEwald(InputMap &in) : Tbase(in) {
            Tbase::name += " (Ewald)";
            updateBoxDimensions(typename Tspace::GeometryType(in).len);

            in.cd("energy/nonbonded/ewald");
            updates = 0;
            values_alpha_ion.empty();
            values_alpha_dip.empty();
            values_kcc_x_ion.empty();
	    values_kcc_y_ion.empty();
	    values_kcc_z_ion.empty();
            values_kcc_x_dip.empty();
	    values_kcc_y_dip.empty();
	    values_kcc_z_dip.empty();
            values_wavefunctions_ion.empty();
            values_wavefunctions_dip.empty();
            lB = Tbase::pairpot.first.bjerrumLength();
            eps_surf = in("eps_surf",pc::infty);
            const_inf = 0.0;    // \varepsilon_{Surface} = \infty
            if(eps_surf < 1e10)  
              const_inf = 1.0;
            spherical_sum = in("spherical_sum", true );
	    delta = in("delta", 5e-5 );
	    
            // Initiate alpha-values
            user_alpha_ion = false;
            user_alpha_dip = false;
            if (in("alpha", -1.0  ) > 0.0 ) {
              user_alpha_ion = true;
              user_alpha_dip = true;
              alpha_ion = in("alpha", -1.0  );
              alpha2_ion = alpha_ion*alpha_ion;
              alpha_dip = alpha_ion;
              alpha2_dip = alpha2_ion;
            } else {
              if (in("alpha_ion", -1.0  ) > 0.0 ) {
                user_alpha_ion = true;
                alpha_ion = in("alpha_ion", -1.0  );
                alpha2_ion = alpha_ion*alpha_ion;
              } else {
                calcAlphaIon();
              }
              if (in("alpha_dip", -1.0  ) > 0.0 ) {
                user_alpha_dip = true;
                alpha_dip = in("alpha_dip", -1.0  );
                alpha2_dip = alpha_dip*alpha_dip;
              } else {
                calcAlphaDip();
              }
            }
            Tbase::pairpot.first.updateAlpha(alpha_ion, alpha_dip);
	    
            // Initiate real cut-off
            user_rc_ion = false;
            user_rc_dip = false;
            if (in("cutoff", -1.0  ) > 0.0 ) {
              user_rc_ion = true;
              user_rc_dip = true;
              rc_ion = in("cutoff", -1.0  );
	      if (rc_ion > minL/2.0)
		rc_ion = minL/2.0;
              rc_dip = rc_ion;
            } else {
              if (in("cutoff_ion", -1.0  ) > 0.0 ) {
                user_rc_ion = true;
                rc_ion = in("cutoff_ion", -1.0  );
		if (rc_ion > minL/2.0)
		  rc_ion = minL/2.0;
              } else {
                calcRcIon(alpha_ion);
              }
              if (in("cutoff_dip", -1.0  ) > 0.0 ) {
                user_rc_dip = true;
                rc_dip = in("cutoff_dip", -1.0  );
		if (rc_dip > minL/2.0)
		  rc_dip = minL/2.0;
              } else {
		calcRcIon(alpha_ion,delta);
              }
            }
            Tbase::pairpot.first.updateRcut(rc_ion, rc_dip);
            
            // Initiate reciprocal cut-off
            user_kc_ion = false;
            user_kc_dip = false;
            if (in("cutoffK", -1.0  ) > 0.99 ) {
              user_kc_ion = true;
              user_kc_dip = true;
	      kc_x_ion = in("cutoffK", -1.0  );
	      kc_y_ion = kc_x_ion;
	      kc_z_ion = kc_x_ion;
	      kc_x_dip = kc_x_ion;
	      kc_y_dip = kc_x_ion;
	      kc_z_dip = kc_x_ion;
            } else {
              if (in("cutoffK_ion", -1.0  ) > 0.99 ) {
                user_kc_ion = true;
		kc_x_ion = in("cutoffK_ion", -1.0  );
		kc_y_ion = kc_x_ion;
		kc_z_ion = kc_x_ion;
              } else {
		if(in("cutoffK_x_ion", -1.0  ) > 0.99 || in("cutoffK_y_ion", -1.0  ) > 0.99 || in("cutoffK_z_ion", -1.0  ) > 0.99) {
		  user_kc_ion = true;
		  kc_x_ion = in("cutoffK_x_ion", -1.0  );
		  kc_y_ion = in("cutoffK_y_ion", kc_x_ion  );
		  kc_z_ion = in("cutoffK_z_ion", kc_z_ion  );
		  
		  // Assures that kc_x_ion, kc_y_ion and kc_z_ion willhave non-negative values
		  if(kc_x_ion < 0.0) {
		    if(kc_y_ion < 0.0) {
		      kc_x_ion = kc_z_ion;
		      kc_y_ion = kc_z_ion;
		    } else {
		      kc_x_ion = kc_y_ion;
		      if(kc_z_ion < 0.0)
			kc_z_ion = kc_y_ion;
		    }
		  } else {
		    if(kc_y_ion < 0.0)
		      kc_y_ion = kc_x_ion;
		    if(kc_z_ion < 0.0)
		      kc_z_ion = kc_x_ion;
		  }
		}
              }
              if (in("cutoffK_dip", -1.0  ) > 0.99 ) {
                user_kc_dip = true;
		kc_x_dip = in("cutoffK_dip", -1.0  );
		kc_y_dip = kc_x_dip;
		kc_z_dip = kc_x_dip;
              } else {
		if(in("cutoffK_x_dip", -1.0  ) > 0.99 || in("cutoffK_y_dip", -1.0  ) > 0.99 || in("cutoffK_z_dip", -1.0  ) > 0.99) {
		  user_kc_dip = true;
		  kc_x_dip = in("cutoffK_x_dip", -1.0  );
		  kc_y_dip = in("cutoffK_y_dip", kc_x_dip  );
		  kc_z_dip = in("cutoffK_z_dip", kc_y_dip  );
		  
		  // Assures that kc_x_dip, kc_y_dip and kc_z_dip willhave non-negative values
		  if(kc_x_dip < 0.0) {
		    if(kc_y_dip < 0.0) {
		      kc_x_dip = kc_z_dip;
		      kc_y_dip = kc_z_dip;
		    } else {
		      kc_x_dip = kc_y_dip;
		      if(kc_z_dip < 0.0)
			kc_z_dip = kc_y_dip;
		    }
		  } else {
		    if(kc_y_dip < 0.0)
		      kc_y_dip = kc_x_dip;
		    if(kc_z_dip < 0.0)
		      kc_z_dip = kc_x_dip;
		  }
		}
              }
            }
            
            kc2_x_ion = kc_x_ion*kc_x_ion;
	    kc2_y_ion = kc_y_ion*kc_y_ion;
	    kc2_z_ion = kc_z_ion*kc_z_ion;
	    kc2_x_dip = kc_x_dip*kc_x_dip;
	    kc2_y_dip = kc_y_dip*kc_y_dip;
	    kc2_z_dip = kc_z_dip*kc_z_dip;
	    kcc_x_ion = ceil(kc_x_ion);
	    kcc_y_ion = ceil(kc_y_ion);
	    kcc_z_ion = ceil(kc_z_ion);
	    kcc_x_dip = ceil(kc_x_dip);
	    kcc_y_dip = ceil(kc_y_dip);
	    kcc_z_dip = ceil(kc_z_dip);       
	    
	    calcParameters(delta);
	    
	    if(user_kc_ion)
	      kVectorChangeIon();
	    if(user_kc_dip)
	      kVectorChangeDip();
          }

          double i_external(const Tpvec &p, int i) FOVERRIDE {
            Group g = Group(i,i);
            return g_external(p,g);
          }

          double g_external(const Tpvec &p, Group &g) FOVERRIDE {
	    std::chrono::steady_clock::time_point init, final;
            init = std::chrono::steady_clock::now();
            double reciEnergyIonIon = getReciEnergyIonIon(p,g);
            final = std::chrono::steady_clock::now();
            timeReciIonIon += double(std::chrono::duration_cast<std::chrono::milliseconds >(final-init).count());
            init = std::chrono::steady_clock::now();
            double reciEnergyDipDip = getReciEnergyDipDip(p,g);
            final = std::chrono::steady_clock::now();
            timeReciDipDip += double(std::chrono::duration_cast<std::chrono::milliseconds >(final-init).count());
            return lB * ( getSelfEnergy(p,g) + reciEnergyIonIon + reciEnergyDipDip );
          }

          double external(const Tpvec &p) FOVERRIDE {
            Group g(0, p.size());
            return lB * getSurfEnergy(p,g);
          }

          /**
           * @brief Set space and update k-vectors and alpha
           */
          void setSpace(Tspace &s) FOVERRIDE {
            Tbase::setSpace(s);
            N = s.p.size();
            V = spc->geo.getVolume();
            updateBoxDimensions(spc->geo.len);
	    calcParameters();
            updates++;
          }
      };

  }//namespace
}//namespace
#endif
