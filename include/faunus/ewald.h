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
        
	 /**
	  * @param M2V Input \f$ \frac{<M^2>}{9V\epsilon_0k_BT} \f$
	  * @brief Returns dielectric constant using Tinfoil conditions.
	  * \f$ 1 + \frac{<M^2>}{3V\epsilon_0k_BT} \f$
	  */
          double dielectric_constant(double M2V) const override { 
	    return (1.0 + 3.0*M2V); 
	  }

        template<class Tparticle>
          double operator() (const Tparticle &a, const Tparticle &b, const Point &r) const {
            if(!useIonIon && !useIonDipole && !useDipoleDipole)
              return 0.0;

            double r2i = 1.0/r.squaredNorm();
            if (r2i < rc2i)
              return 0.0;

            double E = 0.0;
            double r1i = sqrt(r2i);
            double r1i_d = erfc_x(alpha/r1i)*r1i;

            if(useIonIon)
              E += r1i_d*a.charge*b.charge;  // return Tbase::bjerrumLength() * a.charge * b.charge * (1.0/r1 - 1.0/rc + (r1 - rc)*rc2i );
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

	int kcc;
	double alpha, alpha2, rc, kc, kc2, maxL, check_k2_zero;
	Point L;
	
        EwaldParameters() { }

          void update(Point L_in) {
	    L = L_in;
            maxL = L.x();
            if(L.y() > maxL)
              maxL = L.y();
            if(L.z() > maxL)
              maxL = L.z();
            check_k2_zero = 0.1*(4*pc::pi*pc::pi)/(maxL*maxL);
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
	  typedef typename Tbase::Tgeometry Tgeometry;

	  EwaldParameters<useIonIon,useIonDipole,useDipoleDipole> parameters, parameters_trial;
          int kVectorsInUse, kVectorsInUse_trial, N, cnt_accepted, update_frequency;
          double V, V_trial, surfaceEnergy, surfaceEnergyTrial, reciprocalEnergy, reciprocalEnergyTrial, eps_surf, const_inf, lB, update_drift; 
          bool spherical_sum;
          vector<complex<double>> Q_ion_tot, Q_dip_tot, Q_ion_tot_trial, Q_dip_tot_trial;
          typename Tspace::Change change;

          Eigen::MatrixXd kVectors, kVectors_trial;  // Matrices with k-vectors
          Eigen::VectorXd Aks, Aks_trial;  // Stores values based on k-vectors in order to minimize computational effort. (See Eq.24 in DOI: 10.1063/1.481216)
          
          /**
           * @brief Returns Ewald self energy in kT for ions and dipoles.
           * @param p Vector with partciles
           * @param g Group to calculate from
           */
          template<class Tpvec, class Tgroup>
            double getSelfEnergy(const Tpvec &p, const Tgroup &g, EwaldParameters<useIonIon,useIonDipole,useDipoleDipole> &parameters_in) const { 
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
              return (-parameters_in.alpha*( Eq + parameters_in.alpha2*(2.0/3.0)*Emu ) / sqrt(pc::pi))*lB;
            }

          /**
           * @brief Returns Ewald surface energy in kT for ions and dipoles.
           * @param p Vector with partciles
           * @param g Group to calculate from
           */
          template<class Tpvec, class Tgroup>
            double getSurfaceEnergy(const Tpvec &p, const Tgroup &g, double V_in) const { 
	      if(const_inf < 0.5)
		return 0.0;
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
           * @brief Returns reciprocal space energy in kT.
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
           * @brief Returns real space energy in kT.
           * @param p Vector with partciles
           */
          template<class Tpvec>
	    double getRealEnergy(const Tpvec &p) const {
	      double E = 0.0;
	      for(int i = 0; i < p.size()-1; i++)
		for(int j = i+1; j < p.size(); j++)
		  E += Tbase::pairpot.first(p[i],p[j],Tbase::geo->vdist(p[i],p[j]));
	      return E;
	    }

          /**
           * @brief Updates all vectors and matrices which depends on the number of k-space vectors.
           * @note Needs to be called whenever 'kcc_x', 'kcc_y' or 'kcc_z' has been updated
           */
          void kVectorChange(Eigen::MatrixXd &kVectors_in, Eigen::VectorXd &Aks_in, vector<complex<double>> &Q_ion_tot_in, vector<complex<double>> &Q_dip_tot_in, int &kVectorsInUse_in, EwaldParameters<useIonIon,useIonDipole,useDipoleDipole> &parameters_in) const {
	    int kVectorsLength = (2*parameters_in.kcc + 1)*(2*parameters_in.kcc + 1)*(2*parameters_in.kcc + 1) - 1;
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
            for (int kx = 0; kx <= parameters_in.kcc; kx++) {
              if(kx > 0)
                factor = 2.0;
              double dkx2 = double(kx*kx);
              for (int ky = -parameters_in.kcc; ky <= parameters_in.kcc; ky++) {
                double dky2 = double(ky*ky);
                for (int kz = -parameters_in.kcc; kz <= parameters_in.kcc; kz++) {
                  double dkz2 = double(kz*kz);
                  Point kv = 2*pc::pi*Point(kx/parameters_in.L.x(),ky/parameters_in.L.y(),kz/parameters_in.L.z());
                  double k2 = kv.dot(kv);
                  if (k2 < parameters_in.check_k2_zero) // Check if k2 != 0
                    continue;
                  if(spherical_sum)
                    if( (dkx2/parameters_in.kc2) + (dky2/parameters_in.kc2) + (dkz2/parameters_in.kc2) > 1.0)
                      continue;
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
           * @brief Re-calculates the vectors of complex numbers used in getQ2
           * @param p Particle vector
	   * @param Q_ion_tot_in Vector of complex numbers for ions
	   * @param Q_dip_tot_in Vector of complex numbers for dipoles
	   * @param kVectors_in k-vectors
	   * @param kVectorsInUse_in Number of k-vectors (not necessarily the same as the length of 'kVectors_in')
           */
          void updateAllComplexNumbers(const Tpvec &p, vector<complex<double>> &Q_ion_tot_in, vector<complex<double>> &Q_dip_tot_in, const Eigen::MatrixXd &kVectors_in, int kVectorsInUse_in) const {
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

          string _info() override {
	    // Estimate the real number of wave-functions used, i.e. also those who by symmetry is implicitly accounted for
	    int realKvectors = 0;
            for (int kx = -parameters.kcc; kx <= parameters.kcc; kx++) {
	      double dkx2 = double(kx*kx);
              for (int ky = -parameters.kcc; ky <= parameters.kcc; ky++) {
		double dky2 = double(ky*ky);
                for (int kz = -parameters.kcc; kz <= parameters.kcc; kz++) {
		  double dkz2 = double(kz*kz);
		  if(spherical_sum) {
                    if( (dkx2/parameters.kc2) + (dky2/parameters.kc2) + (dkz2/parameters.kc2) > 1.0)
                      continue;
		    realKvectors++;
		  } else {
		    realKvectors++;
		  }
                }
              }
            }
	    
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
	    double charge = 0.0;
	    for(int i = 0; i < spc->p.size(); i++)
	      charge += spc->p[i].charge;
	    if(fabs(charge) > 1e-10 && (useIonIon || useIonDipole) )
	      o << pad(SUB,w, "    WARNING") << "Total charge is not 0 but " << charge << endl;
	    
            o << pad(SUB,w, "Reciprocal cut-off") << "(" << parameters.kc << "," << parameters.kc << "," << parameters.kc << ")" << endl;
	    if(spherical_sum) {
	      o << pad(SUB,w, "    Spherical") << "True" << endl;
	    } else {
	      o << pad(SUB,w, "    Spherical") << "False" << endl;
	    }
	    if(parameters.kcc == 0 && parameters.kcc == 0 && parameters.kcc == 0) {
	      o << pad(SUB,w, "Wavefunctions") << 0 << endl;
	    } else {
	      o << pad(SUB,w, "Wavefunctions") << kVectorsInUse << " (" << realKvectors << ")" << endl;
	    }
            o << pad(SUB,w, "alpha") << parameters.alpha << endl;
            o << pad(SUB,w, "Real cut-off") << parameters.rc << endl;
            if(const_inf < 0.5) {
              o << pad(SUB,w+1, epsilon_m+"(Surface)") << infinity << endl;
            } else {
              o << pad(SUB,w+1, epsilon_m+"(Surface)") << eps_surf << endl;
            }
	    o << pad(SUB,w, "Drift") << update_drift << endl;
	    
            o << pad(SUB,w+4, bracket("Energies / kT")) << endl;
            o << pad(SUB,w+4, bracket("Self energy")) << selfEnergyAverage.avg() << endl;
            o << pad(SUB,w+1, "    "+sigma) << selfEnergyAverage.std() << endl;
            o << pad(SUB,w+4, bracket("Surf energy")) << surfaceEnergyAverage.avg() << endl;
            o << pad(SUB,w+1, "    "+sigma) << surfaceEnergyAverage.std() << endl;
            o << pad(SUB,w+4, bracket("Real energy")) << realEnergyAverage.avg() << endl;
            o << pad(SUB,w+1, "    "+sigma) << realEnergyAverage.std() << endl;
            o << pad(SUB,w+4, bracket("Reci energy")) << reciprocalEnergyAverage.avg() << endl;
            o << pad(SUB,w+1, "    "+sigma) << reciprocalEnergyAverage.std() << endl;
	    o << pad(SUB,w+4, bracket("Total energy")) << selfEnergyAverage.avg()+surfaceEnergyAverage.avg()+realEnergyAverage.avg()+reciprocalEnergyAverage.avg() << endl;
	    o << pad(SUB,w+1, "    "+sigma) << selfEnergyAverage.std()+surfaceEnergyAverage.std()+realEnergyAverage.std()+reciprocalEnergyAverage.std() << endl;
            
            return o.str();
          }

        public:
          MeanValue<double> selfEnergyAverage, surfaceEnergyAverage, realEnergyAverage, reciprocalEnergyAverage;

          NonbondedEwald(Tmjson &j, const string &sec="nonbonded") : Tbase(j,sec) , selfEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) ,surfaceEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) , realEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) , reciprocalEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100)  {
	    Tbase::name += " (Ewald)";
            auto _j = j["energy"]["nonbonded"]["ewald"];
            cnt_accepted = 0;
            update_drift = 0.0;
            lB = Tbase::pairpot.first.bjerrumLength();
            eps_surf = ( _j["eps_surf"] | 0.0 );
            const_inf = (eps_surf < 1) ? 0.0 : 1.0; // if the value is unphysical (< 1) then we set infinity as the dielectric sonatant of the surronding medium
            spherical_sum = _j["spherical_sum"] | true; // specifies if Spherical or Cubical summation should be used in reciprocal space
            update_frequency = _j["update_frequency"] | -1;
	    
            parameters.alpha = ( _j["alpha"] | -1.0 );
	    parameters.alpha2 = parameters.alpha*parameters.alpha;
            parameters.rc = ( _j["cutoff"] | -1.0 );
            parameters.kc = ( _j["cutoffK"] | -1.0 );
            parameters.kc2 = parameters.kc*parameters.kc;
            parameters.kcc = ceil(parameters.kc);
	    
	    parameters_trial = parameters;
            Tbase::pairpot.first.updateAlpha(parameters.alpha);
            Tbase::pairpot.first.updateRcut(parameters.rc);
	    kVectorChange(kVectors,Aks,Q_ion_tot,Q_dip_tot,kVectorsInUse,parameters);
          }
          
          /**
	   * @brief Replaces all trial-entities with the old ones
	   */
          void undo() {
	    V_trial = V;
	    kVectorsInUse_trial = kVectorsInUse;
	    surfaceEnergyTrial = surfaceEnergy;
	    reciprocalEnergyTrial = reciprocalEnergy;
	    parameters_trial = parameters;
            Tbase::pairpot.first.updateAlpha(parameters.alpha);
            Tbase::pairpot.first.updateRcut(parameters.rc);
            Q_ion_tot_trial.resize(kVectorsInUse);
            Q_dip_tot_trial.resize(kVectorsInUse);
            kVectors_trial.resize(3, kVectorsInUse); 
	    Aks_trial.resize(kVectorsInUse); 
            for (int k=0; k < kVectorsInUse; k++) {
              Q_ion_tot_trial.at(k) = Q_ion_tot.at(k);
              Q_dip_tot_trial.at(k) = Q_dip_tot.at(k);
              kVectors_trial.col(k) = kVectors.col(k);
              Aks_trial[k] = Aks[k];
            }
	  }
	  
          /**
	   * @brief Replaces all old-entities with the trial ones
	   */
          void accept() {
	    V = V_trial;
	    kVectorsInUse = kVectorsInUse_trial;
	    surfaceEnergy = surfaceEnergyTrial;
	    reciprocalEnergy = reciprocalEnergyTrial;
	    parameters = parameters_trial;
            Tbase::pairpot.first.updateAlpha(parameters.alpha);
            Tbase::pairpot.first.updateRcut(parameters.rc);
            Q_ion_tot.resize(kVectorsInUse);
            Q_dip_tot.resize(kVectorsInUse);
            kVectors.resize(3, kVectorsInUse);
	    Aks.resize(kVectorsInUse); 
            for (int k=0; k < kVectorsInUse; k++) {
              Q_ion_tot.at(k) = Q_ion_tot_trial.at(k);
              Q_dip_tot.at(k) = Q_dip_tot_trial.at(k);
              kVectors.col(k) = kVectors_trial.col(k);
              Aks[k] = Aks_trial[k];
            }
	  }

          /**
           * @brief Updates vectors of complex numbers if a move has been accepted. These vectors are implemented in order to increase computational speed. 
           * After 'update_frequency' number of non-isobaric updates the entirety of the complex vectors are recalculated in order to avoid numerical errors.
           * @param move_accepted True if a move is accepted, otherwise false
           * 
           * @note Needs to be called before particle-vectors are updated after acception/rejection.
	   * @note Assumes energy and trial-energy has been calculated consecutively
           */
          double update(bool move_accepted) override {
	    assert(!change.empty());
            if (!move_accepted ) {
	      // Move has been declined
	      undo();
	      Group g(0, spc->p.size()-1);
	      selfEnergyAverage += getSelfEnergy(spc->p,g,parameters);
	      surfaceEnergyAverage += surfaceEnergy;
	      reciprocalEnergyAverage += reciprocalEnergy;
	      //realEnergyAverage += getRealEnergy(spc->p); // Takes a lot of time
              return 0.0;
	    }
	    // Move has been accepted
	    Group g(0, spc->trial.size()-1);
	    selfEnergyAverage += getSelfEnergy(spc->trial,g,parameters_trial);
	    surfaceEnergyAverage += surfaceEnergy;
	    reciprocalEnergyAverage += reciprocalEnergy;
	    //realEnergyAverage += getRealEnergy(spc->trial); // Takes a lot of time
	    
	    if(++cnt_accepted > update_frequency - 1) {
	      double duB = getReciprocalEnergy(Q_ion_tot_trial,Q_dip_tot_trial,Aks_trial,V_trial);                        // Calulate with old vectors/matrices
	      updateAllComplexNumbers(spc->trial, Q_ion_tot_trial, Q_dip_tot_trial, kVectors_trial, kVectorsInUse_trial); // Re-calculate the vectors/matrices
	      double duA = getReciprocalEnergy(Q_ion_tot_trial,Q_dip_tot_trial,Aks_trial,V_trial);                        // Calulate with new vectors/matrices
	      accept(); 
	      cnt_accepted = 0;
	      update_drift += fabs(duA - duB);
	      return (duA - duB);
	    }
	    accept();
	    return 0.0;
          }

          /** @brief Update energy function due to Change. 
	   *  @warning Need to update parameters in pairpot.
           */
          void setChange(const typename Tspace::Change &c) override {
            change = c;

            if(c.geometryChange) {
              updateAllComplexNumbers(spc->trial, Q_ion_tot_trial, Q_dip_tot_trial,kVectors_trial,kVectorsInUse_trial);
              V_trial = V + c.dV;
	      parameters_trial.update(spc->geo_trial.len);
              return;
            }

            // If the volume has not changed
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
          
          double i_external(const Tpvec &p, int i) override {
            Group g = Group(i,i);
	    return g_external(p,g);
          }

          double g_external(const Tpvec &p, Group &g) override {
	    if (Tbase::isTrial(p))
	      return getSelfEnergy(p,g,parameters_trial);
	    return getSelfEnergy(p,g,parameters);
          }

          double external(const Tpvec &p) override {
            Group g(0, p.size()-1);
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
            return total;
          }
    
          void setGeometry(typename Tspace::GeometryType &g) override {
            Tbase::setGeometry(g);
	    if(Tbase::isGeometryTrial(g)) {
	      V_trial = g.getVolume();
	      parameters_trial.update(g.len);
	      kVectorChange(kVectors_trial,Aks_trial,Q_ion_tot_trial,Q_dip_tot_trial,kVectorsInUse_trial,parameters_trial);
	      updateAllComplexNumbers(spc->trial,Q_ion_tot_trial,Q_dip_tot_trial,kVectors_trial,kVectorsInUse_trial);
	    } else {
	      V = g.getVolume();
	      parameters.update(g.len);
	      kVectorChange(kVectors,Aks,Q_ion_tot,Q_dip_tot,kVectorsInUse,parameters);
	      updateAllComplexNumbers(spc->p,Q_ion_tot,Q_dip_tot,kVectors,kVectorsInUse);
	    }
	  }
	  
          /**
           * @brief Set space and updates parameters (if not set by user)
           */
          void setSpace(Tspace &s) override {
            Tbase::setSpace(s);
	    setGeometry(s.geo);
            N = s.p.size();
	    Group g(0, N-1);
	    surfaceEnergy = getSurfaceEnergy(s.p,g,V);
	    reciprocalEnergy = getReciprocalEnergy(Q_ion_tot,Q_dip_tot,Aks,V);
	    undo(); // initialization of trial-entities
            change.clear();
          }
      };

  }//namespace
}//namespace
#endif


extern Faunus::InputMap in;