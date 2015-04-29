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
     */
    template<bool useIonIon=true, bool useIonDipole=true, bool useDipoleDipole=true> 
      struct EwaldReal : public Potential::Coulomb {

        typedef Potential::Coulomb Tbase;
        double alpha, alpha2, constant, rc2i;

        EwaldReal(InputMap &in, string pfx="coulomb") : Tbase(in,pfx), alpha(0), rc2i(0) {
          Tbase::name="Ewald Real";
        }
        
        void updateAlpha(double alpha_in) {
	  alpha = alpha_in;
	  alpha2 = alpha*alpha;
	  constant = 2*alpha/sqrt(pc::pi);
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
              double expK = constant*exp(-alpha2/r2i);
              T1 = (expK + r1i_d)*r2i;
              E += a.charge*r.dot(b.mu)*b.muscalar*T1;
              E -= b.charge*r.dot(a.mu)*a.muscalar*T1;
            }
            if(useDipoleDipole) {
              if(!useIonDipole) {
                double expK = constant*exp(-alpha2/r2i);
                T1 = (expK + r1i_d)*r2i;
              }
              double expK = constant*exp(-alpha2/r2i);
              double t3 = -b.mu.dot(a.mu)*T1;
              double t5 = b.mu.dot(r)*a.mu.dot(r)*(3.*r1i_d*r2i + (3.*r2i + 2.*alpha2)*expK)*r2i;
              E += -(t5 + t3)*b.muscalar*a.muscalar;
            }
#endif
            return Tbase::bjerrumLength() * E; 
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
     * `precision`       |  Precision in energy per particle in \f$ k_BT \f$ for purly ionic system.  (Default: 0.01)
     * `cutoff`          |  Real space cut-off                                                        (Default: half the box-length)                         
     * `alpha`           |  Damping parameter                                                         (Default: accorind to DOI: ...)
     * `maxK`            |  A measure of the maximum number of vectors in k-space                     (Default: according to DOI: ...)
     * `eps_surf`        |  Dielectric constant of the surroundings                                   (Default: \f$ \varepsilon_r = \infty \f$)
     * `lamell`          |  Not implemented!                                                          (Optional)
     *
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

	  Average<double> timeReci, timeReal;
	  
          int maxK, maxKz, updates, kVectorsLength, N;
          double alpha, alpha2, V, eps_surf, eps_r, Rc, rc2, rc2i, lB, const_inf, check_k2_zero, ewapre, k2max;
          Point Li, L, lamsep, lamell; // Inverse box-length, Box-length, ???, ???
          bool user_alpha, user_maxK; 
          vector<double> values_alpha, values_maxK, values_maxK_used, values_Rc;

          Eigen::MatrixXd kVectors;
          Eigen::VectorXd k2s, Aks;  // Stores values based on k-vectors in order to minimize computational effort.

          void updateKVectors() {
            kVectors.setZero();
            k2s.setZero();
            Aks.setZero();
            int wavefunctions = 0;
            
            for (int kx = -maxK; kx <= maxK; kx++) 
              for (int ky = -maxK; ky <= maxK; ky++) 
                for (int kz = -maxK; kz <= maxK; kz++) {
                  Point kv = 2*pc::pi*Point(kx*Li.x(),ky*Li.y(),kz*Li.z());
                  double k2 = kv.dot(kv);
                  if (k2 < check_k2_zero)  { // Check if k2 != 0
                    continue;
                  }
                  //if(k2 > k2max)   // For later optimization (Spherical summation)
                  //  continue;
                  kVectors.col(wavefunctions) = kv; 
                  k2s[wavefunctions] = k2;
                  Aks[wavefunctions] = exp(-k2/(4*alpha2))/k2;
                  wavefunctions++;
                }
             values_maxK_used.push_back(double(wavefunctions));
             
          }

          template<class Tpvec, class Tgroup>
            double getSelfEnergy(const Tpvec &p, const Tgroup &g) {  // const
              double Eq = 0;
              double Emu = 0;
              for (auto i : g) {
                if (useIonIon || useIonDipole)
                  Eq += p[i].charge * p[i].charge;
#ifdef DIPOLEPARTICLE
                if (useIonDipole || useDipoleDipole)
                  Emu += p[i].mu.dot(p[i].mu) * p[i].muscalar * p[i].muscalar;
#endif
              }
              selfEA += -alpha/sqrt(pc::pi) * ( Eq + 2*alpha2/3*Emu );
              return -alpha/sqrt(pc::pi) * ( Eq + 2*alpha2/3*Emu );
            }

          /**
           * @warning SOMETIMES TROUBLE WITH FEW ( 2 ) PARTICLES ! 
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
           * @brief Support-function to 'getReciEnergy'
           */
          template<class Tpvec, class Tgroup>
            double getQ2(const Tpvec &p, const Tgroup &g, const Point &kv) const {
              complex<double> Q2(0,0);
              for (auto i : g) {
                double dot = kv.dot(p[i]);
                double cosdot = cos(dot);
                double sindot = sin(dot);
                if (useIonIon || useIonDipole)
                  Q2 += p[i].charge * complex<double>(cosdot,sindot);
#ifdef DIPOLEPARTICLE
                if (useIonDipole || useDipoleDipole) {
                  Q2 += kv.dot(p[i].mu) * p[i].muscalar * complex<double>(-sindot,cosdot);
                }
#endif
              }
              double dQ2 = std::abs(Q2);
              return dQ2*dQ2;
            }

          template<class Tpvec, class Tgroup>
            double getReciEnergy(const Tpvec &p, const Tgroup &g) { // const
              if (abs(alpha2) < 1e-10)
                return 0.0;
              double E = 0.0;
              for (int i=0; i<kVectors.cols(); i++)
                E += Aks[i] * getQ2(p,g,kVectors.col(i));
              
              reciEA += (2*pc::pi/V)*E;
              return (2*pc::pi/V)*E;
            }

          void calcMaxK(double alpha2_in, Point L_in1) {
            double L_in = L_in1.norm();  // FIX!!!
            
            maxK  = std::ceil(1.0 + alpha2_in*Rc*L_in/pc::pi);
            maxKz = std::ceil(1.0 + alpha2_in*Rc*lamsep.norm()/pc::pi);   // CHECK IF NORM IS OK!!!
            
            
            kVectorsLength = (2*maxK + 1)*(2*maxK + 1)*(2*maxK + 1) - 1;  // Can be optimized !!!
            
            kVectors.resize(3, kVectorsLength); 
            k2s.resize(kVectorsLength);
            Aks.resize(kVectorsLength);
            
            values_maxK.push_back(maxK);
          }
          
          
          /**
	   * @brief Calculates damping parameter and reciprocal space cut-off for the Ewald method given a real space cut-off.
	   * @param Rc Real space cut-off. If larger than half the box-length it will reset to half the box-length!
	   * @param delta Accuracy in dipolar systems for energy (\f$ e^2/\AA \f$), force(\f$ e^2/\AA^2 \f$) and torque(\f$ e^2/\AA \f$)
	   * 
	   * Ionic systems according to DOI: 10.1080/08927029208049126  (Linse DOI: 10.1007/b136795)
	   * Dipolar systems according to DOI: 10.1063/1.1398588
	   * 
	   * @warning Not fully implemented of tried! tau_R and tau_F have to be updated when the Ewald routine changes
	   */
          void calcAlpha0(double Rc, double delta=5e-5) {
	    double tau_R = 1.0;  // Have to be updated then Ewald routine changes!
	    double tau_F = timeReci.avg();  // Have to be updated then Ewald routine changes!
	    
	    if(Rc > L.x()/2.0)
	      Rc = L.x()/2.0;
	    
	    // Ionic
	    alpha = pow(tau_R*N/tau_F,1.0/6.0)*sqrt(pc::pi)/L.x();
	    maxK = 1 + ceil(alpha*alpha*Rc*L.x()/pc::pi);
	    
	    // Dipolar
	    double A = 1.0; // Check article for true value!
	    double B = 1.0; // Check article for true value!
	    alpha = -A*sqrt(log(delta))/Rc;
	    maxK = -B*sqrt(log(delta))*alpha;
	    
	    Tbase::pairpot.first.updateAlpha(alpha);
	  }
          

          /**
           * @brief Calculates damping parameter and cut-off in reciprocal space.
	   * @param inewaldpression Precision for particle energy
           * @param rcut2 Real space cut-off
           */
          void calcAlpha(double inewaldpression, double rcut2) {
            
            double ewaldpres = inewaldpression;
            int i = 0;
            double incr = 1e-5;
            double A = 10.0;
            alpha = 0.0001;

            while( A > (1.+1e-10) || A < (1-1e-10)) {
              A=sqrt(alpha) * ewaldpres / exp(-alpha * rcut2);

              if(A>1){
                while(A > 1.){
                  alpha-= incr;
                  A = sqrt(alpha) * ewaldpres / exp(-alpha*rcut2);
                  i++;
                  if( i > 10000){
                    std::cerr << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
                    exit(1);
                    return;
                  };
                };
                if( i > 10000){
                  std::cerr << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
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
                    std::cerr << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
                    exit(1);
                    return;
                  };
                };
                if( i > 10000){
                  std::cerr << "ERROR ewaldpres = " << 1/sqrt(alpha)*exp(-alpha*rcut2) << endl;
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
            
            values_alpha.push_back(alpha);
            values_Rc.push_back(sqrt(rcut2));
            Tbase::pairpot.first.updateAlpha(alpha);
          }
          
          string _info() {
          using namespace Faunus::textio;
	  char w=25;
          std::ostringstream o;
            o << Tbase::_info();
	    o << header("Ewald summation");
	    o << pad(SUB,w, "maxK") << maxK << endl;
	    o << pad(SUB,w, "Wavefunctions") << kVectors.cols() << endl;
	    o << pad(SUB,w, "alpha") << alpha << endl;
	    o << pad(SUB,w, "Cut-off") << Rc << endl;
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
	    L = typename Tspace::GeometryType(in).len;
	    in.cd("energy/nonbonded/ewald");
            
            updates = 0;
            ewapre = in("precision",0.01);
            values_alpha.empty();
            values_maxK.empty();
            values_Rc.empty();

            lB = Tbase::pairpot.first.bjerrumLength();
            Li.x() = 1.0/L.x();
            Li.y() = 1.0/L.y();
            Li.z() = 1.0/L.z();
            lamell.x() = in("lamell_x",0.0);
            lamell.y() = in("lamell_y",0.0);
            lamell.z() = in("lamell_z",0.0);
            
            double minL = L.x();
            if(L.y() < minL)
              minL = L.y();
            if(L.z() < minL)
              minL = L.z();
            
            lamsep = lamell + L;
            Rc = in("cutoff", minL/2.0);
            if (Rc > minL/2.0)
              Rc = minL/2.0;
            rc2 = Rc*Rc;
            rc2i = 1.0/(Rc*Rc);
            Tbase::pairpot.first.rc2i=rc2i;
            eps_surf = in("eps_surf",pc::infty);
            const_inf = 0.0;    // \varepsilon_{Surface} = \infty
            if(eps_surf < 1e10)  
              const_inf = 1.0;

            alpha = 0.01;
            user_alpha = false;
            if (in("alpha", -1.0  ) != -1.0 ) {
              user_alpha = true;
              alpha = in("alpha", -1.0  );
              alpha2 = alpha*alpha;
            }
            Tbase::pairpot.first.updateAlpha(alpha);

            user_maxK = false;
            if (in("maxK", -1.0  ) != -1.0 ) {
              user_maxK = true;
              maxK = in("maxK", -1.0  );
              maxKz = maxK;
              kVectorsLength = (2*maxK + 1)*(2*maxK + 1)*(2*maxK + 1) - 1;
              kVectors.resize(3, kVectorsLength); 
              k2s.resize(kVectorsLength);
              Aks.resize(kVectorsLength);
            } else {
              calcMaxK(alpha*alpha,L);
            }
          }
          
          double i_external(const Tpvec &p, int i) FOVERRIDE {
            Group g = Group(i,i);
            return g_external(p,g);
          }

          double g_external(const Tpvec &p, Group &g) FOVERRIDE {
	    std::chrono::steady_clock::time_point init = std::chrono::steady_clock::now();
	    double reciEnergy = getReciEnergy(p,g);
	    std::chrono::steady_clock::time_point final = std::chrono::steady_clock::now();
	    timeReci += double(std::chrono::duration_cast<std::chrono::milliseconds >(final-init).count());
            return lB * ( getSelfEnergy(p,g) + reciEnergy );
          }

          double external(const Tpvec &p) FOVERRIDE {
            Group g(0, p.size());
            return lB * getSurfEnergy(p,g);
          }

          /**
           * @brief Set space and update k-vectors and alpha
           */
          void setSpace(Tspace &s) FOVERRIDE {
            updates++;
            Tbase::setSpace(s);
	    N = s.p.size();
            V = spc->geo.getVolume();
            L = spc->geo.len;
            
            double minL = L.x();
            if(L.y() < minL)
              minL = L.y();
            if(L.z() < minL)
              minL = L.z();
            if(Rc > minL) {
              Rc = minL;
              rc2 = Rc*Rc;
              rc2i = 1.0/rc2;
            }
            
            if (!user_alpha)
              calcAlpha(ewapre,rc2);
            alpha2 = alpha*alpha;
            if (!user_maxK)  {
              calcMaxK(alpha2,L);
            }
            
            double maxL = L.x();
            if(L.y() > maxL)
              maxL = L.y();
            if(L.z() > maxL)
              maxL = L.z();
            
            check_k2_zero = 0.1*(4*pc::pi*pc::pi)/(maxL*maxL);

            // NEEDED FOR LATER OPTIMIZATION !!!
            //k2max = std::ceil(4.0*pow(((pc::pi/L) + alpha2*Rc),2));  // From Martin's program
            //int k2maxtemp = std::ceil(4*pow(((pc::pi/lamsep.norm()) + alpha2*Rc),2));   // CHECK IF NORM IS OK!!!
            //if (k2maxtemp > k2max) 
            //  k2max = std::ceil(k2maxtemp);

            updateKVectors();
            Tbase::pairpot.first.updateAlpha(alpha);
          }
      };

  }//namespace
}//namespace
#endif
