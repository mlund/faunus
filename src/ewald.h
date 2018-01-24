#pragma once

#include "energy.h"
#include <complex>

namespace Faunus {
    namespace Energy {

        template<bool useIonIon=true, bool useIonDipole=false, bool useDipoleDipole=false>
            struct EwaldParameters {

                int kcc, N;
                double alpha, alpha2, rc, kc, kc2, minL, maxL, check_k2_zero;
                Point L;

                EwaldParameters() { }

                void update(Point L_in) {
                    L = L_in;
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

                /** @note Implemented from  <http://dx.doi.org/10.1080/08927029208049126> */
                void ForLaterUse() {
                    if(pc::pi/std::pow(N,1.0/3.0) > 0.25)
                        cout << "Cutoff is larger than half the box-length!" << endl;
                    rc = minL*sqrt(pc::pi)/std::pow(N,1.0/6.0);
                    alpha = pc::pi/rc;
                    alpha2 = alpha*alpha;
                    kc = alpha*minL;
                    kc2 = kc*kc;
                }
            };

        /**
         * @brief Ewald summation for electrostatic interactions
         * @warning Ewald summation does not work properly at the moment.
         */
        template<class Tspace,bool useIonIon=true, bool useIonDipole=false, bool useDipoleDipole=false>
            class NonbondedEwald : public Tbase {
                private:
                    using Tbase::spc;
                    typedef typename Tspace::Tparticle Tparticle;
                    typedef typename Tspace::Tgeometry Tgeometry;
                    typedef typename Tspace::Tpvec Tpvec;

                    EwaldParameters<useIonIon,useIonDipole,useDipoleDipole> parameters;

                    bool spherical_sum, isotropic_pbc;

                    int kVectorsInUse, kVectorsInUse_trial, N, cnt_accepted, update_frequency;

                    double V, V_trial, surfaceEnergy, surfaceEnergyTrial, reciprocalEnergy,
                           reciprocalEnergyTrial, eps_surf, const_inf, lB, update_drift; 

                    std::vector<std::complex<double>> Q_ion_tot, Q_dip_tot, Q_ion_tot_trial, Q_dip_tot_trial;

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
                                if (useIonDipole || useDipoleDipole)
                                    Emu += p[i].muscalar() * p[i].muscalar();
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
                                if (useIonDipole || useDipoleDipole)
                                    mus = mus + p[i].mu()*p[i].muscalar();
                            }

                            return const_inf * (2*pc::pi/(( 2*eps_surf + 1)*V_in))*( qrs.dot(qrs) + 2*qrs.dot(mus) +  mus.dot(mus) )*lB;
                        }

                    /**
                     * @brief Returns reciprocal space energy in kT.
                     * @param p Vector with partciles
                     * @param g Group to calculate from
                     */
                    double getReciprocalEnergy(const std::vector<std::complex<double>> &Q_ion_tot_in, const std::vector<std::complex<double>> &Q_dip_tot_in, const Eigen::VectorXd &Aks_in, double V_in) const {
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
                    void kVectorChange(
                            Eigen::MatrixXd &kVectors_in,
                            Eigen::VectorXd &Aks_in,
                            std::vector<std::complex<double>> &Q_ion_tot_in,
                            std::vector<std::complex<double>> &Q_dip_tot_in,
                            int &kVectorsInUse_in,
                            EwaldParameters<useIonIon,useIonDipole,useDipoleDipole> &parameters_in) const
                    {
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
                        int startValue = 1 - int(isotropic_pbc);

                        double factor = 1.0;
                        for (int kx = 0; kx <= parameters_in.kcc; kx++) {
                            if(kx > 0)
                                factor = 2.0;
                            double dkx2 = double(kx*kx);
                            for (int ky = -parameters_in.kcc*startValue; ky <= parameters_in.kcc; ky++) {
                                double dky2 = double(ky*ky);
                                for (int kz = -parameters_in.kcc*startValue; kz <= parameters_in.kcc; kz++) {
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
                    void updateAllComplexNumbers(const Tpvec &p,
                            std::vector<std::complex<double>> &Q_ion_tot_in,
                            std::vector<std::complex<double>> &Q_dip_tot_in,
                            const Eigen::MatrixXd &kVectors_in,
                            int kVectorsInUse_in) const
                    {
                        for (int k=0; k<kVectorsInUse_in; k++) {
                            Point kv = kVectors_in.col(k);
                            std::complex<double> Q_temp_ion(0.0,0.0);
                            std::complex<double> Q_temp_dip(0.0,0.0);
                            for (size_t i = 0; i < p.size(); i++) {
                                double dot = kv.dot(p[i]);
                                if( ( useIonIon || useIonDipole ) && !isotropic_pbc ) {
                                    Q_temp_ion += p[i].charge * std::complex<double>(cos(dot),sin(dot));
                                } else if( ( useIonIon || useIonDipole ) && isotropic_pbc ) {
                                    Q_temp_ion += p[i].charge*cos(kv.x()*p[i].x())*cos(kv.y()*p[i].y())*cos(kv.z()*p[i].z()); 
                                }
                                if(useDipoleDipole && !isotropic_pbc) {
                                    Q_temp_dip += kv.dot(p[i].mu()) * p[i].muscalar() * std::complex<double>(-sin(dot),cos(dot));
                                } else if(useDipoleDipole && isotropic_pbc) {
                                    double xx = spc->p[i].x()*kv.x();
                                    double yy = spc->p[i].y()*kv.y();
                                    double zz = spc->p[i].z()*kv.z();
                                    double cosX = cos(xx);
                                    double cosY = cos(yy);
                                    double cosZ = cos(zz);
                                    Q_temp_dip += ( sin(xx)*cosY*cosZ*spc->p[i].mu().x()*kv.x() + cosX*sin(yy)*cosZ*spc->p[i].mu().y()*kv.y() + cosX*cosY*sin(zz)*spc->p[i].mu().z()*kv.z() )*spc->p[i].muscalar();
                                }
                            }
                            Q_ion_tot_in.at(k) = Q_temp_ion;
                            Q_dip_tot_in.at(k) = Q_temp_dip;
                        }
                    }

                    /*
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
                    for(unsigned int i = 0; i < spc->p.size(); i++)
                    charge += spc->p[i].charge;
                    if(std::fabs(charge) > 1e-10 && (useIonIon || useIonDipole) )
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
    */

                public:
        MeanValue<double> selfEnergyAverage, surfaceEnergyAverage, realEnergyAverage, reciprocalEnergyAverage;

        NonbondedEwald(Tmjson &j, const string &sec="nonbonded") : Tbase(j,sec) , selfEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) ,surfaceEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) , realEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100) , reciprocalEnergyAverage(j["energy"]["nonbonded"]["avg_block"] | 100)  {
            Tbase::name += " (Ewald)";
            auto _j = j["energy"]["nonbonded"]["ewald"];
            cnt_accepted = 0;
            update_drift = 0.0;
            lB = Tbase::pairpot.first.bjerrumLength();
            eps_surf = ( _j.at("eps_surf") );
            const_inf = (eps_surf < 1) ? 0.0 : 1.0;                 // if the value is unphysical (< 1) then we set infinity as the dielectric sonatant of the surronding medium
            spherical_sum = ( _j.value("spherical_sum",true) );     // specifies if Spherical or Cubical summation should be used in reciprocal space
            update_frequency = ( _j.value("update_frequency",-1) );
            parameters.alpha = ( _j.at("alpha") );
            parameters.alpha2 = parameters.alpha*parameters.alpha;
            parameters.rc = ( _j.at("cutoff") );
            parameters.kc = ( _j.at("cutoffK") );
            parameters.kc2 = parameters.kc*parameters.kc;
            parameters.kcc = ceil(parameters.kc);
            isotropic_pbc = ( _j.value("isotropic_pbc",false) );
            Tbase::pairpot.first.updateRcut(parameters.rc);
            Tbase::pairpot.first.updateAlpha(parameters.alpha);
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
            assert(!change.empty() && "Change object is empty!");
            if (!move_accepted ) {
                // Move has been declined
                undo();
                Group g(0, spc->p.size()-1);
                selfEnergyAverage += getSelfEnergy(spc->p,g,parameters);
                surfaceEnergyAverage += surfaceEnergy;
                reciprocalEnergyAverage += reciprocalEnergy;
                //realEnergyAverage += getRealEnergy(spc->p); // Takes a lot of time
                change.clear();
                return 0.0;
            }
            // Move has been accepted
            Group g(0, spc->trial.size()-1);
            selfEnergyAverage += getSelfEnergy(spc->trial,g,parameters);
            surfaceEnergyAverage += surfaceEnergy;
            reciprocalEnergyAverage += reciprocalEnergy;
            //realEnergyAverage += getRealEnergy(spc->trial); // Takes a lot of time

            if(++cnt_accepted > update_frequency - 1) {
                double duB = getReciprocalEnergy(Q_ion_tot_trial,Q_dip_tot_trial,Aks_trial,V_trial);                        // Calulate with old vectors/matrices
                updateAllComplexNumbers(spc->trial, Q_ion_tot_trial, Q_dip_tot_trial, kVectors_trial, kVectorsInUse_trial); // Re-calculate the vectors/matrices
                double duA = getReciprocalEnergy(Q_ion_tot_trial,Q_dip_tot_trial,Aks_trial,V_trial);                        // Calulate with new vectors/matrices
                accept(); 
                cnt_accepted = 0;
                update_drift += std::fabs(duA - duB);
                change.clear();
                return (duA - duB);
            }
            accept();
            change.clear();
            return 0.0;
        }

        /** @brief Update energy function due to Change. 
         *  @warning Need to update parameters in pairpot.
         */
        double updateChange(const typename Tspace::Change &c) override {
            change = c;

            if(c.geometryChange) {
                updateAllComplexNumbers(spc->trial, Q_ion_tot_trial, Q_dip_tot_trial,kVectors_trial,kVectorsInUse_trial);
                V_trial = V + c.dV;
                parameters.update(spc->geo_trial.len);
                return 0.0;
            }

            // If the volume has not changed
            for (int k=0; k<kVectorsInUse_trial; k++) {
                std::complex<double> Q2_ion = Q_ion_tot.at(k);
                std::complex<double> Q2_dip = Q_dip_tot.at(k);
                for (auto m : change.mvGroup) {
                    for (auto i : m.second) {
                        double dotTrial = kVectors_trial.col(k).dot(spc->trial[i]);
                        double dot = kVectors.col(k).dot(spc->p[i]);
                        if ( ( useIonIon || useIonDipole ) && !isotropic_pbc ) {
                            Q2_ion += spc->trial[i].charge * std::complex<double>(cos(dotTrial),sin(dotTrial));
                            Q2_ion -= spc->p[i].charge * std::complex<double>(cos(dot),sin(dot));
                        } else if( ( useIonIon || useIonDipole ) && isotropic_pbc ) {
                            Point kv = kVectors_trial.col(k);
                            Q2_ion += spc->trial[i].charge*cos(kv.x()*spc->trial[i].x())*cos(kv.y()*spc->trial[i].y())*cos(kv.z()*spc->trial[i].z()); 
                            kv = kVectors.col(k);
                            Q2_ion -= spc->p[i].charge*cos(kv.x()*spc->p[i].x())*cos(kv.y()*spc->p[i].y())*cos(kv.z()*spc->p[i].z()); 
                        }
                        if ( ( useDipoleDipole || useIonDipole ) && !isotropic_pbc ) {
                            Q2_dip += kVectors_trial.col(k).dot(spc->trial[i].mu()) * spc->trial[i].muscalar() * std::complex<double>(-sin(dotTrial),cos(dotTrial));
                            Q2_dip -= kVectors.col(k).dot(spc->p[i].mu()) * spc->p[i].muscalar() * std::complex<double>(-sin(dot),cos(dot));
                        } else if ( ( useDipoleDipole || useIonDipole ) && isotropic_pbc ) {
                            Point kv = kVectors_trial.col(k);
                            double xx = spc->trial[i].x()*kv.x();
                            double yy = spc->trial[i].y()*kv.y();
                            double zz = spc->trial[i].z()*kv.z();
                            double cosX = cos(xx);
                            double cosY = cos(yy);
                            double cosZ = cos(zz);
                            Q2_dip += ( sin(xx)*cosY*cosZ*spc->trial[i].mu().x()*kv.x() + cosX*sin(yy)*cosZ*spc->trial[i].mu().y()*kv.y() + cosX*cosY*sin(zz)*spc->trial[i].mu().z()*kv.z() )*spc->trial[i].muscalar();
                            kv = kVectors.col(k);
                            xx = spc->p[i].x()*kv.x();
                            yy = spc->p[i].y()*kv.y();
                            zz = spc->p[i].z()*kv.z();
                            cosX = cos(xx);
                            cosY = cos(yy);
                            cosZ = cos(zz);
                            Q2_dip -= ( sin(xx)*cosY*cosZ*spc->p[i].mu().x()*kv.x() + cosX*sin(yy)*cosZ*spc->p[i].mu().y()*kv.y() + cosX*cosY*sin(zz)*spc->p[i].mu().z()*kv.z() )*spc->p[i].muscalar();
                        }
                    }
                }
                Q_ion_tot_trial.at(k) = Q2_ion;
                Q_dip_tot_trial.at(k) = Q2_dip;
            }
            return 0.0;
        }

        double i_external(const Tpvec &p, int i) override {
            Group g = Group(i,i);
            return g_external(p,g);
        }

        double g_external(const Tpvec &p, Group &g) override {
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
            parameters.update(g.len);
            if(Tbase::isGeometryTrial(g)) {
                V_trial = g.getVolume();
                kVectorChange(kVectors_trial,Aks_trial,Q_ion_tot_trial,Q_dip_tot_trial,kVectorsInUse_trial,parameters);
                updateAllComplexNumbers(spc->trial,Q_ion_tot_trial,Q_dip_tot_trial,kVectors_trial,kVectorsInUse_trial);
            } else {
                V = g.getVolume();
                kVectorChange(kVectors,Aks,Q_ion_tot,Q_dip_tot,kVectorsInUse,parameters);
                updateAllComplexNumbers(spc->p,Q_ion_tot,Q_dip_tot,kVectors,kVectorsInUse);
            }
        }

        /**
         * @brief Set space and updates parameters (if not set by user)
         */
        void setSpace(Tspace &s) override {
            Tbase::setSpace(s);
            N = s.p.size();
            Group g(0, N-1);
            surfaceEnergy = getSurfaceEnergy(s.p,g,V);
            reciprocalEnergy = getReciprocalEnergy(Q_ion_tot,Q_dip_tot,Aks,V);
            undo(); // initialization of trial-entities
            //change.clear();
        }
            };

    }//namespace
}//namespace
