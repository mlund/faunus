#pragma once

#include "energy.h"
#include <complex>

namespace Faunus {
    namespace Energy {
        /**
         * This holds Ewald setup and must *not* depend on particle type, nor depend on Space
         */
        struct EwaldData {
            typedef std::complex<double> Tcomplex;
            typedef std::vector<Tcomplex> Tvecx;
            Tvecx Qion;
            Tvecx Qdip;
            Eigen::MatrixXd kVectors; // Matrices with k-vectors
            Eigen::VectorXd Aks;      // values based on k-vectors to minimize computational effort (Eq.24,DOI:10.1063/1.481216)
            double alpha, alpha2, rc, kc, kc2, minL, maxL, check_k2_zero, lB;
            bool ionion;
            bool iondipole;
            bool dipoledipole;
            bool spherical_sum;
            bool ipbc;
            int kcc;
            int kVectorsInUse;
            int N;   //!< is this needed?
            Point L; //!< Box dimensions

            void update(const Point &box) {
                L = box;
                int kVectorsLength = (2*kcc + 1)*(2*kcc + 1)*(2*kcc + 1) - 1;
                if (kVectorsLength == 0) {
                    kVectors.resize(3, 1); 
                    Aks.resize(1);
                    kVectors.col(0) = Point(1,0,0); // Just so it is not the zero-vector
                    Aks[0] = 0.0;
                    Qion.resize(1);
                    Qdip.resize(1);
                } else {
                    kVectors.resize(3, kVectorsLength); 
                    Aks.resize(kVectorsLength);
                    kVectors.setZero();
                    Aks.setZero();
                    int startValue = 1 - int(ipbc);
                    double factor = 1;
                    for (int kx = 0; kx <= kcc; kx++) {
                        if(kx > 0)
                            factor = 2;
                        double dkx2 = double(kx*kx);
                        for (int ky = -kcc*startValue; ky <= kcc; ky++) {
                            double dky2 = double(ky*ky);
                            for (int kz = -kcc*startValue; kz <= kcc; kz++) {
                                double dkz2 = double(kz*kz);
                                Point kv = 2*pc::pi*Point(kx/L.x(),ky/L.y(),kz/L.z());
                                double k2 = kv.dot(kv);
                                if (k2 < check_k2_zero) // Check if k2 != 0
                                    continue;
                                if(spherical_sum)
                                    if( (dkx2/kc2) + (dky2/kc2) + (dkz2/kc2) > 1.0)
                                        continue;
                                kVectors.col(kVectorsInUse) = kv; 
                                Aks[kVectorsInUse] = factor*std::exp(-k2/(4.0*alpha2))/k2;
                                kVectorsInUse++;
                            }
                        }
                    }
                    Qion.resize(kVectorsInUse);
                    Qdip.resize(kVectorsInUse);
                }
            }

        };

        void from_json(const json &j, EwaldData &d) {
            d.alpha = ( j.at("alpha") );
            d.alpha2 = d.alpha*d.alpha;
            d.rc = ( j.at("cutoff") );
            d.kc = ( j.at("cutoffK") );
            d.kc2 = d.kc*d.kc;
            d.kcc = std::ceil(d.kc);
            d.ipbc = j.value("ipbc", false);
            d.spherical_sum = j.value("spherical_sum", true);
            d.lB = pc::lB( j.value("epsr", 1.0) );
            //eps_surf = ( j.at("eps_surf") );
            //const_inf = (eps_surf < 1) ? 0.0 : 1.0;  // if unphysical (< 1), set infinity as the dielectric constant of the surronding medium
        }

        template<class Tspace>
            struct PolicyIonIon {
                Tspace &spc;
                PolicyIonIon(Tspace &spc) : spc(spc) {}
                void updateComplex(EwaldData &data) const {
                    for (int k=0; k<data.kVectorsInUse; k++) {
                        const Point& kv = data.kVectors.col(k);
                        EwaldData::Tcomplex Q(0,0);
                        for (auto &i : spc.p) {
                            double dot = kv.dot(i->pos);
                            Q += i->charge * EwaldData::Tcomplex(std::cos(dot),std::sin(dot));
                        }
                        data.Qion.at(k) = Q;
                    }
                }
                double selfEnergy(EwaldData &data) {
                    double E = 0;
                    for (auto& i : spc.p)
                        E += i.charge * i.charge;
                    return -data.alpha*E / std::sqrt(pc::pi) * data.lB;
                }

                double surfaceEnergy(EwaldData &data) {
                    return 0;
                }

                double reciprocalEnergy(EwaldData &data) {
                    double E = 0;
                    for (size_t k=0; k<data.Qion.size(); k++)
                        E += data.Aks[k] * std::norm( data.Qion[k]);
                    return 2 * pc::pi / spc.geo.getVolume() * E * data.lB;
                }
            }; // recipe or policies for ion-ion ewald

        /**
         * @brief Ewald summation for electrostatic interactions
         * @warning Ewald summation does not work properly at the moment.
         */
        template<class Tspace, class Policy=PolicyIonIon<Tspace>>
            class Ewald : public Energybase {
                private:
                    Tspace& spc;
                    EwaldData data;
                    Policy policy;

                public:
                    Ewald(const json &j, Tspace &spc) : spc(spc) {
                        name = "ewald";
                        data = j;
                        data.update( spc.geo.getLength() );
                        policy.updateComplex(data); // brute force. todo: be selective
                    }

                    double energy(Change &change) override {
                        double u=0;
                        if (!change.empty()) {
                            if (key==NEW) { // we belong to a trial state
                                policy.updateComplex(data); // brute force. todo: be selective
                            }
                            u = policy.selfEnergy(data) + policy.surfaceEnergy(data) + policy.reciprocalEnergy(data); 
                        }
                        return u;
                    }

                    void sync(Energybase *basePtr, Change &change) override {
                        auto other = dynamic_cast<decltype(this)>(basePtr);
                        assert(other);
                        if (other->key==OLD) {} // do specific things if move rejected
                        if (other->key==NEW) {} // do specific things if move accepted

                        data = other->data; // copy everything!

                    } //!< Called after a move is rejected/accepted

            };

    }//namespace
}//namespace
