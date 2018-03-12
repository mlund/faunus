#pragma once
#include "tabulate.h"

namespace Faunus {
    /**
     * @brief Returns ion-dipole interaction.
     * @param QBxMuA Product of ion B:s charge and dipole A:s scalar
     * @param muA Unit dipole moment vector of particel A
     * @param QAxMuB Product of ion A:s charge and dipole B:s scalar
     * @param muB Unit dipole moment vector of particel B
     * @param r Direction \f$ r_A - r_B \f$
     */
    template<class Tvec>
        double q2mu(double QBxMuA, const Tvec &muA, double QAxMuB, const Tvec &muB, const Tvec &r) {
            double r2i = 1/r.squaredNorm();  // B = sol(dip), A = ch(charge)
            double r1i = sqrt(r2i);
            double r3i = r1i*r2i;
            double W1 = QBxMuA*muA.dot(r)*r3i;
            double W2 = QAxMuB*muB.dot(-r)*r3i;
            return (W1+W2);
        }

    /**
     * @brief Returns dipole-dipole interaction
     *
     * @param muA Unit dipole moment vector of particle A
     * @param muB Unit dipole moment vector of particle B
     * @param muAxmuB Product of dipole scalars
     * @param r Direction \f$ r_A - r_B \f$
     */
    template<class Tvec>
        double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r, double a=1.0, double b=0.0) {
#ifdef FAU_APPROXMATH
            double r1i = invsqrtQuake( r.squaredNorm() );
            double r2i = r1i*r1i;
#else
            double r2i = 1/r.squaredNorm();
            double r1i = sqrt(r2i);
#endif
            double dot = muA.dot(muB);
            double T = (3*muA.dot(r)*muB.dot(r)*r2i - dot)*a + dot*b;
            return -muAxmuB*T*r1i*r2i;
        }

    /**
     * @brief Returns ion-quadrupole interaction
     * @param qA Charge of particle A
     * @param quadB Quadrupole moment of particle B
     * @param qB Charge of particle B
     * @param quadA Quadrupole moment of particle A
     * @param r Direction \f$ r_A - r_B \f$
     */
    template<class Tvec, class Tmat>
        double q2quad(double qA, const Tmat &quadB,double qB, const Tmat &quadA, const Tvec &r) {
            double r2i = 1/r.squaredNorm();
            double r1i = sqrt(r2i);
            double r3i = r1i*r2i;
            double r5i = r3i*r2i;
            double WAB = r.transpose()*quadB*r;
            WAB = 3*WAB*r5i - quadB.trace()*r3i;
            double WBA = r.transpose()*quadA*r;
            WBA = 3*WBA*r5i - quadA.trace()*r3i;
            return (qA*WAB + qB*WBA);
        }

    namespace Potential {

        /**
         * @brief Help-function for the q-potential in class CoulombGalore
         */
        inline double qPochhammerSymbol(double q, int k=1, int P=300) {
            //P = 300 gives an error of about 10^-17 for k < 4

            double value = 1.0;
            double temp = pow(q,k);
            for(int i = 0; i < P; i++) {
                value *= (1.0 - temp);
                temp *= q;
            }
            return value;
        }

    } // potential namespace
} // faunus namespace

